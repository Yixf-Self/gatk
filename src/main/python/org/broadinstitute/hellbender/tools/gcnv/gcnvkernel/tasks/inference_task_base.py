import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
import inspect

import logging
import time
import tqdm
import argparse
from pymc3 import Model
from typing import List, Callable, Optional, Set, Tuple, Any, Dict
from abc import abstractmethod
from ..inference.covergence_tracker import NoisyELBOConvergenceTracker
from ..inference.param_tracker import ParamTrackerConfig, ParamTracker
from ..inference.sample_specific_opt import SampleSpecificOptimizer
from ..inference.deterministic_annealing import ADVIDeterministicAnnealing
from .. import types

_logger = logging.getLogger(__name__)


class Sampler:
    def __init__(self, hybrid_inference_params: 'HybridInferenceParameters'):
        self.hybrid_inference_params = hybrid_inference_params

    @abstractmethod
    def update_approximation(self, approx: pm.approximations.MeanField):
        raise NotImplementedError

    @abstractmethod
    def draw(self) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def reset(self):
        raise NotImplementedError

    @abstractmethod
    def increment(self, update):
        raise NotImplementedError

    @abstractmethod
    def get_latest_log_emission_expectation_estimator(self) -> np.ndarray:
        raise NotImplementedError


class Caller:
    @abstractmethod
    def call(self) -> 'CallerUpdateSummary':
        raise NotImplementedError


class CallerUpdateSummary:
    @abstractmethod
    def reduce_to_scalar(self) -> float:
        raise NotImplementedError

    @abstractmethod
    def __repr__(self):  # the summary must be representable
        raise NotImplementedError


class InferenceTask:
    # Lay in a course for starbase one two warp nine point five...
    @abstractmethod
    def engage(self):
        raise NotImplementedError("Core breach imminent!")

    @abstractmethod
    def disengage(self):
        raise NotImplementedError


class GeneralizedContinuousModel(Model):
    def __init__(self):
        self.approx: Optional[pm.MeanField] = None
        self.global_var_registry: Set[str] = set()
        self.sample_specific_var_registry: Set[str] = set()
        super().__init__()

    def register_as_global(self, var):
        self._register(var, self.global_var_registry)

    def register_as_sample_specific(self, var):
        self._register(var, self.sample_specific_var_registry)

    @staticmethod
    def _register(var, registry_set):
        assert hasattr(var, 'name')
        name = var.name
        if hasattr(var, 'transformed'):
            name = var.transformed.name
        assert name not in registry_set
        registry_set.add(name)


class HybridInferenceTask(InferenceTask):
    """
    The "hybrid" inference framework is applicable to a PGM structured as follows:

        +--------------+           +----------------+
        | discrete RVs + --------► + continuous RVs |
        +--------------+           +----------------+

    Note that discrete RVs do not have any continuous parents. The inference is approximately
    performed by factorizing the true posterior into an uncorrelated product of discrete RVs (DRVs)
    and continuous RVs (CRVs):

        p(CRVs, DRVs | observed) ~ q(CRVs) q(DRVs)

    The user must supply the following components:

        (1) a pm.Model that yields the DRV-posterior-expectation of the log joint,
            E_{DRVs ~ q(DRVs)} [log_P(CRVs, DRVs, observed)]

        (2) a "sampler" that provides samples from the log emission, defined as:
            log_emission(DRVs) = E_{CRVs ~ q(CRVs)} [log_P (observed | CRVs, DRV)]

        (3) a "caller" that updates q(DRVs) given log_emission(DRV); it could be as simple as using
            the Bayes rule, or as complicated as doing iterative hierarchical HMM if DRVs are strongly
            correlated.

    The general implementation motif is:

        (a) to store q(DRVs) as a shared theano tensor such that the the model can access it,
        (b) to store log_emission(DRVs) as a shared theano tensor (or ndarray) such that the caller
            can access it, and:
        (c) let the caller directly update the shared q(CRVs).

    This class performs mean-field ADVI to obtain q(CRVs); q(DRV) is handled by the external
    "caller" and is out the scope of this class. This class, however, requires a CallerUpdateSummary
    from the Caller in order to check for convergence.
    """

    task_modes = ['advi', 'hybrid']
    temperature_tolerance = 0.05

    def __init__(self,
                 hybrid_inference_params: 'HybridInferenceParameters',
                 continuous_model: GeneralizedContinuousModel,
                 sampler: Optional[Sampler],
                 caller: Optional[Caller],
                 **kwargs):
        assert hybrid_inference_params is not None
        self.hybrid_inference_params = hybrid_inference_params

        assert continuous_model is not None
        self.continuous_model = continuous_model

        if sampler is None:
            _logger.warning("No discrete emission sampler given; skipping the sampling step")
        self.sampler = sampler

        if caller is None:
            _logger.warning("No discrete caller given; skipping the calling step")
        self.caller = caller

        if self.hybrid_inference_params.track_model_params:
            _logger.info("Instantiating the parameter tracker...")
            self.param_tracker = self._create_param_tracker()
        else:
            self.param_tracker = None

        _logger.info("Instantiating the convergence tracker...")
        self.advi_convergence_tracker = NoisyELBOConvergenceTracker(
            self.hybrid_inference_params.convergence_snr_averaging_window,
            self.hybrid_inference_params.convergence_snr_trigger_threshold,
            self.hybrid_inference_params.convergence_snr_countdown_window)

        _logger.info("Setting up DA-ADVI...")
        with self.continuous_model:
            if not hasattr(self, 'temperature'):
                initial_temperature = self.hybrid_inference_params.initial_temperature
                self.temperature: types.TensorSharedVariable = th.shared(
                    np.asarray([initial_temperature], dtype=types.floatX))
            initial_temperature = self.temperature.get_value()[0]
            if np.abs(initial_temperature - 1.0) < 1e-10:  # no annealing
                temperature_update = None
            else:
                max_thermal_advi_iterations = self.hybrid_inference_params.max_advi_iter_first_epoch
                max_thermal_advi_iterations += ((self.hybrid_inference_params.num_thermal_epochs - 1)
                                                * self.hybrid_inference_params.max_advi_iter_first_epoch)
                temperature_drop_per_iter = (initial_temperature - 1.0) / max_thermal_advi_iterations
                temperature_update = [(self.temperature,
                                       tt.maximum(1.0, self.temperature - temperature_drop_per_iter))]

            self.continuous_model_advi = ADVIDeterministicAnnealing(
                random_seed=self.hybrid_inference_params.random_seed,
                temperature=self.temperature)
            self.continuous_model_approx: pm.MeanField = self.continuous_model_advi.approx
            if 'custom_optimizer' in kwargs.keys():
                opt = kwargs['custom_optimizer']
                assert issubclass(type(opt), SampleSpecificOptimizer)
                self.continuous_model_opt = opt.get_opt(self.continuous_model, self.continuous_model_approx)
            else:
                self.continuous_model_opt = pm.adamax(
                    learning_rate=self.hybrid_inference_params.learning_rate,
                    beta1=self.hybrid_inference_params.adamax_beta1,
                    beta2=self.hybrid_inference_params.adamax_beta2)
            self.continuous_model_step_func = self.continuous_model_advi.objective.step_function(
                score=True,
                obj_optimizer=self.continuous_model_opt,
                total_grad_norm_constraint=self.hybrid_inference_params.total_grad_norm_constraint,
                obj_n_mc=self.hybrid_inference_params.obj_n_mc,
                more_updates=temperature_update)

        if self.sampler is not None:
            self.sampler.update_approximation(self.continuous_model_approx)

        if 'elbo_normalization_factor' in kwargs.keys():
            self.elbo_normalization_factor = kwargs['elbo_normalization_factor']
        else:
            self.elbo_normalization_factor = 1.0

        if 'advi_task_name' in kwargs.keys():
            self.advi_task_name = kwargs['advi_task_name']
        else:
            self.advi_task_name = "ADVI"

        if 'sampling_task_name' in kwargs.keys():
            self.sampling_task_name = kwargs['sampling_task_name']
        else:
            self.sampling_task_name = "sampling"

        if 'calling_task_name' in kwargs.keys():
            self.calling_task_name = kwargs['calling_task_name']
        else:
            self.calling_task_name = "calling_task_name"

        self._t0 = None
        self._t1 = None
        self.elbo_hist: List[float] = []
        self.rls_elbo_hist: List[float] = []
        self.snr_hist: List[float] = []
        self.i_epoch = 1
        self.i_advi = 1
        self.calling_hist: List[Tuple[int, bool, bool]] = []
        self.previous_sampling_rounds = 0
        self.latest_caller_update_summary: Optional[CallerUpdateSummary] = None

    @abstractmethod
    def disengage(self):
        raise NotImplementedError

    def engage(self):
        try:
            all_converged = False
            while self.i_epoch <= self.hybrid_inference_params.max_training_epochs:
                _logger.debug("Starting epoch {0}...".format(self.i_epoch))
                converged_continuous = self._update_continuous_posteriors()
                all_converged = converged_continuous
                if self.sampler is not None:
                    converged_sampling = self._update_log_emission_posterior_expectation()
                    all_converged = all_converged and converged_sampling
                if self.caller is not None:
                    converged_discrete = self._update_discrete_posteriors()
                    all_converged = all_converged and converged_discrete
                self.i_epoch += 1
                if all_converged and not self._premature_convergence():
                    break
            if all_converged:
                _logger.info("Inference task completed successfully and convergence achieved.")
            else:
                _logger.warning("Inference task completed successfully but convergence not achieved.")
        except KeyboardInterrupt:
            pass

    def _log_start(self, task_name: str, i_epoch: int):
        self._t0 = time.time()
        _logger.debug("Starting {0} for epoch {1}...".format(task_name, i_epoch))

    def _log_stop(self, task_name: str, i_epoch: int):
        self._t1 = time.time()
        _logger.debug('The {0} for epoch {1} successfully finished in {2:.2f}s'.format(
            task_name, i_epoch, self._t1 - self._t0))

    def _log_interrupt(self, task_name: str, i_epoch: int):
        _logger.warning('The {0} for epoch {1} was interrupted'.format(task_name, i_epoch))

    def _premature_convergence(self):
        too_few_epochs = self.i_epoch < self.hybrid_inference_params.min_training_epochs
        still_in_annealing = np.abs(self.temperature.get_value()[0] - 1) > self.temperature_tolerance
        return too_few_epochs or still_in_annealing

    def _create_param_tracker(self):
        assert all([param_name in self.continuous_model.vars or
                    param_name in self.continuous_model.deterministics
                    for param_name in self.hybrid_inference_params.param_tracker_config.param_names]),\
            "Some of the parameters chosen to be tracker are not present in the model"
        return ParamTracker(self.hybrid_inference_params.param_tracker_config)

    def _update_continuous_posteriors(self) -> bool:
        self._log_start(self.advi_task_name, self.i_epoch)
        max_advi_iters = self.hybrid_inference_params.max_advi_iter_subsequent_epochs if self.i_epoch > 1 \
            else self.hybrid_inference_params.max_advi_iter_first_epoch
        converged = False
        with tqdm.trange(max_advi_iters, desc="({0}) starting...".format(self.advi_task_name)) as progress_bar:
            try:
                for _ in progress_bar:
                    loss = self.continuous_model_step_func() / self.elbo_normalization_factor
                    self.i_advi += 1
                    try:
                        self.advi_convergence_tracker(self.continuous_model_advi.approx, loss, self.i_advi)
                    except StopIteration:
                        if not self._premature_convergence():  # suppress signal if deemed premature
                            raise StopIteration
                    snr = self.advi_convergence_tracker.snr
                    elbo_mean = self.advi_convergence_tracker.mean
                    elbo_variance = self.advi_convergence_tracker.variance
                    if snr is not None:
                        self.snr_hist.append(snr)
                    self.elbo_hist.append(-loss)
                    self.rls_elbo_hist.append(elbo_mean)
                    progress_bar.set_description("({0} epoch {1}) ELBO: {2}, SNR: {3}, T: {4:.2f}".format(
                        self.advi_task_name,
                        self.i_epoch,
                        "{0:.3f} ± {1:.3f}".format(-elbo_mean, np.sqrt(elbo_variance))
                        if elbo_mean is not None and elbo_variance is not None else "N/A",
                        "{0:.1f}".format(snr) if snr is not None else "N/A",
                        self.temperature.get_value()[0]),
                        refresh=False)
                    if self.param_tracker is not None \
                            and self.i_advi % self.hybrid_inference_params.track_model_params_every == 0:
                        self.param_tracker(self.continuous_model_advi.approx, loss, self.i_advi)

            except StopIteration:
                converged = True
                progress_bar.close()
                self._log_stop(self.advi_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.close()
                self._log_interrupt(self.advi_task_name, self.i_epoch)
                raise KeyboardInterrupt

        return converged

    def _update_log_emission_posterior_expectation(self):
        self._log_start(self.sampling_task_name, self.i_epoch)
        if self.i_epoch == 1:
            self.sampler.reset()  # clear out log emission
        lag = min(self.previous_sampling_rounds, self.hybrid_inference_params.sampler_smoothing_window)
        converged = False
        median_rel_err = np.nan
        with tqdm.trange(self.hybrid_inference_params.log_emission_sampling_rounds,
                         desc="({0} epoch {1})".format(self.sampling_task_name, self.i_epoch)) as progress_bar:
            try:
                for i_round in progress_bar:
                    update_to_estimator = self.sampler.draw()
                    latest_estimator = self.sampler.get_latest_log_emission_expectation_estimator()
                    update_to_estimator = (update_to_estimator - latest_estimator) / (i_round + 1 + lag)
                    self.sampler.increment(update_to_estimator)
                    latest_estimator = self.sampler.get_latest_log_emission_expectation_estimator()
                    median_rel_err = np.median(np.abs(update_to_estimator / latest_estimator).flatten())
                    std_rel_err = np.std(np.abs(update_to_estimator / latest_estimator).flatten())
                    del update_to_estimator
                    progress_bar.set_description("({0} epoch {1}) relative error: {2:2.4f} ± {3:2.4f}".format(
                        self.sampling_task_name, self.i_epoch, median_rel_err, std_rel_err),
                        refresh=False)
                    if median_rel_err < self.hybrid_inference_params.log_emission_sampling_median_rel_error:
                        _logger.debug('{0} converged after {1} rounds with final '
                                      'median relative error {2:.3}.'.format(self.sampling_task_name, i_round + 1,
                                                                            median_rel_err))
                        raise StopIteration

            except StopIteration:
                converged = True
                progress_bar.refresh()
                progress_bar.close()
                self._log_stop(self.sampling_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.close()
                raise KeyboardInterrupt

            finally:
                self.previous_sampling_rounds = i_round + 1
                if not converged:
                    _logger.warning('{0} did not converge (median relative error '
                                    '= {1:.3}). Increase sampling rounds (current: {2}) if this behavior persists.'
                                    .format(self.sampling_task_name, median_rel_err,
                                            self.hybrid_inference_params.log_emission_sampling_rounds))

        return converged

    def _update_discrete_posteriors(self):
        self._log_start(self.calling_task_name, self.i_epoch)
        first_call_converged = False  # if convergence is achieved on the first call (stronger)
        iters_converged = False  # if internal loop is converged (weaker, does not imply global convergence)
        with tqdm.trange(self.hybrid_inference_params.max_calling_iters,
                         desc="({0} epoch {1})".format(self.calling_task_name, self.i_epoch)) as progress_bar:
            try:
                for i_calling_iter in progress_bar:
                    caller_summary = self.caller.call()
                    self.latest_caller_update_summary = caller_summary
                    progress_bar.set_description("({0} epoch {1}) {2}".format(
                        self.calling_task_name, self.i_epoch, repr(caller_summary)), refresh=False)
                    caller_update_size_scalar = caller_summary.reduce_to_scalar()
                    if caller_update_size_scalar < self.hybrid_inference_params.caller_update_convergence_threshold:
                        iters_converged = True
                        if i_calling_iter == 0:
                            first_call_converged = True
                        raise StopIteration

            except StopIteration:
                progress_bar.refresh()
                progress_bar.close()
                self._log_stop(self.calling_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.close()
                self._log_interrupt(self.calling_task_name, self.i_epoch)
                raise KeyboardInterrupt

            finally:
                self.calling_hist.append((self.i_advi, iters_converged, first_call_converged))
                # if there is a self-consistency loop and not converged ...
                if not iters_converged and self.hybrid_inference_params.max_calling_iters > 1:
                    _logger.warning('{0} did not converge. Increase maximum calling rounds (current: {1}) '
                                    'if this behavior persists.'.format(
                        self.calling_task_name, self.hybrid_inference_params.max_calling_iters))

        return first_call_converged


class HybridInferenceParameters:
    """ Hybrid ADVI (for continuous RVs) + external calling (for discrete RVs) inference parameters """
    def __init__(self,
                 learning_rate: float = 0.2,
                 adamax_beta1: float = 0.9,
                 adamax_beta2: float = 0.99,
                 obj_n_mc: int = 1,
                 random_seed: int = 1984,
                 total_grad_norm_constraint: Optional[float] = None,
                 log_emission_samples_per_round: int = 50,
                 log_emission_sampling_median_rel_error: float = 5e-3,
                 log_emission_sampling_rounds: int = 10,
                 max_advi_iter_first_epoch: int = 100,
                 max_advi_iter_subsequent_epochs: int = 100,
                 min_training_epochs: int = 5,
                 max_training_epochs: int = 50,
                 initial_temperature: float = 2.0,
                 num_thermal_epochs: int = 20,
                 track_model_params: bool = False,
                 track_model_params_every: int = 10,
                 param_tracker_config: Optional['ParamTrackerConfig'] = None,
                 convergence_snr_averaging_window: int = 500,
                 convergence_snr_trigger_threshold: float = 0.1,
                 convergence_snr_countdown_window: int = 10,
                 max_calling_iters: int = 10,
                 caller_update_convergence_threshold: float = 1e-3,
                 caller_admixing_rate: float = 0.75,
                 sampler_smoothing_window: int = 0,
                 caller_summary_statistics_reducer: Callable[[np.ndarray], float] = np.mean,
                 disable_sampler: bool = False,
                 disable_caller: bool = False):
        self.learning_rate = learning_rate
        self.adamax_beta1 = adamax_beta1
        self.adamax_beta2 = adamax_beta2
        self.obj_n_mc = obj_n_mc
        self.random_seed = random_seed
        self.total_grad_norm_constraint = total_grad_norm_constraint
        self.log_emission_samples_per_round = log_emission_samples_per_round
        self.log_emission_sampling_median_rel_error = log_emission_sampling_median_rel_error
        self.log_emission_sampling_rounds = log_emission_sampling_rounds
        self.max_advi_iter_first_epoch = max_advi_iter_first_epoch
        self.max_advi_iter_subsequent_epochs = max_advi_iter_subsequent_epochs
        self.min_training_epochs = min_training_epochs
        self.max_training_epochs = max_training_epochs
        self.initial_temperature = initial_temperature
        self.num_thermal_epochs = num_thermal_epochs
        self.track_model_params = track_model_params
        self.track_model_params_every = track_model_params_every
        self.param_tracker_config = param_tracker_config
        self.convergence_snr_averaging_window = convergence_snr_averaging_window
        self.convergence_snr_trigger_threshold = convergence_snr_trigger_threshold
        self.convergence_snr_countdown_window = convergence_snr_countdown_window
        self.max_calling_iters = max_calling_iters
        self.caller_update_convergence_threshold = caller_update_convergence_threshold
        self.caller_admixing_rate = caller_admixing_rate
        self.sampler_smoothing_window = sampler_smoothing_window
        self.caller_summary_statistics_reducer = caller_summary_statistics_reducer
        self.disable_sampler = disable_sampler
        self.disable_caller = disable_caller

        self._assert_params()

    # todo the rest of assertions
    def _assert_params(self):
        assert self.learning_rate >= 0
        assert self.obj_n_mc >= 0
        assert self.log_emission_samples_per_round >= 1
        assert self.log_emission_sampling_rounds >= 1
        assert 0.0 < self.log_emission_sampling_median_rel_error < 1.0

        if self.track_model_params:
            assert self.param_tracker_config is not None

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, override_default: Dict[str, Any] = None, hide: Set[str] = None):
        group = args.add_argument_group(title="Inference parameters")
        if override_default is None:
            override_default = dict()
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(HybridInferenceParameters.__init__).parameters

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            if full_arg in override_default:
                kwargs['default'] = override_default[full_arg]
            else:
                kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        process_and_maybe_add("learning_rate",
                              type=float,
                              help="Adamax optimizer learning rate")

        process_and_maybe_add("adamax_beta1",
                              type=float,
                              help="Adamax first moment estimation forgetting factor")

        process_and_maybe_add("adamax_beta2",
                              type=float,
                              help="Adamax second moment estimation forgetting factor")

        process_and_maybe_add("log_emission_samples_per_round",
                              type=int,
                              help="Number of log emission posterior samples per sampling round")

        process_and_maybe_add("log_emission_sampling_median_rel_error",
                              type=float,
                              help="Maximum tolerated median relative error in log emission posterior sampling")

        process_and_maybe_add("log_emission_sampling_rounds",
                              type=int,
                              help="Maximum log emission posterior sampling rounds")

        process_and_maybe_add("max_advi_iter_first_epoch",
                              type=int,
                              help="Maximum ADVI iterations in the first epoch (before sampling and calling)")

        process_and_maybe_add("max_advi_iter_subsequent_epochs",
                              type=int,
                              help="Maximum ADVI iterations after the first epoch")

        process_and_maybe_add("min_training_epochs",
                              type=int,
                              help="Minimum number of training epochs before evaluating convergence")

        process_and_maybe_add("max_training_epochs",
                              type=int,
                              help="Maximum number of training epochs before stopping")

        process_and_maybe_add("initial_temperature",
                              type=float,
                              help="Initial temperature for deterministic annealing (must be >= 1.0)")

        process_and_maybe_add("num_thermal_epochs",
                              type=int,
                              help="Annealing duration (in the units of training epochs)")

        process_and_maybe_add("convergence_snr_averaging_window",
                              type=int,
                              help="Averaging window for calculating training SNR for evaluating convergence")

        process_and_maybe_add("convergence_snr_trigger_threshold",
                              type=float,
                              help="The SNR threshold to be reached for triggering convergence")

        process_and_maybe_add("convergence_snr_countdown_window",
                              type=int,
                              help="The number of ADVI iterations during which the SNR is required to stay below the "
                                   "set threshold for convergence")

        process_and_maybe_add("max_calling_iters",
                              type=int,
                              help="Maximum number of calling internal self-consistency iterations")

        process_and_maybe_add("caller_update_convergence_threshold",
                              type=float,
                              help="Maximum tolerated calling update size for convergence")

        process_and_maybe_add("caller_admixing_rate",
                              type=float,
                              help="Admixing ratio of new and old caller posteriors (between 0 and 1; higher means "
                                   "more of the new posterior)")

        process_and_maybe_add("disable_sampler",
                              type=bool,
                              help="Disable sampler (advanced)")

        process_and_maybe_add("disable_caller",
                              type=bool,
                              help="Disable caller (advanced)")

    @staticmethod
    def from_args_dict(args_dict: Dict):
        relevant_keys = set(inspect.getfullargspec(HybridInferenceParameters.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return HybridInferenceParameters(**relevant_kwargs)
