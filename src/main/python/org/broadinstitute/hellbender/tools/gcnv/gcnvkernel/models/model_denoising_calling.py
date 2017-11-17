import collections
import logging
import argparse
import inspect
from abc import abstractmethod
from typing import List, Tuple, Set, Dict

import numpy as np
import pymc3 as pm
import scipy.sparse as sp
import theano as th
import theano.sparse as tst
import theano.tensor as tt
from pymc3 import Normal, Deterministic, DensityDist, Lognormal, Uniform, Exponential

from ..tasks.inference_task_base import HybridInferenceParameters, GeneralizedContinuousModel
from . import commons
from .dists import HalfFlat
from .theano_hmm import TheanoForwardBackward
from .. import config, types
from ..structs.interval import Interval, GCContentAnnotation

# import pathos.multiprocessing as mp

_logger = logging.getLogger(__name__)


class DenoisingModelConfig:
    # approximation schemes for calculating posterior expectations with respect to copy number posteriors
    _q_c_expectation_modes = ['map', 'exact', 'hybrid']

    """ Configuration of the denoising model """
    def __init__(self,
                 max_bias_factors: int = 5,
                 mapping_error_rate: float = 0.01,
                 psi_t_scale: float = 0.001,
                 psi_s_scale: float = 0.0001,
                 depth_correction_tau: float = 10000.0,
                 log_mean_bias_std: float = 0.1,
                 init_ard_rel_unexplained_variance: float = 0.1,
                 num_gc_bins: int = 20,
                 gc_curve_sd: float = 1.0,
                 q_c_expectation_mode: str = 'hybrid',
                 enable_bias_factors: bool = True,
                 enable_explicit_gc_bias_modeling: bool = False,
                 disable_bias_factors_in_flat_class: bool = False):
        """ Constructor
        :param max_bias_factors: maximum number of bias factors
        :param mapping_error_rate: typical mapping error rate
        :param psi_t_scale: scale of target-specific unexplained variance
        :param psi_s_scale: scale of sample-specific unexplained variance
        :param depth_correction_tau: precision for pinning the read depth for each denoising task to its global value
        :param log_mean_bias_std: standard deviation of the log mean bias
        :param init_ard_rel_unexplained_variance: initial ARD precision relative to unexplained variance
        :param num_gc_bins: number of GC bins (if enable_explicit_gc_bias_modeling is True)
        :param gc_curve_sd: standard deviation of each knob in the GC curve
        :param q_c_expectation_mode: approximation scheme to use for calculating posterior expectations
                                     with respect to the copy number posteriors
        :param enable_bias_factors: enable bias factor discovery
        :param enable_explicit_gc_bias_modeling: enable explicit GC bias modeling
        :param disable_bias_factors_in_flat_class:
        """

        self.max_bias_factors = max_bias_factors
        self.mapping_error_rate = mapping_error_rate
        self.psi_t_scale = psi_t_scale
        self.psi_s_scale = psi_s_scale
        self.depth_correction_tau = depth_correction_tau
        self.log_mean_bias_std = log_mean_bias_std
        self.init_ard_rel_unexplained_variance = init_ard_rel_unexplained_variance
        self.num_gc_bins = num_gc_bins
        self.gc_curve_sd = gc_curve_sd
        self.q_c_expectation_mode = q_c_expectation_mode
        self.enable_bias_factors = enable_bias_factors
        self.enable_explicit_gc_bias_modeling = enable_explicit_gc_bias_modeling
        self.disable_bias_factors_in_flat_class = disable_bias_factors_in_flat_class

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, hide: Set[str] = None):
        group = args.add_argument_group(title="Coverage denoising model parameters")
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(DenoisingModelConfig.__init__).parameters

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        process_and_maybe_add("max_bias_factors",
                              type=int,
                              help="Maximum number of bias factors")

        process_and_maybe_add("mapping_error_rate",
                              type=float,
                              help="Typical mapping error rate")

        process_and_maybe_add("psi_t_scale",
                              type=float,
                              help="Typical scale of target-specific unexplained variance")

        process_and_maybe_add("psi_s_scale",
                              type=float,
                              help="Typical scale of sample-specific unexplained variance")

        process_and_maybe_add("depth_correction_tau",
                              type=float,
                              help="Precision of read depth pinning to its global value")

        process_and_maybe_add("log_mean_bias_std",
                              type=float,
                              help="Standard deviation of log mean bias")

        process_and_maybe_add("init_ard_rel_unexplained_variance",
                              type=float,
                              help="Initial value of ARD prior precision relative to the typical target-specific "
                                   "unexplained variance scale")

        process_and_maybe_add("num_gc_bins",
                              type=int,
                              help="Number of knobs on the GC curve")

        process_and_maybe_add("gc_curve_sd",
                              type=int,
                              help="Prior standard deviation of the GC curve from flat")

        process_and_maybe_add("q_c_expectation_mode",
                              type=str,
                              choices=DenoisingModelConfig._q_c_expectation_modes,
                              help="The strategy for calculating copy number posterior expectations in the denoising "
                                   "model")

        process_and_maybe_add("enable_bias_factors",
                              type=bool,
                              help="Enable discovery of bias factors")

        process_and_maybe_add("enable_explicit_gc_bias_modeling",
                              type=bool,
                              help="Enable explicit modeling of GC bias")

        process_and_maybe_add("disable_bias_factors_in_flat_class",
                              type=bool,
                              help="Disable bias factor discovery in intervals in CNV-active regions")

    @staticmethod
    def from_args_dict(args_dict: Dict):
        relevant_keys = set(inspect.getfullargspec(DenoisingModelConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return DenoisingModelConfig(**relevant_kwargs)


class CopyNumberCallingConfig:
    """ Configuration of the copy number caller """
    def __init__(self,
                 p_alt: float = 1e-6,
                 p_flat: float = 1e-3,
                 cnv_coherence_length: float = 10000.0,
                 class_coherence_length: float = 10000.0,
                 max_copy_number: int = 5,
                 initialize_to_flat_class: bool = True,
                 num_calling_processes: int = 1):
        """
        :param p_alt: prior probability of assigning an alt copy number (with respect to the contig baseline copy
                      number) for an arbitrary target
        :param p_flat: prior probability of assigning a flat copy number prior to a target
        :param cnv_coherence_length: coherence length of copy number states (in bp units)
        :param class_coherence_length: coherence length of target class (ref, flat) states (in bp units)
        :param max_copy_number: maximum allowed copy number
        :param initialize_to_flat_class: initialize to flat class on all targets
        """
        assert 0.0 <= p_alt <= 1.0
        assert 0.0 <= p_flat <= 1.0
        assert cnv_coherence_length > 0.0
        assert class_coherence_length > 0.0
        assert max_copy_number > 0
        assert max_copy_number * p_alt < 1.0
        assert num_calling_processes > 0

        self.p_alt = p_alt
        self.p_flat = p_flat
        self.cnv_coherence_length = cnv_coherence_length
        self.class_coherence_length = class_coherence_length
        self.max_copy_number = max_copy_number
        self.initialize_to_flat_class = initialize_to_flat_class
        self.num_calling_processes = num_calling_processes

        self.num_copy_number_states = max_copy_number + 1
        self.num_copy_number_classes = 2

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, hide: Set[str] = None):
        group = args.add_argument_group(title="Copy number calling parameters")
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(CopyNumberCallingConfig.__init__).parameters

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        process_and_maybe_add("p_alt",
                              type=float,
                              help="Prior probability of alt copy number with respect to contig baseline state "
                                   "in the reference copy number class")

        process_and_maybe_add("p_flat",
                              type=float,
                              help="Prior probability of using flat copy number distribution as prior")

        process_and_maybe_add("cnv_coherence_length",
                              type=float,
                              help="Coherence length of CNV events (in the units of bp)")

        process_and_maybe_add("class_coherence_length",
                              type=float,
                              help="Coherence length of copy number classes (in the units of bp)")

        process_and_maybe_add("max_copy_number",
                              type=int,
                              help="Highest considered copy number")

        process_and_maybe_add("initialize_to_flat_class",
                              type=bool,
                              help="Initialize with flat copy number prior everywhere")

        process_and_maybe_add("num_calling_processes",
                              type=int,
                              help="Number of concurrent forward-backward threads (not implemented yet)")

    @staticmethod
    def from_args_dict(args_dict: Dict):
        relevant_keys = set(inspect.getfullargspec(CopyNumberCallingConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return CopyNumberCallingConfig(**relevant_kwargs)


class PosteriorInitializer:
    """ Base class for posterior initializers """
    @staticmethod
    @abstractmethod
    def initialize_posterior(denoising_config: DenoisingModelConfig,
                             calling_config: CopyNumberCallingConfig,
                             shared_workspace: 'DenoisingCallingWorkspace') -> None:
        raise NotImplementedError


class DefaultPosteriorInitializer(PosteriorInitializer):
    """ Initialize posteriors to priors """
    @staticmethod
    def initialize_posterior(denoising_config: DenoisingModelConfig,
                             calling_config: CopyNumberCallingConfig,
                             shared_workspace: 'DenoisingCallingWorkspace'):
        # class log posterior probs
        if calling_config.initialize_to_flat_class:
            log_q_tau_tk = (-np.log(calling_config.num_copy_number_classes - 1)
                            * np.ones((shared_workspace.num_targets, calling_config.num_copy_number_classes),
                                      dtype=types.floatX))
            log_q_tau_tk[:, 0] = -np.inf
        else:
            log_q_tau_tk = np.tile(np.log(shared_workspace.class_probs_k.get_value(borrow=True)),
                                   (shared_workspace.num_targets, 1))
        shared_workspace.log_q_tau_tk = th.shared(log_q_tau_tk, name="log_q_tau_tk", borrow=config.borrow_numpy)

        # copy number log posterior probs
        log_q_c_stc = np.zeros((shared_workspace.num_samples, shared_workspace.num_targets,
                                calling_config.num_copy_number_states), dtype=types.floatX)

        # auxiliary variables
        c_map_st = np.zeros((shared_workspace.num_samples, shared_workspace.num_targets), dtype=np.int)

        log_p_alt = np.log(calling_config.p_alt)
        log_p_baseline = np.log(1.0 - calling_config.max_copy_number * calling_config.p_alt)
        for si in range(shared_workspace.num_samples):
            sample_baseline_copy_number = shared_workspace.baseline_copy_number_s[si]
            log_q_c_stc[si, :, :] = log_p_alt
            log_q_c_stc[si, :, sample_baseline_copy_number] = log_p_baseline
            c_map_st[si, :] = sample_baseline_copy_number

        shared_workspace.log_q_c_stc = th.shared(log_q_c_stc, name="log_q_c_stc", borrow=config.borrow_numpy)
        shared_workspace.c_map_st = th.shared(c_map_st, name="c_map_st", borrow=config.borrow_numpy)


class DenoisingCallingWorkspace:
    """ This class contains objects (numpy arrays, theano tensors, etc) shared between the denoising model
    and the copy number caller """
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 targets_interval_list: List[Interval],
                 n_st: np.ndarray,
                 baseline_copy_number_s: np.ndarray,
                 global_read_depth_s: np.ndarray,
                 posterior_initializer: PosteriorInitializer = DefaultPosteriorInitializer):
        """
        Note:
            - sample_names must be in the same order as the rows of n_st
            - sample names in sample_metadata_collection must match those in sample_names

        """
        assert n_st.ndim == 2, "read counts matrix must be a dim=2 ndarray with shape (num_samples, num_targets)"

        self.num_samples: int = n_st.shape[0]
        self.num_targets: int = n_st.shape[1]

        assert len(targets_interval_list) == self.num_targets,\
            "the length of the targets interval list is incompatible with the shape of the read counts matrix"
        self.targets_interval_list = targets_interval_list

        assert baseline_copy_number_s.ndim == 1, "TODO informative message"
        assert baseline_copy_number_s.size == self.num_samples, "TODO informative message"
        assert baseline_copy_number_s.dtype in types.uint_dtypes, "TODO informative message"
        assert np.max(baseline_copy_number_s) <= calling_config.max_copy_number, "TODO informative message"
        self.baseline_copy_number_s = baseline_copy_number_s.astype(types.small_uint)

        assert global_read_depth_s.ndim == 1, "TODO informative message"
        assert global_read_depth_s.size == self.num_samples, "TODO informative message"
        self.global_read_depth_s = global_read_depth_s.astype(types.floatX)

        # todo note the type change!
        # read counts array as a shared theano tensor
        self.n_st: types.TensorSharedVariable = th.shared(
            n_st.astype(types.floatX), name="n_st", borrow=config.borrow_numpy)

        # distance between subsequent targets
        self.dist_t: types.TensorSharedVariable = th.shared(
            np.asarray([self.targets_interval_list[ti + 1].distance(self.targets_interval_list[ti])
                        for ti in range(self.num_targets - 1)], dtype=types.floatX),
            borrow=config.borrow_numpy)

        # copy number values for each copy number state
        copy_number_values_c = np.arange(0, calling_config.num_copy_number_states, dtype=types.small_uint)
        self.copy_number_values_c = th.shared(copy_number_values_c, name='copy_number_values_c',
                                              borrow=config.borrow_numpy)

        # copy number log posterior and derived quantities
        #   initialized by PosteriorInitializer.initialize(),
        #   subsequently updated by HHMMClassAndCopyNumberCaller.update_copy_number_log_posterior()
        self.log_q_c_stc: types.TensorSharedVariable = None

        # auxiliary
        self.c_map_st: types.TensorSharedVariable = None

        # copy number emission log posterior
        #   updated by LogEmissionPosteriorSampler.update_log_copy_number_emission_posterior()
        log_copy_number_emission_stc = np.zeros(
            (self.num_samples, self.num_targets, calling_config.num_copy_number_states), dtype=types.floatX)
        self.log_copy_number_emission_stc: types.TensorSharedVariable = th.shared(
            log_copy_number_emission_stc, name="log_copy_number_emission_stc", borrow=config.borrow_numpy)

        # class log posterior
        #   initialized by PosteriorInitializer.initialize()
        #   subsequently updated by HHMMClassAndCopyNumberCaller.update_class_log_posterior()
        self.log_q_tau_tk: types.TensorSharedVariable = None

        # class emission log posterior
        #   updated by HHMMClassAndCopyNumberCaller.update_log_class_emission_tk()
        log_class_emission_tk = np.zeros(
            (self.num_targets, calling_config.num_copy_number_classes), dtype=types.floatX)
        self.log_class_emission_tk: types.TensorSharedVariable = th.shared(
            log_class_emission_tk, name="log_class_emission_tk", borrow=True)

        # class assignment prior probabilities
        # Note:
        #   - The first class is the "ref" class (highly biased toward the baseline copy number)
        #     The second class is a "flat" class (all copy number states are equally probable)
        class_probs_k = np.asarray([1.0 - calling_config.p_flat, calling_config.p_flat], dtype=types.floatX)
        self.class_probs_k: types.TensorSharedVariable = th.shared(
            class_probs_k, name='class_probs_k', borrow=config.borrow_numpy)

        # class Markov chain log prior
        #   initialized here and remains constant throughout
        self.log_prior_k: np.ndarray = np.log(class_probs_k)

        # class Markov chain log transition
        #   initialized here and remains constant throughout
        self.log_trans_tkk: np.ndarray = self._get_log_trans_tkk(
            self.dist_t.get_value(borrow=True),
            calling_config.class_coherence_length,
            calling_config.num_copy_number_classes,
            class_probs_k)

        # GC bias factors
        self.W_gc_tg: tst.SparseConstant = None
        if denoising_config.enable_explicit_gc_bias_modeling:
            self.W_gc_tg = self._create_sparse_gc_bin_tensor_tg(
                self.targets_interval_list, denoising_config.num_gc_bins)

        # initialize posterior
        posterior_initializer.initialize_posterior(denoising_config, calling_config, self)

    @staticmethod
    def _get_log_trans_tkk(dist_t: np.ndarray,
                           class_coherence_length: float,
                           num_copy_number_classes: int,
                           class_probs_k: np.ndarray) -> np.ndarray:
        """ Calculates the log transition probability between copy number classes """
        stay_t = np.exp(-dist_t / class_coherence_length)
        not_stay_t = np.ones_like(stay_t) - stay_t
        delta_kl = np.eye(num_copy_number_classes, dtype=types.floatX)
        trans_tkl = (not_stay_t[:, None, None] * class_probs_k[None, None, :]
                     + stay_t[:, None, None] * delta_kl[None, :, :])
        return np.log(trans_tkl)

    @staticmethod
    def _create_sparse_gc_bin_tensor_tg(targets_interval_list: List[Interval], num_gc_bins: int) -> tst.SparseConstant:
        """ Creates a sparse 2d theano tensor with shape (num_targets, gc_bin). The sparse tensor represents a
        1-hot mapping of each target to its GC bin index. The range [0, 1] is uniformly divided into num_gc_bins.
        """
        assert all([GCContentAnnotation.get_key() in interval.annotations.keys()
                    for interval in targets_interval_list]), "explicit GC bias modeling is enabled, however, " \
                                                             "some or all targets lack the GC_CONTENT annotation."

        def get_gc_bin_idx(gc_content):
            return min(int(gc_content * num_gc_bins), num_gc_bins - 1)

        num_targets = len(targets_interval_list)
        data = np.ones((num_targets,))
        indices = [get_gc_bin_idx(interval.get_annotation(GCContentAnnotation.get_key()))
                   for interval in targets_interval_list]
        indptr = np.arange(0, num_targets + 1)
        scipy_gc_matrix = sp.csr_matrix((data, indices, indptr), shape=(num_targets, num_gc_bins),
                                        dtype=types.small_uint)
        theano_gc_matrix: tst.SparseConstant = tst.as_sparse(scipy_gc_matrix)
        return theano_gc_matrix


class InitialModelParametersSupplier:
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace):
        self.denoising_model_config = denoising_model_config
        self.calling_config = calling_config
        self.shared_workspace = shared_workspace

    @abstractmethod
    def get_init_psi_t(self) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def get_init_log_mean_bias_t(self) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def get_init_ard_u(self) -> np.ndarray:
        raise NotImplementedError


class DefaultInitialModelParametersSupplier(InitialModelParametersSupplier):
    """ TODO """
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace):
        super().__init__(denoising_model_config, calling_config, shared_workspace)

    def get_init_psi_t(self) -> np.ndarray:
        return self.denoising_model_config.psi_t_scale * np.ones(
            (self.shared_workspace.num_targets,), dtype=types.floatX)

    # todo better initialization?
    def get_init_log_mean_bias_t(self) -> np.ndarray:
        return np.zeros((self.shared_workspace.num_targets,), dtype=types.floatX)

    def get_init_ard_u(self) -> np.ndarray:
        fact = self.denoising_model_config.psi_t_scale * self.denoising_model_config.init_ard_rel_unexplained_variance
        return fact * np.ones((self.denoising_model_config.max_bias_factors,), dtype=types.floatX)


class DenoisingModel(GeneralizedContinuousModel):
    """ The gCNV coverage denoising model """
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 shared_workspace: DenoisingCallingWorkspace,
                 test_value_supplier: InitialModelParametersSupplier):
        super().__init__()
        self.shared_workspace = shared_workspace
        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        eps = denoising_model_config.mapping_error_rate

        # target-specific unexplained variance
        psi_t = Exponential(name='psi_t', lam=1.0 / denoising_model_config.psi_t_scale,
                            shape=(shared_workspace.num_targets,),
                            broadcastable=(False,))
        register_as_global(psi_t)

        # sample-specific unexplained variance
        psi_s = Exponential(name='psi_s', lam=1.0 / denoising_model_config.psi_s_scale,
                            shape=(shared_workspace.num_samples,),
                            broadcastable=(False,))
        register_as_sample_specific(psi_s)

        # convert "unexplained variance" to negative binomial over-dispersion
        alpha_st = tt.inv((tt.exp(psi_t.dimshuffle('x', 0) + psi_s.dimshuffle(0, 'x')) - 1.0))

        # target-specific mean log bias
        log_mean_bias_t = Normal(name='log_mean_bias_t', mu=0.0, sd=denoising_model_config.log_mean_bias_std,
                                 shape=(shared_workspace.num_targets,),
                                 broadcastable=(False,),
                                 testval=test_value_supplier.get_init_log_mean_bias_t())
        register_as_global(log_mean_bias_t)

        # log-normal read depth centered at the global read depth
        read_depth_mu_s = (np.log(shared_workspace.global_read_depth_s)
                           - 0.5 / denoising_model_config.depth_correction_tau)
        read_depth_s = Lognormal(name='read_depth_s',
                                 mu=read_depth_mu_s,
                                 tau=denoising_model_config.depth_correction_tau,
                                 shape=(shared_workspace.num_samples,),
                                 broadcastable=(False,),
                                 testval=shared_workspace.global_read_depth_s)
        register_as_sample_specific(read_depth_s)

        # log bias modelling, starting with the log mean bias
        log_bias_st = tt.tile(log_mean_bias_t, (shared_workspace.num_samples, 1))

        if denoising_model_config.enable_bias_factors:
            # ARD prior precisions
            ard_u = HalfFlat(name='ard_u',
                             shape=(denoising_model_config.max_bias_factors,),
                             broadcastable=(False,),
                             testval=test_value_supplier.get_init_ard_u())
            register_as_global(ard_u)

            # bias factors
            W_tu = Normal(name='W_tu', mu=0.0, tau=ard_u.dimshuffle('x', 0),
                          shape=(shared_workspace.num_targets, denoising_model_config.max_bias_factors),
                          broadcastable=(False, False))
            register_as_global(W_tu)

            # sample-specific bias factor loadings
            z_su = Normal(name='z_su', mu=0.0, sd=1.0,
                          shape=(shared_workspace.num_samples, denoising_model_config.max_bias_factors),
                          broadcastable=(False, False))
            register_as_sample_specific(z_su)

            # add contribution to total log bias
            if denoising_model_config.disable_bias_factors_in_flat_class:
                prob_ref_state_t = tt.exp(shared_workspace.log_q_tau_tk[:, 0])
                log_bias_st += (prob_ref_state_t.dimshuffle('x', 0) * tt.dot(W_tu, z_su.T).T)
            else:
                log_bias_st += tt.dot(W_tu, z_su.T).T

        # GC bias
        if denoising_model_config.enable_explicit_gc_bias_modeling:
            # sample-specific GC bias factor loadings
            z_sg = Normal(name='z_sg', mu=0.0, sd=denoising_model_config.gc_curve_sd,
                          shape=(shared_workspace.num_samples, denoising_model_config.num_gc_bins),
                          broadcastable=(False, False))
            register_as_sample_specific(z_sg)

            # add contribution to total log bias
            log_bias_st += tst.dot(shared_workspace.W_gc_tg, z_sg.T).T

        # useful expressions
        bias_st = tt.exp(log_bias_st)

        mean_mapping_error_correction_s: types.TheanoVector = eps * read_depth_s

        mu_stc = ((1.0 - eps) * read_depth_s.dimshuffle(0, 'x', 'x')
                  * bias_st.dimshuffle(0, 1, 'x')
                  * shared_workspace.copy_number_values_c.dimshuffle('x', 'x', 0)
                  + mean_mapping_error_correction_s.dimshuffle(0, 'x', 'x'))

        Deterministic(name='log_copy_number_emission_stc',
                      var=commons.negative_binomial_logp(
                          mu_stc, alpha_st.dimshuffle(0, 1, 'x'), shared_workspace.n_st.dimshuffle(0, 1, 'x')))

        # n_st (observed)
        if denoising_model_config.q_c_expectation_mode == 'map':
            def _copy_number_emission_logp(_n_st):
                mu_st = ((1.0 - eps) * read_depth_s.dimshuffle(0, 'x') * bias_st
                         * shared_workspace.c_map_st + mean_mapping_error_correction_s.dimshuffle(0, 'x'))
                log_copy_number_emission_st = commons.negative_binomial_logp(
                    mu_st, alpha_st, _n_st)
                return log_copy_number_emission_st

        elif denoising_model_config.q_c_expectation_mode == 'exact':
            def _copy_number_emission_logp(_n_st):
                _log_copy_number_emission_stc = commons.negative_binomial_logp(
                    mu_stc,
                    alpha_st.dimshuffle(0, 1, 'x'),
                    _n_st.dimshuffle(0, 1, 'x'))
                log_q_c_stc = shared_workspace.log_q_c_stc
                q_c_stc = tt.exp(log_q_c_stc)
                return tt.sum(q_c_stc * (_log_copy_number_emission_stc - log_q_c_stc), axis=2)

        elif denoising_model_config.q_c_expectation_mode == 'hybrid':
            def _copy_number_emission_logp(_n_st):
                flat_class_bitmask = tt.lt(self.shared_workspace.log_q_tau_tk[:, 0], -tt.log(2))
                flat_class_indices = flat_class_bitmask.nonzero()[0]
                ref_class_indices = (1 - flat_class_bitmask).nonzero()[0]

                # for flat classes, calculate exact posterior expectation
                mu_flat_stc = ((1.0 - eps) * read_depth_s.dimshuffle(0, 'x', 'x')
                               * bias_st.dimshuffle(0, 1, 'x')[:, flat_class_indices, :]
                               * shared_workspace.copy_number_values_c.dimshuffle('x', 'x', 0)
                               + mean_mapping_error_correction_s.dimshuffle(0, 'x', 'x'))
                alpha_flat_stc = tt.inv((tt.exp(psi_t.dimshuffle('x', 0)[:, flat_class_indices]
                                                + psi_s.dimshuffle(0, 'x')) - 1.0)).dimshuffle(0, 1, 'x')
                n_flat_stc = _n_st.dimshuffle(0, 1, 'x')[:, flat_class_indices, :]
                flat_class_logp_stc = commons.negative_binomial_logp(mu_flat_stc, alpha_flat_stc, n_flat_stc)
                log_q_c_flat_stc = shared_workspace.log_q_c_stc[:, flat_class_indices, :]
                q_c_flat_stc = tt.exp(log_q_c_flat_stc)
                flat_class_logp = tt.sum(q_c_flat_stc * (flat_class_logp_stc - log_q_c_flat_stc))

                # for ref classes, use MAP copy number state
                mu_ref_st = ((1.0 - eps) * read_depth_s.dimshuffle(0, 'x') * bias_st[:, ref_class_indices]
                             * shared_workspace.c_map_st[:, ref_class_indices]
                             + mean_mapping_error_correction_s.dimshuffle(0, 'x'))
                alpha_ref_st = alpha_st[:, ref_class_indices]
                n_ref_st = _n_st[:, ref_class_indices]
                ref_class_logp = tt.sum(commons.negative_binomial_logp(mu_ref_st, alpha_ref_st, n_ref_st))

                return flat_class_logp + ref_class_logp

        else:
            raise Exception("Unknown q_c expectation mode")

        DensityDist(name='n_st_obs',
                    logp=_copy_number_emission_logp,
                    observed=shared_workspace.n_st)

    def reset(self):
        return
        # self.shared_workspace.log_copy_number_emission_stc.zero()

    # @property
    # def more_updates(self):
    #     if not self.approx_available:
    #         return None
    #     return [(self.shared_workspace.log_copy_number_emission_stc,
    #              self.shared_workspace.log_copy_number_emission_stc
    #              + self.approx.sample_node(self['log_copy_number_emission_stc'], size=1)[0, ...])]


class CopyNumberEmissionBasicSampler:
    """ Draws posterior samples from the log copy number emission probability for a given variational approximation to
    the denoising model parameters """
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel):
        self.model_config = denoising_model_config
        self.calling_config = calling_config
        self.inference_params = inference_params
        self.shared_workspace = shared_workspace
        self.denoising_model = denoising_model
        self._simultaneous_log_copy_number_emission_sampler = None

    def update_approximation(self, approx: pm.approximations.MeanField):
        self._simultaneous_log_copy_number_emission_sampler =\
            self._get_compiled_simultaneous_log_copy_number_emission_sampler(approx)

    @property
    def is_sampler_initialized(self):
        return self._simultaneous_log_copy_number_emission_sampler is not None

    def draw(self):
        assert self.is_sampler_initialized, "posterior approximation is not provided yet"
        return self._simultaneous_log_copy_number_emission_sampler()

    # todo the code duplication here can be reduced by making log_copy_number_emission_stc a deterministic
    # todo in the denoising model declaration
    @th.configparser.change_flags(compute_test_value="off")
    def _get_compiled_simultaneous_log_copy_number_emission_sampler(self, approx: pm.approximations.MeanField):
        """ For a given variational approximation, returns a compiled theano function that draws posterior samples
        from the log copy number emission """
        log_copy_number_emission_stc = commons.node_posterior_mean_symbolic(
            approx, self.denoising_model['log_copy_number_emission_stc'],
            size=self.inference_params.log_emission_samples_per_round)
        return th.function(inputs=[], outputs=log_copy_number_emission_stc)


class HHMMClassAndCopyNumberBasicCaller:
    """ This class updates copy number and class posteriors.

        class_prior_k --► (tau_1) --► (tau_2) --► (tau_3) --► ...
                             |           |           |
                             |           |           |
                             ▼           ▼           ▼
                           (c_s1) --►  (c_s2) --►  (c_s3) --► ...
                             |           |           |
                             |           |           |
                             ▼           ▼           ▼
                            n_s1        n_s2        n_s3

        We assume the variational ansatz \prod_s p(tau, c_s | n) ~ q(tau) \prod_s q(c_s)
        Accordingly, q(tau) and q(c_s) are determined by minimizing the KL divergence w.r.t. the true
        posterior. This yields the following iterative scheme:

        - Given q(tau), the (variational) copy number prior for the first state and the copy number
          transition probabilities are determined (see _get_update_copy_number_hmm_specs_compiled_function).
          Along with the given emission probabilities to sample read counts, q(c_s) is updated using the
          forward-backward algorithm for each sample (see _update_copy_number_log_posterior)

        - Given q(c_s), the emission probability of each copy number class (tau) is determined
          (see _get_update_log_class_emission_tk_theano_func). The class prior and transition probabilities
          are fixed hyperparameters. Therefore, q(tau) can be updated immediately using a single run
          of forward-backward algorithm (see _update_class_log_posterior).
    """
    CopyNumberForwardBackwardResult = collections.namedtuple(
        'CopyNumberForwardBackwardResult',
        'sample_index, new_log_posterior_tc, copy_number_update_size, log_likelihood')

    def __init__(self,
                 calling_config: CopyNumberCallingConfig,
                 inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 disable_class_update: bool,
                 temperature: types.TensorSharedVariable):
        self.calling_config = calling_config
        self.inference_params = inference_params
        self.shared_workspace = shared_workspace
        self.disable_class_update = disable_class_update
        self.temperature = temperature

        # generate sample-specific inventory of copy number priors according to their baseline copy number state
        pi_skc = np.zeros((shared_workspace.num_samples, calling_config.num_copy_number_classes,
                           calling_config.num_copy_number_states), dtype=types.floatX)
        p_baseline = 1.0 - calling_config.max_copy_number * calling_config.p_alt
        for si in range(shared_workspace.num_samples):
            # the ref class
            pi_skc[si, 0, :] = calling_config.p_alt
            pi_skc[si, 0, shared_workspace.baseline_copy_number_s[si]] = p_baseline
            # the flat class
            pi_skc[si, 1, :] = 1.0 / calling_config.num_copy_number_states
        self.pi_skc: types.TensorSharedVariable = th.shared(pi_skc, name='pi_skc', borrow=config.borrow_numpy)

        # compiled function for forward-backward updatesof copy number posterior
        self._hmm_q_copy_number = TheanoForwardBackward(None, self.inference_params.caller_admixing_rate)

        if not disable_class_update:
            # compiled function for forward-backward update of class posterior
            self._hmm_q_class = TheanoForwardBackward(shared_workspace.log_q_tau_tk,
                                                      self.inference_params.caller_admixing_rate,
                                                      resolve_nans=(calling_config.p_flat == 0))
            # compiled function for update of class log emission
            self._update_log_class_emission_tk_theano_func = self._get_update_log_class_emission_tk_theano_func()

        else:
            self._hmm_q_class = None
            self._update_log_class_emission_tk_theano_func = None

        # compiled function for variational update of copy number HMM specs
        self._get_copy_number_hmm_specs_theano_func = self._get_compiled_copy_number_hmm_specs_theano_func()

        # compiled function for update of auxiliary variables
        self._update_aux_theano_func = self._get_update_aux_theano_func()

    def call(self,
             copy_number_update_summary_statistic_reducer,
             class_update_summary_statistic_reducer) -> Tuple[np.ndarray, np.ndarray, float, float]:
        """
        todo
        :param copy_number_update_summary_statistic_reducer:
        :param class_update_summary_statistic_reducer:
        :return:
        """
        # copy number posterior update
        copy_number_update_s, copy_number_log_likelihoods_s = self._update_copy_number_log_posterior(
            copy_number_update_summary_statistic_reducer)

        if not self.disable_class_update:
            # class posterior update
            self._update_log_class_emission_tk()
            class_update, class_log_likelihood = self._update_class_log_posterior(
                class_update_summary_statistic_reducer)
        else:
            class_update = None
            class_log_likelihood = None

        # update auxiliary variables
        self._update_aux()

        return copy_number_update_s, copy_number_log_likelihoods_s, class_update, class_log_likelihood

    def _update_copy_number_log_posterior(self, copy_number_update_summary_statistic_reducer)\
            -> Tuple[np.ndarray, np.ndarray]:
        ws = self.shared_workspace
        copy_number_update_s = np.zeros((ws.num_samples,), dtype=types.floatX)
        copy_number_log_likelihoods_s = np.zeros((ws.num_samples,), dtype=types.floatX)
        num_calling_processes = self.calling_config.num_calling_processes

        def _run_single_sample_fb(_sample_index: int):
            pi_kc = self.pi_skc.get_value(borrow=True)[_sample_index, ...]
            hmm_spec = self._get_copy_number_hmm_specs_theano_func(pi_kc)
            log_prior_c = hmm_spec[0]
            log_trans_tcc = hmm_spec[1]
            prev_log_posterior_tc = ws.log_q_c_stc.get_value(borrow=True)[_sample_index, ...]
            log_copy_number_emission_tc = ws.log_copy_number_emission_stc.get_value(borrow=True)[_sample_index, ...]

            # step 2. run forward-backward and update copy number posteriors
            fb_result = self._hmm_q_copy_number.perform_forward_backward(
                self.calling_config.num_copy_number_states,
                self.temperature.get_value()[0],
                log_prior_c, log_trans_tcc,
                log_copy_number_emission_tc, prev_log_posterior_tc)
            del log_prior_c
            del log_trans_tcc
            new_log_posterior_tc = fb_result[0]
            copy_number_update_size = copy_number_update_summary_statistic_reducer(fb_result[1])
            log_likelihood = float(fb_result[2])
            return self.CopyNumberForwardBackwardResult(
                _sample_index, new_log_posterior_tc, copy_number_update_size, log_likelihood)

        def _update_log_q_c_stc_inplace(log_q_c_stc, _sample_index, new_log_posterior_tc):
            log_q_c_stc[_sample_index, :, :] = new_log_posterior_tc[:, :]
            return log_q_c_stc

        max_chunks = ws.num_samples // num_calling_processes + 1
        for chunk_index in range(max_chunks):
            begin_index = chunk_index * num_calling_processes
            end_index = min((chunk_index + 1) * num_calling_processes, ws.num_samples)
            if begin_index >= ws.num_samples:
                break
            # todo multiprocessing
            # with mp.Pool(processes=num_calling_processes) as pool:
            #     for fb_result in pool.map(_run_single_sample_fb, range(begin_index, end_index)):
            for fb_result in [_run_single_sample_fb(sample_index)
                              for sample_index in range(begin_index, end_index)]:
                    # update log posterior in the workspace
                    ws.log_q_c_stc.set_value(
                        _update_log_q_c_stc_inplace(
                            ws.log_q_c_stc.get_value(borrow=True),
                            fb_result.sample_index, fb_result.new_log_posterior_tc),
                        borrow=True)
                    # update summary stats
                    copy_number_update_s[fb_result.sample_index] = fb_result.copy_number_update_size
                    copy_number_log_likelihoods_s[fb_result.sample_index] = fb_result.log_likelihood

        return copy_number_update_s, copy_number_log_likelihoods_s

    def _update_log_class_emission_tk(self):
        self._update_log_class_emission_tk_theano_func()

    def _update_class_log_posterior(self, class_update_summary_statistic_reducer) -> Tuple[float, float]:
        """
        todo
        :param class_update_summary_statistic_reducer:
        :return:
        """
        fb_result = self._hmm_q_class.perform_forward_backward(
            self.calling_config.num_copy_number_classes,
            self.temperature.get_value()[0],
            self.shared_workspace.log_prior_k,
            self.shared_workspace.log_trans_tkk,
            self.shared_workspace.log_class_emission_tk.get_value(borrow=True),
            self.shared_workspace.log_q_tau_tk.get_value(borrow=True))
        class_update_size = class_update_summary_statistic_reducer(fb_result[0])
        log_likelihood = float(fb_result[1])
        return class_update_size, log_likelihood

    def _update_aux(self):
        self._update_aux_theano_func()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_compiled_copy_number_hmm_specs_theano_func(self):
        """ Returns a compiled function that calculates the class-averaged and probability-sum-normalized log copy
        number transition matrix and log copy number prior for the first state

        Note:
            In the following, we use "a" and "b" subscripts in the variable names to refer to the departure
            and destination states, respectively. Like before, "t" and "k" denote target and class.
        """
        # shorthands
        pi_kc = tt.matrix('pi_kc')
        dist_t = self.shared_workspace.dist_t
        log_q_tau_tk = self.shared_workspace.log_q_tau_tk
        cnv_coherence_length = self.calling_config.cnv_coherence_length
        num_copy_number_states = self.calling_config.num_copy_number_states

        # log prior probability for the first target
        log_prior_c_first_state = tt.dot(tt.log(pi_kc.T), tt.exp(log_q_tau_tk[0, :]))
        log_prior_c_first_state -= pm.logsumexp(log_prior_c_first_state)

        # log transition matrix
        stay_t = tt.exp(-dist_t / cnv_coherence_length)  # todo can be cached in the workspace
        not_stay_t = tt.ones_like(stay_t) - stay_t
        delta_ab = tt.eye(num_copy_number_states)
        # todo use logaddexp
        log_trans_tkab = tt.log(not_stay_t.dimshuffle(0, 'x', 'x', 'x') * pi_kc.dimshuffle('x', 0, 'x', 1)
                                + stay_t.dimshuffle(0, 'x', 'x', 'x') * delta_ab.dimshuffle('x', 'x', 0, 1))
        q_tau_tkab = tt.exp(log_q_tau_tk[1:, :]).dimshuffle(0, 1, 'x', 'x')
        log_trans_tab = tt.sum(q_tau_tkab * log_trans_tkab, axis=1)
        log_trans_tab -= pm.logsumexp(log_trans_tab, axis=2)

        return th.function(inputs=[pi_kc], outputs=[log_prior_c_first_state, log_trans_tab])

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_log_class_emission_tk_theano_func(self):
        """ Returns a compiled function that calculates the log class emission and updates the
        corresponding placeholder in the shared workspace

        Note:

            xi_tab ~ posterior copy number probability of two subsequent targets

            correlations are ignored, i.e. we assume:

              xi_st(a, b) \equiv q_c(c_{s,t} = a, c_{s,t+1} = b)
                         \approx q_c(c_{s,t} = a) q_c(c_{s,t+1} = b)

            if needed, xi can be calculated exactly from the forward-backward tables

        """
        # shorthands
        dist_t = self.shared_workspace.dist_t
        q_c_stc = tt.exp(self.shared_workspace.log_q_c_stc)
        pi_skc = self.pi_skc
        cnv_coherence_length = self.calling_config.cnv_coherence_length
        num_copy_number_states = self.calling_config.num_copy_number_states

        # log copy number transition matrix for each class
        stay_t = tt.exp(-dist_t / cnv_coherence_length)
        not_stay_t = tt.ones_like(stay_t) - stay_t
        delta_ab = tt.eye(num_copy_number_states)

        # calculate log class emission by reducing over samples; see below
        log_class_emission_cum_sum_tk = tt.zeros((self.shared_workspace.num_targets - 1,
                                                  self.calling_config.num_copy_number_classes), dtype=types.floatX)

        def inc_log_class_emission_tk_except_for_first_target(pi_kc, q_c_tc, cum_sum_tk):
            """
            Adds the contribution of a given sample to the log class emission
            :param pi_kc: copy number prior inventory for the sample
            :param q_c_tc: copy number posteriors for the sample
            :param cum_sum_tk: current cumulative sum of log class emission
            :return: updated cumulative sum of log class emission
            """
            # todo use logaddexp
            log_trans_tkab = tt.log(
                not_stay_t.dimshuffle(0, 'x', 'x', 'x') * pi_kc.dimshuffle('x', 0, 'x', 1)
                + stay_t.dimshuffle(0, 'x', 'x', 'x') * delta_ab.dimshuffle('x', 'x', 0, 1))
            xi_tab = q_c_tc[:-1, :].dimshuffle(0, 1, 'x') * q_c_tc[1:, :].dimshuffle(0, 'x', 1)
            current_log_class_emission_tk = tt.sum(tt.sum(
                xi_tab.dimshuffle(0, 'x', 1, 2) * log_trans_tkab, axis=-1), axis=-1)
            return cum_sum_tk + current_log_class_emission_tk

        reduce_output = th.reduce(inc_log_class_emission_tk_except_for_first_target,
                                  sequences=[pi_skc, q_c_stc],
                                  outputs_info=[log_class_emission_cum_sum_tk])
        log_class_emission_tk_except_for_first_target = reduce_output[0]

        log_class_emission_k_first = tt.sum(tt.sum(
            tt.log(pi_skc) * q_c_stc[:, 0, :].dimshuffle(0, 'x', 1), axis=0), axis=-1)

        log_class_emission_tk = tt.concatenate((log_class_emission_k_first.dimshuffle('x', 0),
                                                log_class_emission_tk_except_for_first_target))

        return th.function(inputs=[], outputs=[], updates=[
            (self.shared_workspace.log_class_emission_tk, log_class_emission_tk)])

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_aux_theano_func(self):
        c_map_st = tt.argmax(self.shared_workspace.log_q_c_stc, axis=2)

        return th.function(inputs=[], outputs=[], updates=[
            (self.shared_workspace.c_map_st, c_map_st)
        ])
