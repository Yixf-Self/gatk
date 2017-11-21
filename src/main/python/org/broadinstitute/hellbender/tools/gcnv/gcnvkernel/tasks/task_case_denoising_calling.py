import numpy as np
import pymc3 as pm
import logging
from typing import Callable
import theano as th

from .inference_task_base import Sampler, Caller, CallerUpdateSummary, HybridInferenceTask, HybridInferenceParameters
from .. import config, types
from ..models.model_denoising_calling import DenoisingModel, DenoisingModelConfig, \
    CopyNumberEmissionBasicSampler, InitialModelParametersSupplier, \
    DenoisingCallingWorkspace, CopyNumberCallingConfig, HHMMClassAndCopyNumberBasicCaller
from ..io.io_denoising_calling import DenoisingModelImporter
from ..inference.sample_specific_opt import SampleSpecificAdamax

_logger = logging.getLogger(__name__)


class HMMCopyNumberCaller(Caller):
    """ This is a wrapper around HHMMClassAndCopyNumberBasicCaller to be used in a case-sample calling task
    Update of copy number class is disabled (as it is a global part of the model and is set by the imported
    model)

    """
    def __init__(self,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 temperature: types.TensorSharedVariable):
        self.hybrid_inference_params = hybrid_inference_params
        self.copy_number_basic_caller = HHMMClassAndCopyNumberBasicCaller(
            calling_config, hybrid_inference_params, shared_workspace, True, temperature)

    def call(self) -> 'HMMCopyNumberCallerUpdateSummary':
        copy_number_update_s, copy_number_log_likelihoods_s, _, _ = self.copy_number_basic_caller.call(
            self.hybrid_inference_params.caller_summary_statistics_reducer,
            self.hybrid_inference_params.caller_summary_statistics_reducer)
        return HMMCopyNumberCallerUpdateSummary(
            copy_number_update_s, copy_number_log_likelihoods_s,
            self.hybrid_inference_params.caller_summary_statistics_reducer)


class HMMCopyNumberCallerUpdateSummary(CallerUpdateSummary):
    def __init__(self,
                 copy_number_update_s: np.ndarray,
                 copy_number_log_likelihoods_s: np.ndarray,
                 reducer: Callable[[np.ndarray], float]):
        self.copy_number_update_s = copy_number_update_s
        self.copy_number_log_likelihoods_s = copy_number_log_likelihoods_s
        self.copy_number_update_reduced = reducer(copy_number_update_s)

    def __repr__(self):
        return "d_q_c: {0:2.6f}".format(self.copy_number_update_reduced)

    def reduce_to_scalar(self) -> float:
        return self.copy_number_update_reduced


class CopyNumberEmissionSampler(Sampler):
    """ This is a wrapper around CopyNumberEmissionBasicSampler to be used in a HybridInferenceTask """
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel):
        super().__init__(hybrid_inference_params)
        self.shared_workspace = shared_workspace
        self.calling_config = calling_config
        self.copy_number_emission_basic_sampler = CopyNumberEmissionBasicSampler(
            denoising_config, calling_config, hybrid_inference_params, shared_workspace, denoising_model)

    def update_approximation(self, approx: pm.approximations.MeanField):
        self.copy_number_emission_basic_sampler.update_approximation(approx)

    def draw(self) -> np.ndarray:
        return self.copy_number_emission_basic_sampler.draw()

    def reset(self):
        self.shared_workspace.log_copy_number_emission_stc.set_value(
            np.zeros((self.shared_workspace.num_samples,
                      self.shared_workspace.num_targets,
                      self.calling_config.num_copy_number_states),
                     dtype=types.floatX), borrow=config.borrow_numpy)

    def increment(self, update):
        self.shared_workspace.log_copy_number_emission_stc.set_value(
            self.shared_workspace.log_copy_number_emission_stc.get_value(borrow=True) + update,
            borrow=True)

    def get_latest_log_emission_expectation_estimator(self) -> np.ndarray:
        return self.shared_workspace.log_copy_number_emission_stc.get_value(borrow=True)


class CaseDenoisingCallingTask(HybridInferenceTask):
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_param_supplier: InitialModelParametersSupplier,
                 input_model_path: str):
        _logger.info("Instantiating the denoising model...")
        denoising_model = DenoisingModel(
            denoising_config, shared_workspace, initial_param_supplier)

        if hybrid_inference_params.disable_sampler:
            copy_number_emission_sampler = None
        else:
            _logger.info("Instantiating the sampler...")
            copy_number_emission_sampler = CopyNumberEmissionSampler(
                hybrid_inference_params, denoising_config, calling_config, shared_workspace, denoising_model)

        if hybrid_inference_params.disable_caller:
            copy_number_caller = None
        else:
            _logger.info("Instantiating the copy number caller...")
            initial_temperature = hybrid_inference_params.initial_temperature
            self.temperature: types.TensorSharedVariable = th.shared(
                np.asarray([initial_temperature], dtype=types.floatX))
            copy_number_caller = HMMCopyNumberCaller(
                calling_config, hybrid_inference_params, shared_workspace, self.temperature)

        elbo_normalization_factor = shared_workspace.num_samples * shared_workspace.num_targets
        opt = SampleSpecificAdamax(hybrid_inference_params.learning_rate,
                                   hybrid_inference_params.adamax_beta1,
                                   hybrid_inference_params.adamax_beta2)

        super().__init__(hybrid_inference_params, denoising_model, copy_number_emission_sampler, copy_number_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising",
                         calling_task_name="CNV calling",
                         custom_optimizer=opt)

        _logger.info("Loading the model and updating the instantiated model and workspace...")
        DenoisingModelImporter(denoising_config, calling_config, shared_workspace, denoising_model,
                               self.continuous_model_approx, input_model_path)()

    def disengage(self):
        pass

