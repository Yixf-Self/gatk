import numpy as np
import pymc3 as pm
import logging
from typing import Callable
from ..models import commons

from .inference_task_base import Sampler, Caller, CallerUpdateSummary,\
    HybridInferenceTask, HybridInferenceParameters
from .. import config, types
from ..models.model_ploidy import PloidyModelConfig, PloidyModel,\
    PloidyWorkspace, PloidyEmissionBasicSampler, PloidyBasicCaller
from ..structs import metadata

_logger = logging.getLogger(__name__)


class PloidyCaller(Caller):
    """ This is a wrapper around PloidyBasicCaller to be used in a HybridInferenceTask """
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_workspace: PloidyWorkspace):
        self.hybrid_inference_params = hybrid_inference_params
        self.ploidy_basic_caller = PloidyBasicCaller(hybrid_inference_params, ploidy_workspace)

    def call(self) -> 'PloidyCallerUpdateSummary':
        update_norm_sj = self.ploidy_basic_caller.call()
        return PloidyCallerUpdateSummary(
            update_norm_sj, self.hybrid_inference_params.caller_summary_statistics_reducer)


class PloidyCallerUpdateSummary(CallerUpdateSummary):
    def __init__(self,
                 update_norm_sj: np.ndarray,
                 reducer: Callable[[np.ndarray], float]):
        self.scalar_update = reducer(update_norm_sj)

    def __repr__(self):
        return "ploidy update size: {0:2.6}".format(self.scalar_update)

    def reduce_to_scalar(self) -> float:
        return self.scalar_update


class PloidyEmissionSampler(Sampler):
    """ This is a wrapper around PloidyEmissionBasicSampler to be used in a HybridInferenceTask """
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_model: PloidyModel,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__(hybrid_inference_params)
        self.ploidy_workspace = ploidy_workspace
        self.ploidy_emission_basic_sampler = PloidyEmissionBasicSampler(
            ploidy_model, self.hybrid_inference_params.log_emission_samples_per_round)

    def update_approximation(self, approx: pm.approximations.MeanField):
        self.ploidy_emission_basic_sampler.update_approximation(approx)

    def draw(self) -> np.ndarray:
        return self.ploidy_emission_basic_sampler.draw()

    def reset(self):
        self.ploidy_workspace.log_ploidy_emission_sjk.set_value(
            np.zeros((self.ploidy_workspace.num_samples,
                      self.ploidy_workspace.num_contigs,
                      self.ploidy_workspace.num_ploidy_states),
                     dtype=types.floatX), borrow=config.borrow_numpy)

    def increment(self, update):
        self.ploidy_workspace.log_ploidy_emission_sjk.set_value(
            self.ploidy_workspace.log_ploidy_emission_sjk.get_value(borrow=True) + update)

    def get_latest_log_emission_expectation_estimator(self) -> np.ndarray:
        return self.ploidy_workspace.log_ploidy_emission_sjk.get_value(borrow=True)


class PloidyInferenceTask(HybridInferenceTask):
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        _logger.info("Instantiating the contig-level coverage model...")
        ploidy_model = PloidyModel(ploidy_config, ploidy_workspace)

        _logger.info("Instantiating the ploidy emission sampler...")
        ploidy_emission_sampler = PloidyEmissionSampler(hybrid_inference_params, ploidy_model, ploidy_workspace)

        _logger.info("Instantiating the ploidy caller...")
        ploidy_caller = PloidyCaller(hybrid_inference_params, ploidy_workspace)

        elbo_normalization_factor = ploidy_workspace.num_samples * ploidy_workspace.num_contigs
        super().__init__(hybrid_inference_params, ploidy_model, ploidy_emission_sampler, ploidy_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising",
                         calling_task_name="ploidy calling")

        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

    # todo warn if ploidy genotyping quality is low
    # todo warn if ploidy genotyping is incompatible with a given list of sex genotypes
    def disengage(self):
        log_q_ploidy_sjk = self.ploidy_workspace.log_q_ploidy_sjk.get_value(borrow=True)

        for si, sample_name in enumerate(self.ploidy_workspace.sample_names):
            ploidy_j = np.zeros((self.ploidy_workspace.num_contigs,), dtype=types.small_uint)
            ploidy_genotyping_quality_j = np.zeros((self.ploidy_workspace.num_contigs,), dtype=types.floatX)
            for j in range(self.ploidy_workspace.num_contigs):
                ploidy_j[j], ploidy_genotyping_quality_j[j] = commons.perform_genotyping(log_q_ploidy_sjk[si, j, :])

            # generate sample ploidy metadata
            sample_ploidy_metadata = metadata.SamplePloidyMetadata(
                sample_name, ploidy_j, ploidy_genotyping_quality_j, self.ploidy_workspace.targets_metadata.contig_list)

            # generate sample read depth metadata
            sample_read_depth_metadata = metadata.SampleReadDepthMetadata.generate_sample_read_depth_metadata(
                self.ploidy_workspace.sample_metadata_collection.get_sample_coverage_metadata(sample_name),
                sample_ploidy_metadata,
                self.ploidy_workspace.targets_metadata)

            # add to the collection
            self.ploidy_workspace.sample_metadata_collection.add_sample_ploidy_metadata(sample_ploidy_metadata)
            self.ploidy_workspace.sample_metadata_collection.add_sample_read_depth_metadata(sample_read_depth_metadata)
