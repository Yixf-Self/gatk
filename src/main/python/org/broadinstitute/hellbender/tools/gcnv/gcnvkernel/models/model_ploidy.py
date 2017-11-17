import logging
import argparse
import inspect
from typing import List, Dict, Set, Tuple

import numpy as np
import pandas as pd
import pymc3 as pm
import theano as th
import theano.tensor as tt
from pymc3 import Normal, Deterministic, DensityDist, Bound, Exponential

from ..tasks.inference_task_base import HybridInferenceParameters, GeneralizedContinuousModel
from . import commons
from .. import config, types
from ..structs.interval import Interval
from ..structs.metadata import TargetsIntervalListMetadata, SampleMetadataCollection

_logger = logging.getLogger(__name__)


class PloidyModelConfig:
    _contig_name_column = 'CONTIG_NAME'
    _ploidy_prior_prefix = 'PLOIDY_PRIOR_'

    def __init__(self,
                 contig_ploidy_prior_map: Dict[str, np.ndarray] = None,
                 mean_bias_sd: float = 1e-2,
                 psi_j_scale: float = 1e-3,
                 psi_s_scale: float = 1e-4,
                 mapping_error_rate: float = 1e-2):
        """
        todo
        mention that good priors are important
        """
        assert contig_ploidy_prior_map is not None
        self.mean_bias_sd = mean_bias_sd
        self.psi_j_scale = psi_j_scale
        self.psi_s_scale = psi_s_scale
        self.mapping_error_rate = mapping_error_rate
        self.contig_ploidy_prior_map, self.num_ploidy_states = self._get_validated_contig_ploidy_prior_map(
            contig_ploidy_prior_map)
        self.contig_set = set(contig_ploidy_prior_map.keys())

    @staticmethod
    def _get_validated_contig_ploidy_prior_map(given_contig_ploidy_prior_map: Dict[str, np.ndarray],
                                               min_prob: float = 0) -> Tuple[Dict[str, np.ndarray], int]:
        given_contigs = set(given_contig_ploidy_prior_map.keys())
        num_ploidy_states: int = 0
        for contig in given_contigs:
            num_ploidy_states = max(num_ploidy_states, given_contig_ploidy_prior_map[contig].size)
        validated_contig_ploidy_prior_map: Dict[str, np.ndarray] = dict()
        for contig in given_contigs:
            padded_validated_prior = np.zeros((num_ploidy_states,), dtype=types.floatX) + min_prob
            given_prior = given_contig_ploidy_prior_map[contig].flatten()
            padded_validated_prior[:given_prior.size] = padded_validated_prior[:given_prior.size] + given_prior
            padded_validated_prior = commons.get_normalized_prob_vector(padded_validated_prior, config.prob_sum_tol)
            validated_contig_ploidy_prior_map[contig] = padded_validated_prior
        return validated_contig_ploidy_prior_map, num_ploidy_states

    @staticmethod
    def get_contig_ploidy_prior_map_from_tsv_file(input_path: str):
        contig_ploidy_prior_pd = pd.read_csv(input_path, delimiter='\t', comment='#')
        columns = [str(x) for x in contig_ploidy_prior_pd.columns.values]
        assert len(columns) > 1
        assert columns[0] == PloidyModelConfig._contig_name_column
        contig_list = [str(x) for x in contig_ploidy_prior_pd['CONTIG_NAME'].values]
        assert all([len(column) > len(PloidyModelConfig._ploidy_prior_prefix)
                    and column[:len(PloidyModelConfig._ploidy_prior_prefix)] == PloidyModelConfig._ploidy_prior_prefix
                    for column in columns[1:]])
        ploidy_values = [int(column[len(PloidyModelConfig._ploidy_prior_prefix):]) for column in columns[1:]]
        num_ploidy_states = np.max(ploidy_values) + 1
        contig_ploidy_prior_map: Dict[str, np.ndarray] = dict()
        for contig in contig_list:
            contig_ploidy_prior_map[contig] = np.zeros((num_ploidy_states,), dtype=types.floatX)
        for ploidy in range(num_ploidy_states):
            column_name = PloidyModelConfig._ploidy_prior_prefix + str(ploidy)
            if column_name in columns:
                values = [float(x) for x in contig_ploidy_prior_pd[column_name].values]
                for j, contig in enumerate(contig_list):
                    contig_ploidy_prior_map[contig][ploidy] = values[j]
        return contig_ploidy_prior_map

    @staticmethod
    def expose_args(args: argparse.ArgumentParser):
        group = args.add_argument_group(title="Germline contig ploidy determination model hyperparameters")
        initializer_params = inspect.signature(PloidyModelConfig.__init__).parameters

        group.add_argument("--mean_bias_sd",
                           type=float,
                           help="Contig-level mean bias standard deviation",
                           default=initializer_params['mean_bias_sd'].default)

        group.add_argument("--mapping_error_rate",
                           type=float,
                           help="Typical mapping error rate",
                           default=initializer_params['mapping_error_rate'].default)

        group.add_argument("--psi_j_scale",
                           type=float,
                           help="Global contig-level unexplained variance scale",
                           default=initializer_params['psi_j_scale'].default)

        group.add_argument("--psi_s_scale",
                           type=float,
                           help="Sample-specific contig-level unexplained variance scale",
                           default=initializer_params['psi_s_scale'].default)

    @staticmethod
    def from_args_dict(args_dict: Dict):
        relevant_keys = set(inspect.getfullargspec(PloidyModelConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return PloidyModelConfig(**relevant_kwargs)


class PloidyWorkspace:
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 targets_metadata: TargetsIntervalListMetadata,
                 sample_names: List[str],
                 sample_metadata_collection: SampleMetadataCollection):
        self.targets_metadata = targets_metadata
        self.sample_metadata_collection = sample_metadata_collection
        self.ploidy_config = ploidy_config
        self.num_contigs = targets_metadata.num_contigs
        self.sample_names = sample_names
        self.num_samples: int = len(sample_names)
        self.num_ploidy_states = ploidy_config.num_ploidy_states
        assert all([contig in ploidy_config.contig_set for contig in targets_metadata.contig_set]), \
            "Some contigs do not have ploidy priors; cannot continue."
        assert sample_metadata_collection.all_samples_have_coverage_metadata(sample_names), \
            "Some samples do not have coverage metadata; cannot continue."

        # number of targets per contig as a shared theano tensor
        self.t_j: types.TensorSharedVariable = th.shared(targets_metadata.t_j.astype(types.floatX),
                                                         name='t_j', borrow=config.borrow_numpy)

        # count per contig and total count as shared theano tensors
        n_sj = np.zeros((self.num_samples, self.num_contigs), dtype=types.floatX)
        n_s = np.zeros((self.num_samples,), dtype=types.floatX)
        for si, sample_name in enumerate(self.sample_names):
            sample_metadata = sample_metadata_collection.get_sample_coverage_metadata(sample_name)
            n_sj[si, :] = sample_metadata.n_j[:]
            n_s[si] = sample_metadata.n_total
        self.n_sj: types.TensorSharedVariable = th.shared(n_sj, name='n_sj', borrow=config.borrow_numpy)
        self.n_s: types.TensorSharedVariable = th.shared(n_s, name='n_s', borrow=config.borrow_numpy)

        # integer ploidy values
        int_ploidy_values_k = np.arange(0, ploidy_config.num_ploidy_states, dtype=types.small_uint)
        self.int_ploidy_values_k = th.shared(int_ploidy_values_k, name='int_ploidy_values_k',
                                             borrow=config.borrow_numpy)

        # ploidy priors
        p_ploidy_jk = np.zeros((self.num_contigs, self.ploidy_config.num_ploidy_states), dtype=types.floatX)
        for j, contig in enumerate(targets_metadata.contig_list):
            p_ploidy_jk[j, :] = ploidy_config.contig_ploidy_prior_map[contig][:]
        log_p_ploidy_jk = np.log(p_ploidy_jk)
        self.log_p_ploidy_jk: types.TensorSharedVariable = th.shared(log_p_ploidy_jk, name='log_p_ploidy_jk',
                                                                     borrow=config.borrow_numpy)

        # ploidy log posteriors (placeholder)
        #
        log_q_ploidy_sjk = np.tile(log_p_ploidy_jk, (self.num_samples, 1, 1))
        self.log_q_ploidy_sjk: types.TensorSharedVariable = th.shared(
            log_q_ploidy_sjk, name='log_q_ploidy_sjk', borrow=config.borrow_numpy)

        # ploidy log emission (placeholder)
        log_ploidy_emission_sjk = np.zeros(
            (self.num_samples, self.num_contigs, ploidy_config.num_ploidy_states), dtype=types.floatX)
        self.log_ploidy_emission_sjk: types.TensorSharedVariable = th.shared(
            log_ploidy_emission_sjk, name="log_ploidy_emission_sjk", borrow=config.borrow_numpy)

        # exclusion mask; mask(j, k) = 1 - delta(j, k)
        contig_exclusion_mask_jj = (np.ones((self.num_contigs, self.num_contigs), dtype=types.small_uint)
                                    - np.eye(self.num_contigs, dtype=types.small_uint))
        self.contig_exclusion_mask_jj = th.shared(contig_exclusion_mask_jj, name='contig_exclusion_mask_jj')

    @staticmethod
    def _get_contig_set_from_interval_list(targets_interval_list: List[Interval]) -> Set[str]:
        return {target.contig for target in targets_interval_list}


class PloidyModel(GeneralizedContinuousModel):
    PositiveNormal = Bound(Normal, lower=0)  # how cool is this?

    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__()

        # shorthands
        t_j = ploidy_workspace.t_j
        contig_exclusion_mask_jj = ploidy_workspace.contig_exclusion_mask_jj
        n_s = ploidy_workspace.n_s
        n_sj = ploidy_workspace.n_sj
        ploidy_k = ploidy_workspace.int_ploidy_values_k
        q_ploidy_sjk = tt.exp(ploidy_workspace.log_q_ploidy_sjk)
        eps = ploidy_config.mapping_error_rate

        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        # mean per-contig bias
        mean_bias_j = self.PositiveNormal('mean_bias_j',
                                          mu=1.0,
                                          sd=ploidy_config.mean_bias_sd,
                                          shape=(ploidy_workspace.num_contigs,))
        register_as_global(mean_bias_j)

        # contig coverage unexplained variance
        psi_j = Exponential(name='psi_j',
                            lam=1.0 / ploidy_config.psi_j_scale,
                            shape=(ploidy_workspace.num_contigs,))
        register_as_global(psi_j)

        # sample-specific contig unexplained variance
        psi_s = Exponential(name='psi_s',
                            lam=1.0 / ploidy_config.psi_j_scale,
                            shape=(ploidy_workspace.num_samples,))
        register_as_sample_specific(psi_s)

        # convert "unexplained variance" to negative binomial over-dispersion
        alpha_sj = tt.inv((tt.exp(psi_j.dimshuffle('x', 0) + psi_s.dimshuffle(0, 'x')) - 1.0))

        # mean ploidy per contig per sample
        mean_ploidy_sj = tt.sum(tt.exp(ploidy_workspace.log_q_ploidy_sjk)
                                * ploidy_workspace.int_ploidy_values_k.dimshuffle('x', 'x', 0), axis=2)

        # mean-field amplification coefficient per contig
        gamma_sj = mean_ploidy_sj * t_j.dimshuffle('x', 0) * mean_bias_j.dimshuffle('x', 0)

        # gamma_rest_sj \equiv sum_{j' \neq j} gamma_sj
        gamma_rest_sj = tt.dot(gamma_sj, contig_exclusion_mask_jj)

        # NB per-contig counts
        mu_num_sjk = (t_j.dimshuffle('x', 0, 'x') * mean_bias_j.dimshuffle('x', 0, 'x')
                      * ploidy_k.dimshuffle('x', 'x', 0))
        mu_den_sjk = gamma_rest_sj.dimshuffle(0, 1, 'x') + mu_num_sjk
        eps_j = eps * t_j / tt.sum(t_j)  # proportion of fragments erroneously mapped to contig j
        mu_sjk = ((1.0 - eps) * (mu_num_sjk / mu_den_sjk)
                  + eps_j.dimshuffle('x', 0, 'x')) * n_s.dimshuffle(0, 'x', 'x')

        def _get_logp_sjk(_n_sj):
            _logp_sjk = commons.negative_binomial_logp(mu_sjk,  # mean
                                                       alpha_sj.dimshuffle(0, 1, 'x'),  # over-dispersion
                                                       _n_sj.dimshuffle(0, 1, 'x'))  # contig counts
            return _logp_sjk

        DensityDist(name='n_sj_obs',
                    logp=lambda _n_sj: tt.sum(q_ploidy_sjk * _get_logp_sjk(_n_sj)),
                    observed=n_sj)

        # for log ploidy emission sampling
        Deterministic(name='logp_sjk', var=_get_logp_sjk(n_sj))

    def reset(self):
        pass


class PloidyEmissionBasicSampler:
    """ Draws posterior samples from the ploidy log emission probability for a given variational approximation to
    the ploidy determination model posterior """
    def __init__(self, ploidy_model: PloidyModel, samples_per_round: int):
        self.ploidy_model = ploidy_model
        self.samples_per_round = samples_per_round
        self._simultaneous_log_ploidy_emission_sampler = None

    def update_approximation(self, approx: pm.approximations.MeanField):
        self._simultaneous_log_ploidy_emission_sampler = \
            self._get_compiled_simultaneous_log_ploidy_emission_sampler(approx)

    def is_sampler_initialized(self):
        return self._simultaneous_log_ploidy_emission_sampler is not None

    def draw(self):
        return self._simultaneous_log_ploidy_emission_sampler()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_compiled_simultaneous_log_ploidy_emission_sampler(self, approx: pm.approximations.MeanField):
        """ For a given variational approximation, returns a compiled theano function that draws posterior samples
        from the log ploidy emission """
        log_ploidy_emission_sjk = commons.node_posterior_mean_symbolic(
            approx, self.ploidy_model['logp_sjk'], size=self.samples_per_round)
        return th.function(inputs=[], outputs=log_ploidy_emission_sjk)


class PloidyBasicCaller:
    """ Simple Bayesian update of contig ploidy log posteriors """
    def __init__(self,
                 inference_params: HybridInferenceParameters,
                 ploidy_workspace: PloidyWorkspace):
        self.ploidy_workspace = ploidy_workspace
        self.inference_params = inference_params
        self._update_log_q_ploidy_sjk_theano_func = self._get_update_log_q_ploidy_sjk_theano_func()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_log_q_ploidy_sjk_theano_func(self):
        new_log_q_ploidy_sjk = (self.ploidy_workspace.log_p_ploidy_jk.dimshuffle('x', 0, 1)
                                + self.ploidy_workspace.log_ploidy_emission_sjk)
        new_log_q_ploidy_sjk -= pm.logsumexp(new_log_q_ploidy_sjk, axis=2)
        old_log_q_ploidy_sjk = self.ploidy_workspace.log_q_ploidy_sjk
        admixed_new_log_q_ploidy_sjk = commons.safe_logaddexp(
            new_log_q_ploidy_sjk + np.log(self.inference_params.caller_admixing_rate),
            old_log_q_ploidy_sjk + np.log(1.0 - self.inference_params.caller_admixing_rate))
        update_norm_sj = commons.get_hellinger_distance(admixed_new_log_q_ploidy_sjk, old_log_q_ploidy_sjk)
        return th.function(inputs=[],
                           outputs=[update_norm_sj],
                           updates=[(self.ploidy_workspace.log_q_ploidy_sjk, admixed_new_log_q_ploidy_sjk)])

    def call(self) -> np.ndarray:
        return self._update_log_q_ploidy_sjk_theano_func()

