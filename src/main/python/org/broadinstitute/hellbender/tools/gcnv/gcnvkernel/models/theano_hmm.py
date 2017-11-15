import numpy as np
import theano as th
import theano.tensor as tt
import pymc3 as pm
from typing import Optional, Tuple
from .. import types
from . import commons


class TheanoForwardBackward:
    """ Implementation of the forward-backward algorithm using th.scan """
    def __init__(self,
                 log_posterior_output: Optional[types.TensorSharedVariable],
                 admixing_rate: float,
                 include_alpha_beta_output: bool = False,
                 resolve_nans: bool = False):
        """
        :param log_posterior_output:
        :param admixing_rate: a float in range [0, 1] denoting the amount of the new posterior to admix with the
        old posterior
        """
        self.admixing_rate = admixing_rate
        self.include_alpha_beta_output = include_alpha_beta_output
        self.resolve_nans = resolve_nans
        assert 0.0 < admixing_rate <= 1.0, "Admixing rate must be in range (0, 1]"
        self.log_posterior_output = log_posterior_output
        self._forward_backward_theano_func = self._get_compiled_forward_backward_theano_func()

    def perform_forward_backward(self,
                                 num_states: int,
                                 log_prior_c: np.ndarray,
                                 log_trans_tcc: np.ndarray,
                                 log_emission_tc: np.ndarray,
                                 prev_log_posterior_tc: np.ndarray):
        return self._forward_backward_theano_func(
            num_states, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc)

    # todo update docstring
    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_forward_backward_theano_func(self):
        """ Returns a compiled theano function that updates log posterior probabilities.

        The theano function takes 5 inputs:
            num_states (integer scalar),
            log_prior_c (float vector),
            og_trans_tcc (float tensor3),
            log_emission_tc (float matrix)
            prev_log_posterior_tc (float matrix)

        If a log_posterior_output is provided in the class initializer, the return will be:
            [update_norm_t, log_data_likelihood]
        otherwise, the return will be:
            [admixed_log_posterior_tc, update_norm_t, log_data_likelihood]

        :return: a theano function
        """
        num_states = tt.iscalar('num_states')
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')
        prev_log_posterior_tc = tt.matrix('prev_log_posterior_tc')

        new_log_posterior_tc, log_data_likelihood_t, alpha_tc, beta_tc = self._get_symbolic_log_posterior(
            num_states, log_prior_c, log_trans_tcc, log_emission_tc, self.resolve_nans)

        admixed_log_posterior_tc = commons.safe_logaddexp(
            new_log_posterior_tc + np.log(self.admixing_rate),
            prev_log_posterior_tc + np.log(1.0 - self.admixing_rate))

        log_data_likelihood = log_data_likelihood_t[-1]  # in theory, they are all the same
        update_norm_t = commons.get_jensen_shannon_divergence(admixed_log_posterior_tc, prev_log_posterior_tc)

        ext_output = [alpha_tc, beta_tc] if self.include_alpha_beta_output else []
        inputs = [num_states, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc]
        if self.log_posterior_output is not None:
            return th.function(inputs=inputs,
                               outputs=[update_norm_t, log_data_likelihood] + ext_output,
                               updates=[(self.log_posterior_output, admixed_log_posterior_tc)])
        else:
            return th.function(inputs=inputs,
                               outputs=[admixed_log_posterior_tc, update_norm_t, log_data_likelihood] + ext_output)

    @staticmethod
    def _get_symbolic_log_posterior(num_states: tt.iscalar,
                                    log_prior_c: types.TheanoVector,
                                    log_trans_tcc: types.TheanoTensor3,
                                    log_emission_tc: types.TheanoMatrix,
                                    resolve_nans: bool):
        """ Returns a symbolic tensor for log posterior and log data likelihood
        :return: tuple of (log_posterior_probs, log_data_likelihood)
        """

        def calculate_next_alpha(c_log_trans_mat: types.TheanoMatrix, c_log_emission_vec: types.TheanoVector,
                                 p_alpha_vec: types.TheanoVector):
            """ Calculates the next entry on the forward table, alpha(t), from alpha(t-1)
            :param c_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                                    from the previous state at position t-1 and to the next state at position t,
                                    respectively
            :param c_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
            :param p_alpha_vec: a 1D tensor representing alpha(t-1)
            :return: a 1d tensor representing alpha(t)
            """
            mu = tt.tile(p_alpha_vec, (num_states, 1)) + c_log_trans_mat.T
            n_alpha_vec = c_log_emission_vec + pm.math.logsumexp(mu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(n_alpha_vec), -np.inf, n_alpha_vec)
            else:
                return n_alpha_vec

        def calculate_prev_beta(n_log_trans_mat: types.TheanoMatrix, n_log_emission_vec: types.TheanoVector,
                                n_beta_vec: types.TheanoVector):
            """ Calculates the previous entry on the backward table, beta(t-1), from beta(t)
            :param n_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                                    from the previous state at position t-1 and to the next state at position t,
                                    respectively
            :param n_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
            :param n_beta_vec: a 1d tensor representing beta(t)
            :return: a 1d tensor representing beta(t-1)
            """
            nu = tt.tile(n_beta_vec + n_log_emission_vec, (num_states, 1)) + n_log_trans_mat
            p_beta_vec = pm.math.logsumexp(nu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(p_beta_vec), -np.inf, p_beta_vec)
            else:
                return p_beta_vec

        # first entry of the forward table
        alpha_first = log_prior_c + log_emission_tc[0, :]

        # the rest of the forward table
        alpha_outputs, alpha_updates = th.scan(
            fn=calculate_next_alpha,
            sequences=[log_trans_tcc, log_emission_tc[1:, :]],
            outputs_info=[alpha_first])

        # concatenate with the first alpha
        alpha_full_outputs = tt.concatenate((alpha_first.dimshuffle('x', 0), alpha_outputs))

        # last entry of the backward table (zero for all states)
        beta_last = tt.zeros_like(log_prior_c)

        # the rest of the backward table
        beta_outputs, beta_updates = th.scan(
            fn=calculate_prev_beta,
            sequences=[log_trans_tcc, log_emission_tc[1:, :]],
            go_backwards=True,
            outputs_info=[beta_last])

        # concatenate with the last beta and reverse
        beta_full_outputs = tt.concatenate((beta_last.dimshuffle('x', 0), beta_outputs))[::-1, :]
        log_unnormalized_posterior_probs = alpha_full_outputs + beta_full_outputs
        log_data_likelihood = pm.math.logsumexp(log_unnormalized_posterior_probs, axis=1)
        log_posterior_probs = log_unnormalized_posterior_probs - log_data_likelihood

        return log_posterior_probs, log_data_likelihood.dimshuffle(0), alpha_full_outputs, beta_full_outputs
