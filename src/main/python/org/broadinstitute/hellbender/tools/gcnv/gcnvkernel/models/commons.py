import numpy as np
import logging
import theano.tensor as tt
import pymc3 as pm
import theano as th
import pymc3.distributions.dist_math as pm_dist_math
from typing import Tuple

_logger = logging.getLogger(__name__)

_log_2_pi = 1.837877066409345  # np.log(2 * np.pi)
_10_inv_log_10 = 4.342944819032518  # 10 / np.log(10)


def get_normalized_prob_vector(prob_vector: np.ndarray, prob_sum_tol: float) -> np.ndarray:
    """
    todo
    :param prob_vector:
    :param prob_sum_tol:
    :return:
    """
    assert all(prob_vector >= 0), "probabilities must be non-negative"
    prob_sum = np.sum(prob_vector)
    if np.abs(prob_sum - 1.0) < prob_sum_tol:
        return prob_vector
    else:
        _logger.warning("The given probability vector ({0}) was not normalized to unity within the provided "
                        "tolerance ({1}); sum = {2}; normalizing and proceeding.".format(
            prob_vector, prob_sum_tol, prob_sum))
        return prob_vector / prob_sum


def poisson_logp(mu, value):
    """
    Poisson log probability

    :param mu: poisson mean
    :param value: observed
    :return: theano tensor
    """
    return pm_dist_math.bound(
        pm_dist_math.logpow(mu, value) - pm_dist_math.factln(value) - mu,
        mu > 0, value >= 0)


def negative_binomial_logp(mu, alpha, value):
    """
    Negative binomial log probability

    :param mu: mean
    :param alpha: inverse over-dispersion
    :param value: observed
    :return: theano tensor
    """
    return pm_dist_math.bound(pm_dist_math.binomln(value + alpha - 1, value)
                              + pm_dist_math.logpow(mu / (mu + alpha), value)
                              + pm_dist_math.logpow(alpha / (mu + alpha), alpha),
                              mu > 0, value >= 0, alpha > 0)


def negative_binomial_normal_approx_logp(mu, alpha, value):
    """
    Normal approximation to Negative binomial log probability

    :param mu: mean
    :param alpha: inverse over-dispersion
    :param value: observed
    :return: theano tensor
    """
    tau = alpha / (mu * (alpha + mu))  # precision
    return pm_dist_math.bound(0.5 * (tt.log(tau) - _log_2_pi - tau * tt.square(value - mu)),
                              mu > 0, value >= 0, alpha > 0)


# todo a lazy switch between normal approximation and exact
def negative_binomial_smart_approx_logp(mu, alpha, value, where):
    raise NotImplementedError


def centered_heavy_tail_logp(mu, value):
    """
    This distribution is obtained by taking X ~ Exp and performing a Bose transformation
    Y = (exp(X) - 1)^{-1}. The result is:

        p(y) = (1 + 2 \mu) y^{2\mu} (1 + y)^{-2(1 + \mu)}

    It is a heavy-tail distribution with non-existent first moment.

    :param mu: mode of the distribution
    :param value: observed
    :return: theano tensor
    """
    return pm_dist_math.bound(tt.log(1.0 + 2.0 * mu) + 2.0 * mu * tt.log(value)
                              - 2.0 * (1.0 + mu) * tt.log(1.0 + value),
                              mu >= 0, value > 0)


def safe_logaddexp(a, b):
    diff = b - a
    safe_diff = tt.switch(tt.isnan(diff), 0, diff)
    return tt.switch(safe_diff >= 0,
                     b + tt.log1p(tt.exp(-safe_diff)),
                     a + tt.log1p(tt.exp(safe_diff)))


def get_jensen_shannon_divergence(log_p_1, log_p_2):
    """
    todo
    :param log_p_1:
    :param log_p_2:
    :return:
    """
    p_1 = tt.exp(log_p_1)
    p_2 = tt.exp(log_p_2)
    diff_12 = p_1 - p_2
    log_diff_12 = log_p_1 - log_p_2
    safe_log_diff_12 = tt.switch(tt.isnan(log_diff_12), 0, log_diff_12)
    return 0.5 * tt.sum(diff_12 * safe_log_diff_12, axis=-1)


def get_hellinger_distance(log_p_1, log_p_2):
    p_1 = tt.exp(log_p_1)
    p_2 = tt.exp(log_p_2)
    return tt.sqrt(tt.sum(tt.square(tt.sqrt(p_1) - tt.sqrt(p_2)), axis=-1)) / tt.sqrt(2)


def perform_genotyping(log_p: np.ndarray) -> Tuple[int, float]:
    """
    Takes an array of log probs and perform genotyping

    Note:
        log_p must be properly normalized, i.e. np.sum(np.exp(log_p)) == 1
        (this is _not_ checked for performance)

    :param log_p: an array of log probs
    :return: a tuple (most likely genotype, phred-scaled genotyping quality; see below)
    """
    assert log_p.ndim == 1, "the log_p array must be a vector"
    assert log_p.size >= 2, "at least two states are required"
    sorted_log_p = sorted(enumerate(log_p), key=lambda x: -x[1])
    max_likely_genotype_idx = sorted_log_p[0][0]
    phred_genotype_quality = _10_inv_log_10 * (sorted_log_p[0][1] - sorted_log_p[1][1])
    return max_likely_genotype_idx, phred_genotype_quality


def node_posterior_mean_symbolic(approx: pm.MeanField, node, size=100,
                                 more_replacements=None, shape=None, dtype=None):
    """ Samples a given node over shared posterior and calculates the posterior mean
    :param approx: a pymc3 variational posterior
    :param node: a theano tensor or expression
    :param size: number of samples
    :param more_replacements: nodes to replace before sampling
    :param shape:
    :param dtype:
    :return: symbolic posterior mean
    """

    assert size > 0

    if shape is not None and dtype is not None:
        cum_sum = tt.zeros(shape, dtype)
    elif hasattr(node.tag, 'test_value') and node.tag.test_value is not None:
        cum_sum = tt.zeros(node.tag.test_value.shape, node.tag.test_value.dtype)
    else:
        raise Exception("Can not determine the shape of the node to be sampled")

    if more_replacements is not None:
        node = th.clone(node, more_replacements, strict=False)
    posterior_samples = approx.random(size)
    node = approx.to_flat_input(node)

    def add_sample_to_cum_sum(posterior_sample, _cum_sum):
        new_sample = th.clone(node, {approx.input: posterior_sample}, strict=False)
        return _cum_sum + tt.patternbroadcast(new_sample, _cum_sum.broadcastable)

    outputs, _ = th.scan(add_sample_to_cum_sum,
                         sequences=[posterior_samples],
                         outputs_info=[cum_sum],
                         n_steps=size)

    return outputs[-1] / size
