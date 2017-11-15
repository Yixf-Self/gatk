from pymc3.distributions.continuous import PositiveContinuous,\
    assert_negative_support, bound, gammaln, get_variable_name
import theano.tensor as tt


class PositiveFlatTop(PositiveContinuous):
    R"""
    Positive flat-top log-likelihood.

    ========  ============================
    Support   :math:`x \in [0, \infty)`
    ========  ============================

    Parameters
    ----------
    u : float
        upper bound (u > 0)
    k : float
        decay slope (k > 0)
    """

    def __init__(self, u, k, *args, **kwargs):
        super(PositiveFlatTop, self).__init__(*args, **kwargs)
        self.u = u = tt.as_tensor_variable(u)
        self.k = k = tt.as_tensor_variable(k)

        assert_negative_support(u, 'u', 'PositiveFlatTop')
        assert_negative_support(k, 'k', 'PositiveFlatTop')

    def random(self, point=None, size=None, repeat=None):
        raise NotImplementedError

    def logp(self, value):
        u = self.u
        k = self.k
        logp = -tt.pow(value / u, k) - tt.log(u) - gammaln(1 + tt.inv(k))
        return bound(logp, value > 0, u > 0, k > 0)

    def _repr_latex_(self, name=None, dist=None):
        if dist is None:
            dist = self
        u = dist.u
        k = dist.k
        return r'${} \sim \text{{FlatTop}}(\mathit{{u}}={0}, \mathit{{k}}={1})$'.format(
            name, get_variable_name(u), get_variable_name(k))


# todo borrowed from pymc3 3.2
class HalfFlat(PositiveContinuous):
    """Improper flat prior over the positive reals."""

    def __init__(self, *args, **kwargs):
        self._default = 1
        super(HalfFlat, self).__init__(defaults=('_default',), *args, **kwargs)

    def random(self, point=None, size=None, repeat=None):
        raise ValueError('Cannot sample from HalfFlat distribution')

    def logp(self, value):
        return bound(tt.zeros_like(value), value > 0)

    def _repr_latex_(self, name=None, dist=None):
        return r'${} \sim \text{{HalfFlat}}()$'.format(name)

