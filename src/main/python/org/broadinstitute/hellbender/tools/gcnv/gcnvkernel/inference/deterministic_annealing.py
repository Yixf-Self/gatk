import pymc3 as pm

Operator = pm.operators.Operator
Inference = pm.Inference
MeanField = pm.MeanField


class KLThermal(Operator):
    """ Operator based on Kullback-Leibler Divergence with Temperature """

    def __init__(self, approx, temperature=None):
        super().__init__(approx)
        assert temperature is not None
        self.temperature = temperature

    def apply(self, f):
        z = self.input
        return self.temperature * self.logq_norm(z) - self.logp_norm(z)


class ADVIDeterministicAnnealing(Inference):
    def __init__(self, local_rv=None, model=None,
                 cost_part_grad_scale=1,
                 scale_cost_to_minibatch=False,
                 random_seed=None, start=None,
                 temperature=None):

        assert temperature is not None
        super().__init__(
            KLThermal, MeanField, None,
            local_rv=local_rv,
            model=model,
            cost_part_grad_scale=cost_part_grad_scale,
            scale_cost_to_minibatch=scale_cost_to_minibatch,
            random_seed=random_seed,
            start=start,
            op_kwargs={'temperature': temperature})


