import numpy as np
import pymc3 as pm
from typing import Callable


class ParamTrackerConfig:
    def __init__(self):
        self.param_names = []
        self.inv_trans_list = []
        self.inv_trans_param_names = []

    def add(self,
            param_name: str,
            inv_trans: Callable[[np.ndarray], np.ndarray],
            inv_trans_param_name: str):
        self.param_names.append(param_name)
        self.inv_trans_list.append(inv_trans)
        self.inv_trans_param_names.append(inv_trans_param_name)


class ParamTracker:
    def __init__(self, param_tracker_config: ParamTrackerConfig):
        self._param_trans_dict = {}
        self.tracked_param_values_dict = {}
        for param_name, inv_trans, inv_trans_param_name in zip(
                param_tracker_config.param_names,
                param_tracker_config.inv_trans_list,
                param_tracker_config.inv_trans_param_names):
            self._param_trans_dict[param_name] = (inv_trans, inv_trans_param_name)
            self.tracked_param_values_dict[inv_trans_param_name] = []

    def _extract_param_mean(self, approx: pm.approximations.MeanField):
        mu_flat_view = approx.mean.get_value(borrow=True)
        vmap_list = approx.bij.ordering.vmap
        out = dict()
        for vmap in vmap_list:
            param_name = vmap.var
            if param_name in self._param_trans_dict.keys():
                bare_param_mean = mu_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
                inv_trans = self._param_trans_dict[param_name][0]
                inv_trans_param_name = self._param_trans_dict[param_name][1]
                if inv_trans is None:
                    out[inv_trans_param_name] = bare_param_mean
                else:
                    out[inv_trans_param_name] = inv_trans(bare_param_mean)
        return out

    def record(self, approx, _loss, _i):
        out = self._extract_param_mean(approx)
        for key in self.tracked_param_values_dict.keys():
            self.tracked_param_values_dict[key].append(out[key])

    __call__ = record

    def clear(self):
        for key in self.tracked_param_values_dict.keys():
            self.tracked_param_values_dict[key] = []

    def __getitem__(self, key):
        return self.tracked_param_values_dict[key]
