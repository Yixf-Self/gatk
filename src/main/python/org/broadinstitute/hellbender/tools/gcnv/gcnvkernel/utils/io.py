import logging
from typing import List, Optional, Tuple, Set
from ast import literal_eval as make_tuple

import numpy as np
import pandas as pd
import pymc3 as pm
import os
import json
from .._version import __version__ as gcnv_kernel_version
from .. import config

from ..structs.interval import Interval, IntervalAnnotation, interval_annotations_dict, interval_annotations_dtypes
from ..models.model_denoising_calling import DenoisingCallingWorkspace, DenoisingModel
from ..models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig
from .. import types
from collections import namedtuple

_logger = logging.getLogger(__name__)

# standard read counts and target interval list files data types
interval_dtypes_dict = {'CONTIG': np.str, 'START': types.med_uint, 'END': types.med_uint}
read_count_dtypes_dict = {**interval_dtypes_dict, 'COUNT': types.med_uint}


def load_read_counts_tsv_file(read_counts_tsv_file: str,
                              max_rows: Optional[int] = None,
                              output_targets: bool = False)\
        -> Tuple[np.ndarray, Optional[List[Interval]]]:
    counts_pd = pd.read_csv(read_counts_tsv_file, delimiter='\t', nrows=max_rows,
                            dtype={**read_count_dtypes_dict})
    if output_targets:
        targets_pd = counts_pd[list(interval_dtypes_dict.keys())]
        targets_interval_list = _convert_targets_pd_to_interval_list(targets_pd)
        return counts_pd['COUNT'].as_matrix(), targets_interval_list
    else:
        return counts_pd['COUNT'], None


def load_targets_tsv_file(targets_tsv_file: str) -> List[Interval]:
    targets_pd = pd.read_csv(targets_tsv_file, delimiter='\t',
                             dtype={**interval_dtypes_dict, **interval_annotations_dtypes})
    return _convert_targets_pd_to_interval_list(targets_pd)


def load_truncated_counts(sample_counts_table_tsv_file: str,
                          first_target_index: int,
                          last_target_index: int,
                          assert_targets_are_equal: bool = False) -> Tuple[List[str], np.ndarray, List[Interval]]:
    assert first_target_index >= 0
    assert last_target_index > first_target_index

    sample_counts_table_pd = pd.read_csv(sample_counts_table_tsv_file, delimiter='\t',
                                         dtype={'SAMPLE_NAME': str, 'PATH': str})
    columns = [str(x) for x in sample_counts_table_pd.columns.values]
    assert 'SAMPLE_NAME' in columns
    assert 'PATH' in columns

    sample_names = [x for x in sample_counts_table_pd['SAMPLE_NAME']]
    paths = [x for x in sample_counts_table_pd['PATH']]

    num_samples = len(sample_names)
    num_targets = last_target_index - first_target_index
    n_st = np.zeros((num_samples, num_targets), dtype=types.med_uint)

    targets_interval_list = None
    for si, sample_coverage_path in enumerate(paths):
        sample_coverage_t, sample_targets_interval_list = load_read_counts_tsv_file(
            sample_coverage_path, output_targets=targets_interval_list is None or assert_targets_are_equal)
        if assert_targets_are_equal and targets_interval_list is not None:
            assert targets_interval_list == sample_targets_interval_list[first_target_index:last_target_index]
        if targets_interval_list is None:
            targets_interval_list = sample_targets_interval_list[first_target_index:last_target_index]
        n_st[si, :] = sample_coverage_t[first_target_index:last_target_index]

    return sample_names, n_st, targets_interval_list


def _convert_targets_pd_to_interval_list(targets_pd: pd.DataFrame) -> List[Interval]:
    """
    Converts a pandas dataframe targets intervals to list(Interval). Annotations will be parsed
    and added to the intervals as well.
    """
    interval_list: List[Interval] = []
    columns = [str(x) for x in targets_pd.columns.values]
    assert all([required_column in columns
                for required_column in interval_dtypes_dict.keys()]), "Some columns missing"
    for contig, start, end in zip(targets_pd['CONTIG'], targets_pd['START'], targets_pd['END']):
        interval = Interval(contig, start, end)
        interval_list.append(interval)

    for annotation_key in set(columns).intersection(interval_annotations_dict.keys()):
        bad_annotations_found = False
        for ti, raw_value in enumerate(targets_pd[annotation_key]):
            try:
                annotation: IntervalAnnotation = interval_annotations_dict[annotation_key](raw_value)
                interval_list[ti].add_annotation(annotation_key, annotation)
            except ValueError:
                bad_annotations_found = True
        if bad_annotations_found:
            _logger.warning("Some of the annotations for {0} contained bad values and were ignored".format(
                annotation_key))

    return interval_list


def assert_output_path_writable(output_path: str, try_creating_output_path: bool = True):
    if os.path.exists(output_path):
        if not os.path.isdir(output_path):
            raise AssertionError("the provided output path ({0}) is not a directory")
    elif try_creating_output_path:
        try:
            os.makedirs(output_path)
        except IOError:
            raise AssertionError("the provided output path ({0}) does not exist and can not be created")
    tmp_prefix = "write_tester"
    count = 0
    filename = os.path.join(output_path, tmp_prefix)
    while os.path.exists(filename):
        filename = "{}.{}".format(os.path.join(output_path, tmp_prefix), count)
        count = count + 1
    try:
        filehandle = open(filename, 'w')
        filehandle.close()
        os.remove(filename)
    except IOError:
        raise AssertionError("the output path ({0}) is not writeable".format(output_path))


def write_ndarray_to_npy(output_file: str, array: np.ndarray) -> None:
    np.save(output_file, array, allow_pickle=False, fix_imports=False)


def read_ndarray_from_npy(input_file: str) -> np.ndarray:
    return np.load(input_file, allow_pickle=False, fix_imports=False)


def write_ndarray_to_tsv(output_file: str, array: np.ndarray, comment='#', delimiter='\t') -> None:
    array = np.asarray(array)
    assert array.ndim <= 2
    shape = array.shape
    dtype = array.dtype
    if array.ndim == 2:
        array_matrix = array
    else:
        array_matrix = array.reshape((array.size, 1))

    with open(output_file, 'w') as f:
        f.write(comment + ' shape=' + repr(shape) + '\n')
        f.write(comment + ' dtype=' + str(dtype) + '\n')
        for i_row in range(array_matrix.shape[0]):
            row = array_matrix[i_row, :]
            row_repr = delimiter.join([repr(x) for x in row])
            f.write(row_repr + '\n')


def read_ndarray_from_tsv(input_file: str, comment='#', delimiter='\t') -> np.ndarray:
    dtype = None
    shape = None
    rows: List[np.ndarray] = []

    def _get_value(key: str, _line: str):
        key_loc = _line.find(key)
        if key_loc >= 0:
            val_loc = _line.find('=')
            return _line[val_loc + 1:].strip()
        else:
            return None

    with open(input_file, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            elif stripped_line[0] == comment:
                if dtype is None:
                    dtype = _get_value('dtype', stripped_line)
                if shape is None:
                    shape = _get_value('shape', stripped_line)
            else:
                assert dtype is not None and shape is not None,\
                    "shape and dtype information could not be found in the comments; " \
                    "cannot continue loading {0}".format(input())
                row = np.asarray(stripped_line.split(delimiter), dtype=dtype)
                rows.append(row)

    return np.vstack(rows).reshape(make_tuple(shape))


def write_interval_list_to_tsv_file(output_file: str, interval_list: List[Interval]):
    assert len(interval_list) > 0, "can not write an empty interval list to disk"
    annotation_found_keys: Set[str] = {}
    for interval in interval_list:
        for key in interval.annotations.keys():
            annotation_found_keys.add(key)
    mutual_annotation_key_list: List[str] = []
    for key in annotation_found_keys:
        if all(key in interval.annotations.keys() for interval in interval_list):
            mutual_annotation_key_list.append(key)
        else:
            _logger.warning("Only some targets have annotation ({0}) and others do not; "
                            "cannot write this annotation to disk; proceeding...")
    with open(output_file, 'w') as out:
        header = '\t'.join(['CONTIG', 'START', 'END'] + mutual_annotation_key_list)
        out.write(header + '\n')
        for interval in interval_list:
            row = '\t'.join([interval.contig, repr(interval.start), repr(interval.end)] +
                            [repr(interval.annotations[key]) for key in mutual_annotation_key_list])
            out.write(row + '\n')


def _get_vmap_from_approx(approx: pm.MeanField):
    if pm.__version__ == "3.1":
        return approx.gbij.ordering.vmap
    elif pm.__version__ == "3.2":
        return approx.bij.ordering.vmap
    else:
        raise Exception("Unsupported PyMC3 version")


def _extract_meanfield_posterior_parameters(approx: pm.MeanField):
    mu_flat_view = approx.mean.get_value()
    std_flat_view = approx.std.eval()
    mu_map = dict()
    std_map = dict()
    var_set = set()
    for vmap in _get_vmap_from_approx(approx):
        var_set.add(vmap.var)
        mu_map[vmap.var] = mu_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
        std_map[vmap.var] = std_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
    return var_set, mu_map, std_map


class DenoisingModelExporter:
    """ Writes denoising model parameters to disk """
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 output_path: str):
        assert_output_path_writable(output_path)
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path
        self._approx_var_set, self._approx_mu_map, self._approx_std_map = _extract_meanfield_posterior_parameters(
            self.denoising_model_approx)

    @staticmethod
    def _export_class_log_posterior(output_path, log_q_tau_tk):
        write_ndarray_to_tsv(os.path.join(output_path, "log_q_tau_tk.tsv"), log_q_tau_tk)

    @staticmethod
    def _export_dict_to_json_file(output_file, raw_dict, blacklisted_keys):
        filtered_dict = {k: v for k, v in raw_dict.items() if k not in blacklisted_keys}
        with open(output_file, 'w') as fp:
            json.dump(filtered_dict, fp, indent=1)

    def __call__(self):
        # export gcnvkernel version
        self._export_dict_to_json_file(
            os.path.join(self.output_path, "gcnv_kernel_version.json"), {'version': gcnv_kernel_version}, {})

        # export denoising config
        self._export_dict_to_json_file(
            os.path.join(self.output_path, "denoising_config.json"), self.denoising_config.__dict__, {})

        # export calling config
        self._export_dict_to_json_file(
            os.path.join(self.output_path, "calling_config.json"), self.calling_config.__dict__, {})

        # export interval list
        write_interval_list_to_tsv_file(os.path.join(self.output_path, "interval_list.tsv"),
                                        self.denoising_calling_workspace.targets_interval_list)

        # export global variables in the workspace
        self._export_class_log_posterior(
            self.output_path, self.denoising_calling_workspace.log_q_tau_tk.get_value(borrow=True))

        # export global variables in the posterior
        for var_name in self.denoising_model.global_var_registry:
            assert var_name in self._approx_var_set,\
                "a variable named {0} does not exist in the approximation".format(var_name)
            _logger.info("exporting {0}...".format(var_name))
            var_mu = self._approx_mu_map[var_name]
            var_std = self._approx_std_map[var_name]
            var_mu_out_path = os.path.join(self.output_path, 'mu_' + var_name + '.tsv')
            write_ndarray_to_tsv(var_mu_out_path, var_mu)
            var_std_out_path = os.path.join(self.output_path, 'std_' + var_name + '.tsv')
            write_ndarray_to_tsv(var_std_out_path, var_std)


class DenoisingModelImporter:
    """ Reads denoising model parameters from disk """
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 input_path: str):
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.input_path = input_path

    def __call__(self):
        # import kernel version
        with open(os.path.join(self.input_path, "gcnv_kernel_version.json"), 'r') as fp:
            imported_gcnv_kernel_version = json.load(fp)['version']

        if imported_gcnv_kernel_version != gcnv_kernel_version:
            _logger.warning("The exported model is created with a different gcnvkernel version {exported: (0}, "
                            "current: {1}); backwards compatibility is not guaranteed; proceed at your own "
                            "risk!").format(imported_gcnv_kernel_version, gcnv_kernel_version)

        # import denoising config
        with open(os.path.join(self.input_path, "denoising_config.json"), 'r') as fp:
            imported_denoising_config_dict = json.load(fp)
            # todo do we want to override?

        # import calling config
        with open(os.path.join(self.input_path, "calling_config.json"), 'r') as fp:
            imported_calling_config_dict = json.load(fp)
            # todo do we want to override?

        # import interval list
        imported_interval_list = load_targets_tsv_file(os.path.join(self.input_path, "interval_list.tsv"))
        assert imported_interval_list == self.denoising_calling_workspace.targets_interval_list,\
            "The interval list in the exported model is different from the one in the workspace; cannot continue"
        # Note: unless the exported model has been hampered with, we can be sure that number of targets
        # will match hereafter for all model parameters

        # import global workspace variables
        self.denoising_calling_workspace.log_q_tau_tk.set_value(
            read_ndarray_from_tsv(os.path.join(self.input_path, "log_q_tau_tk.tsv")), borrow=config.borrow_numpy)

        # import global posterior parameters
        vmap_list = _get_vmap_from_approx(self.denoising_model_approx)
        model_var_set = {vmap.var for vmap in vmap_list}

        def _update_param_inplace(param, slc, dtype, new_value):
            param[slc] = new_value.astype(dtype).flatten()
            return param

        model_mu = self.denoising_model_approx.params[0]
        model_rho = self.denoising_model_approx.params[1]

        for var_name in self.denoising_model.global_var_registry:
            assert var_name in model_var_set,\
                "A variable named ({0}) could not be found in the model vmap; cannot continue"
            var_mu_input_file = os.path.join(self.input_path, 'mu_' + var_name + '.tsv')
            var_std_input_file = os.path.join(self.input_path, 'std_' + var_name + '.tsv')
            if not os.path.exists(var_mu_input_file) or not os.path.exists(var_std_input_file):
                _logger.warning("Model parameter values for ({0}) could not be found in the provided exported model "
                                "path; ignoring and proceeding...".format(var_name))
            var_mu = read_ndarray_from_tsv(var_mu_input_file)
            var_std = read_ndarray_from_tsv(var_std_input_file)

            # convert std to rho, see pymc3.dist_math.sd2rho
            var_rho = np.log(np.exp(var_std) - 1)
            del var_std

            for vmap in vmap_list:
                if vmap.var == var_name:
                    assert var_mu.shape == vmap.shp
                    assert var_rho.shape == vmap.shp
                    model_mu.set_value(_update_param_inplace(
                        model_mu.get_value(borrow=True), vmap.slc, vmap.dtyp, var_mu), borrow=True)
                    model_rho.set_value(_update_param_inplace(
                        model_rho.get_value(borrow=True), vmap.slc, vmap.dtyp, var_rho), borrow=True)


class SamplePosteriorsExporter:
    """ Performs writing and reading of calls """
    def __init__(self,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 sample_names: List[str],
                 output_path: str):
        assert_output_path_writable(output_path)
        assert len(sample_names) == denoising_calling_workspace.num_samples
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path
        self.sample_names = sample_names
        self._approx_var_set, self._approx_mu_map, self._approx_std_map = _extract_meanfield_posterior_parameters(
            self.denoising_model_approx)

    _ModelExportRecipe = namedtuple('ModelExportRecipe', 'var_name, output_filename, slicer')

    _model_sample_specific_var_export_recipes = [
        _ModelExportRecipe('read_depth_s_log__', 'log_read_depth', lambda si, array: array[si]),
        _ModelExportRecipe('gamma_s_log__', 'log_sample_specific_variance', lambda si, array: array[si]),
        _ModelExportRecipe('z_su', 'bias_factor_loadings', lambda si, array: array[si, ...]),
        _ModelExportRecipe('z_sg', 'gc_bias_factor_loadings', lambda si, array: array[si, ...]),
    ]

    @staticmethod
    def _export_sample_copy_number_log_posterior(sample_posterior_path, log_q_c_tc):
        write_ndarray_to_tsv(os.path.join(sample_posterior_path, 'log_q_c_tc.tsv'), log_q_c_tc)

    def __call__(self):
        for si, sample_name in enumerate(self.sample_names):
            sample_posterior_path = os.path.join(self.output_path, sample_name)
            assert_output_path_writable(sample_posterior_path, try_creating_output_path=True)

            for var_name in self._approx_var_set:
                for export_recipe in self._model_sample_specific_var_export_recipes:
                    if export_recipe.var_name == var_name:
                        mean_out_file_name = os.path.join(
                            sample_posterior_path, export_recipe.output_filename + "_mean.tsv")
                        mean_array = export_recipe.slicer(si, self._approx_mu_map[var_name])
                        write_ndarray_to_tsv(mean_out_file_name, mean_array)

                        std_out_file_name = os.path.join(
                            sample_posterior_path, export_recipe.output_filename + "_std.tsv")
                        std_array = export_recipe.slicer(si, self._approx_std_map[var_name])
                        write_ndarray_to_tsv(std_out_file_name, std_array)
                        break

            self._export_sample_copy_number_log_posterior(
                sample_posterior_path,
                self.denoising_calling_workspace.log_q_c_stc.get_value(borrow=True)[si, ...])
