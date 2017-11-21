import logging
from typing import List, Optional, Tuple, Set, Dict
from ast import literal_eval as make_tuple

import numpy as np
import pymc3 as pm
import os
import json
import re

from .._version import __version__ as gcnvkernel_version
from . import io_consts
from ..tasks.inference_task_base import GeneralizedContinuousModel
from collections import namedtuple

_logger = logging.getLogger(__name__)


def extract_sample_name_from_header(input_file: str,
                                    max_scan_lines: int = 100) -> str:
    with open(input_file, 'r') as f:
        for _ in range(max_scan_lines):
            line = f.readline()
            match = re.search(io_consts.sample_name_header_regexp, line, re.M)
            if match is None:
                continue
            groups = match.groups()
            return groups[0]
    raise Exception("Sample name could not be found in \"{0}\"".format(input_file))


def assert_output_path_writable(output_path: str, try_creating_output_path: bool = True):
    if os.path.exists(output_path):
        if not os.path.isdir(output_path):
            raise AssertionError("the provided output path \"{0}\" is not a directory")
    elif try_creating_output_path:
        try:
            os.makedirs(output_path)
        except IOError:
            raise AssertionError("the provided output path \"{0}\" does not exist and can not be created")
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
        raise AssertionError("the output path \"{0}\" is not writeable".format(output_path))


def write_ndarray_to_tsv(output_file: str, array: np.ndarray, comment='#', delimiter='\t',
                         extra_comment_lines: Optional[List[str]] = None) -> None:
    array = np.asarray(array)
    assert array.ndim <= 2
    shape = array.shape
    dtype = array.dtype
    if array.ndim == 2:
        array_matrix = array
    else:
        array_matrix = array.reshape((array.size, 1))

    with open(output_file, 'w') as f:
        f.write(comment + 'shape=' + repr(shape) + '\n')
        f.write(comment + 'dtype=' + str(dtype) + '\n')
        if extra_comment_lines is not None:
            for comment_line in extra_comment_lines:
                f.write(comment + comment_line + '\n')
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


def get_var_map_list_from_approx(approx: pm.MeanField):
    if pm.__version__ == "3.1":
        return approx.gbij.ordering.vmap
    elif pm.__version__ == "3.2":
        return approx.bij.ordering.vmap
    else:
        raise Exception("Unsupported PyMC3 version")


def extract_meanfield_posterior_parameters(approx: pm.MeanField)\
        -> Tuple[Set[str], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    mu_flat_view = approx.mean.get_value()
    std_flat_view = approx.std.eval()
    mu_map = dict()
    std_map = dict()
    var_set = set()
    for vmap in get_var_map_list_from_approx(approx):
        var_set.add(vmap.var)
        mu_map[vmap.var] = mu_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
        std_map[vmap.var] = std_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
    return var_set, mu_map, std_map


def export_dict_to_json_file(output_file, raw_dict, blacklisted_keys):
    filtered_dict = {k: v for k, v in raw_dict.items() if k not in blacklisted_keys}
    with open(output_file, 'w') as fp:
        json.dump(filtered_dict, fp, indent=1)


def check_gcnvkernel_version(gcnvkernel_version_json_file: str):
    with open(gcnvkernel_version_json_file, 'r') as fp:
        imported_gcnvkernel_version = json.load(fp)['version']
        if imported_gcnvkernel_version != gcnvkernel_version:
            _logger.warning("The exported model is created with a different gcnvkernel version {exported: (0}, "
                            "current: {1}\"; backwards compatibility is not guaranteed; proceed at your own "
                            "risk!").format(imported_gcnvkernel_version, gcnvkernel_version)


ModelExportRecipe = namedtuple('ModelExportRecipe', 'var_name, output_filename, slicer')


def export_sample_specific_posteriors(sample_index: int,
                                      sample_posterior_path: str,
                                      approx_var_name_set: Set[str],
                                      approx_mu_map: Dict[str, np.ndarray],
                                      approx_std_map: Dict[str, np.ndarray],
                                      export_recipes: List[ModelExportRecipe],
                                      extra_comment_lines: Optional[List[str]] = None):
    for var_name in approx_var_name_set:
        for export_recipe in export_recipes:
            if export_recipe.var_name == var_name:
                mean_out_file_name = os.path.join(
                    sample_posterior_path, export_recipe.output_filename + "_mean.tsv")
                mean_array = export_recipe.slicer(sample_index, approx_mu_map[var_name])
                write_ndarray_to_tsv(mean_out_file_name, mean_array,
                                     extra_comment_lines=extra_comment_lines)

                std_out_file_name = os.path.join(
                    sample_posterior_path, export_recipe.output_filename + "_std.tsv")
                std_array = export_recipe.slicer(sample_index, approx_std_map[var_name])
                write_ndarray_to_tsv(std_out_file_name, std_array,
                                     extra_comment_lines=extra_comment_lines)
                break


def import_global_posteriors(input_path: str,
                             approx: pm.MeanField,
                             model: GeneralizedContinuousModel):
    # import global posterior parameters
    vmap_list = get_var_map_list_from_approx(approx)
    model_var_set = {vmap.var for vmap in vmap_list}

    def _update_param_inplace(param, slc, dtype, new_value):
        param[slc] = new_value.astype(dtype).flatten()
        return param

    model_mu = approx.params[0]
    model_rho = approx.params[1]

    for var_name in model.global_var_registry:
        assert var_name in model_var_set,\
            "A variable named \"{0}\" could not be found in the model variable map; cannot continue"
        var_mu_input_file = os.path.join(input_path, 'mu_' + var_name + '.tsv')
        var_std_input_file = os.path.join(input_path, 'std_' + var_name + '.tsv')
        if not os.path.exists(var_mu_input_file) or not os.path.exists(var_std_input_file):
            _logger.warning("Model parameter values for \"{0}\" could not be found in the provided model "
                            "path; ignoring and proceeding...".format(var_name))
            continue
        _logger.info("Importing model parameter values for \"{0}\"...".format(var_name))
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
