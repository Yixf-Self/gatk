import logging

import numpy as np
import pymc3 as pm
import os
import json
from typing import List, Optional

from .._version import __version__ as gcnvkernel_version
from .. import config
from ..models.model_denoising_calling import DenoisingCallingWorkspace, DenoisingModel
from ..models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig
from . import io_commons
from . import io_consts

_logger = logging.getLogger(__name__)


class DenoisingModelExporter:
    """ Writes denoising model parameters to disk """
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 output_path: str):
        io_commons.assert_output_path_writable(output_path)
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path

    @staticmethod
    def _export_class_log_posterior(output_path, log_q_tau_tk):
        io_commons.write_ndarray_to_tsv(
            os.path.join(output_path, io_consts.default_class_log_posterior_tsv_filename), log_q_tau_tk)

    @staticmethod
    def _export_dict_to_json_file(output_file, raw_dict, blacklisted_keys):
        filtered_dict = {k: v for k, v in raw_dict.items() if k not in blacklisted_keys}
        with open(output_file, 'w') as fp:
            json.dump(filtered_dict, fp, indent=1)

    def __call__(self):
        # export gcnvkernel version
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_gcnvkernel_version_json_filename),
            {'version': gcnvkernel_version}, {})

        # export denoising config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_denoising_config_json_filename),
            self.denoising_config.__dict__, {})

        # export calling config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_calling_config_json_filename),
            self.calling_config.__dict__, {})

        # export global variables in the workspace
        self._export_class_log_posterior(
            self.output_path, self.denoising_calling_workspace.log_q_tau_tk.get_value(borrow=True))

        # export global variables in the posterior
        io_commons.export_meanfield_global_params(
            self.output_path, self.denoising_model_approx, self.denoising_model)


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
        # check if the model is created with the same gcnvkernel version
        io_commons.check_gcnvkernel_version(
            os.path.join(self.input_path, io_consts.default_gcnvkernel_version_json_filename))

        # import global workspace variables
        self.denoising_calling_workspace.log_q_tau_tk.set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, io_consts.default_class_log_posterior_tsv_filename)),
            borrow=config.borrow_numpy)

        # import global posterior parameters
        io_commons.import_meanfield_global_params(
            self.input_path, self.denoising_model_approx, self.denoising_model)


class SampleDenoisingAndCallingPosteriorsExporter:
    """ Performs writing and reading of calls """
    def __init__(self,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 sample_names: List[str],
                 output_path: str):
        io_commons.assert_output_path_writable(output_path)
        assert len(sample_names) == denoising_calling_workspace.num_samples
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path
        self.sample_names = sample_names
        (self._approx_var_set, self._approx_mu_map,
         self._approx_std_map) = io_commons.extract_meanfield_posterior_parameters(self.denoising_model_approx)

    @staticmethod
    def _export_sample_copy_number_log_posterior(sample_posterior_path: str,
                                                 log_q_c_tc: np.ndarray,
                                                 delimiter='\t',
                                                 comment='#',
                                                 extra_comment_lines: Optional[List[str]] = None):
        assert isinstance(log_q_c_tc, np.ndarray)
        assert log_q_c_tc.ndim == 2
        num_copy_number_states = log_q_c_tc.shape[1]
        copy_number_header_columns = [io_consts.copy_number_column_prefix + str(cn)
                                      for cn in range(num_copy_number_states)]
        with open(os.path.join(sample_posterior_path,
                               io_consts.default_copy_number_log_posterior_tsv_filename), 'w') as f:
            if extra_comment_lines is not None:
                for comment_line in extra_comment_lines:
                    f.write(comment + comment_line + '\n')
            f.write(delimiter.join(copy_number_header_columns) + '\n')
            for ti in range(log_q_c_tc.shape[0]):
                f.write(delimiter.join([repr(x) for x in log_q_c_tc[ti, :]]) + '\n')

    @staticmethod
    def _export_sample_name(sample_posterior_path: str,
                            sample_name: str):
        with open(os.path.join(sample_posterior_path, io_consts.default_sample_name_txt_filename), 'w') as f:
            f.write(sample_name + '\n')

    def __call__(self):
        for si, sample_name in enumerate(self.sample_names):
            sample_name_comment_line = [io_consts.sample_name_header_prefix + sample_name]
            sample_posterior_path = os.path.join(self.output_path, io_consts.sample_folder_prefix + repr(si))
            _logger.info("Saving posteriors for sample \"{0}\" in \"{1}\"...".format(
                sample_name, sample_posterior_path))
            io_commons.assert_output_path_writable(sample_posterior_path, try_creating_output_path=True)

            # export sample-specific posteriors in the approximation
            io_commons.export_meanfield_sample_specific_params(
                si, sample_posterior_path, self._approx_var_set, self._approx_mu_map, self._approx_std_map,
                self.denoising_model, sample_name_comment_line)

            # export sample name
            self._export_sample_name(sample_posterior_path, sample_name)

            # export copy number posterior
            self._export_sample_copy_number_log_posterior(
                sample_posterior_path,
                self.denoising_calling_workspace.log_q_c_stc.get_value(borrow=True)[si, ...],
                extra_comment_lines=sample_name_comment_line)
