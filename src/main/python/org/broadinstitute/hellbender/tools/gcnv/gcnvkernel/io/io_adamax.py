import os
import logging
from ..inference.fancy_optimizers import FancyAdamax
from . import io_commons
from . import io_consts
from .. import config

_logger = logging.getLogger(__name__)


class AdamaxMomentEstimateExporter:
    def __init__(self,
                 fancy_adamax: FancyAdamax,
                 output_path: str):
        self.fancy_adamax = fancy_adamax
        self.output_path = output_path

    def __call__(self):
        _logger.info("Exporting adamax moment estimates...")
        io_commons.assert_output_path_writable(self.output_path)

        mu_m = self.fancy_adamax.get_mu_m().get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "mu_" + io_consts.default_adamax_m_filename), mu_m)

        rho_m = self.fancy_adamax.get_rho_m().get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "rho_" + io_consts.default_adamax_m_filename), rho_m)

        mu_u = self.fancy_adamax.get_mu_u().get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "mu_" + io_consts.default_adamax_u_filename), mu_u)

        rho_u = self.fancy_adamax.get_rho_u().get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "rho_" + io_consts.default_adamax_u_filename), rho_u)

        res = self.fancy_adamax.get_res_tensor().get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, io_consts.default_adamax_res_filename), res)


class AdamaxMomentEstimateImporter:
    def __init__(self,
                 fancy_adamax: FancyAdamax,
                 input_path: str):
        self.fancy_adamax = fancy_adamax
        self.input_path = input_path

    @staticmethod
    def _assert_shape(imported_ndarray, expected_shared_tensor):
        assert imported_ndarray.shape == expected_shared_tensor.get_value(borrow=True).shape, \
            "The imported adamax moments have a different shape ({0}) than expected ({1}). This can " \
            "occur if the imported moments correspond to a different model".format(
                imported_ndarray.shape, expected_shared_tensor.get_value(borrow=True).shape)

    def __call__(self):
        _logger.info("Importing adamax moment estimates...")

        imported_mu_m = io_commons.read_ndarray_from_tsv(
            os.path.join(self.input_path, "mu_" + io_consts.default_adamax_m_filename))
        self._assert_shape(imported_mu_m, self.fancy_adamax.get_mu_m())
        self.fancy_adamax.get_mu_m().set_value(imported_mu_m, borrow=config.borrow_numpy)

        imported_mu_u = io_commons.read_ndarray_from_tsv(
            os.path.join(self.input_path, "mu_" + io_consts.default_adamax_u_filename))
        self._assert_shape(imported_mu_u, self.fancy_adamax.get_mu_u())
        self.fancy_adamax.get_mu_u().set_value(imported_mu_u, borrow=config.borrow_numpy)

        imported_rho_m = io_commons.read_ndarray_from_tsv(
            os.path.join(self.input_path, "rho_" + io_consts.default_adamax_m_filename))
        self._assert_shape(imported_rho_m, self.fancy_adamax.get_rho_m())
        self.fancy_adamax.get_rho_m().set_value(imported_rho_m, borrow=config.borrow_numpy)

        imported_rho_u = io_commons.read_ndarray_from_tsv(
            os.path.join(self.input_path, "rho_" + io_consts.default_adamax_u_filename))
        self._assert_shape(imported_rho_u, self.fancy_adamax.get_rho_u())
        self.fancy_adamax.get_rho_u().set_value(imported_rho_u, borrow=config.borrow_numpy)

        imported_res = io_commons.read_ndarray_from_tsv(
            os.path.join(self.input_path, io_consts.default_adamax_res_filename))
        self._assert_shape(imported_res, self.fancy_adamax.get_res_tensor())
        self.fancy_adamax.get_res_tensor().set_value(imported_res, borrow=config.borrow_numpy)
