import os
import logging
from ..inference.fancy_optimizers import FancyAdamax
from . import io_commons
from . import io_consts

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

        mu_m = self.fancy_adamax.m_adam[0].get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "mu_" + io_consts.default_adamax_m_filename), mu_m)

        mu_u = self.fancy_adamax.u_adam[0].get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "mu_" + io_consts.default_adamax_u_filename), mu_u)

        rho_m = self.fancy_adamax.m_adam[1].get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "rho_" + io_consts.default_adamax_m_filename), rho_m)

        rho_u = self.fancy_adamax.u_adam[1].get_value(borrow=True)
        io_commons.write_ndarray_to_tsv(
            os.path.join(self.output_path, "rho_" + io_consts.default_adamax_u_filename), rho_u)


# todo assertions
class AdamaxMomentEstimateImporter:
    def __init__(self,
                 fancy_adamax: FancyAdamax,
                 input_path: str):
        self.fancy_adamax = fancy_adamax
        self.input_path = input_path

    def __call__(self):
        _logger.info("Importing adamax moment estimates...")

        self.fancy_adamax.m_adam[0].set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, "mu_" + io_consts.default_adamax_m_filename)),
            borrow=True)

        self.fancy_adamax.m_adam[1].set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, "rho_" + io_consts.default_adamax_m_filename)),
            borrow=True)

        self.fancy_adamax.u_adam[0].set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, "mu_" + io_consts.default_adamax_u_filename)),
            borrow=True)

        self.fancy_adamax.u_adam[1].set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, "rho_" + io_consts.default_adamax_u_filename)),
            borrow=True)
