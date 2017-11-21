from pymc3 import __version__ as pymc3_version
from ._version import __version__

# inference tasks
from .tasks.inference_task_base import HybridInferenceParameters
from .tasks.task_cohort_ploidy_determination import CohortPloidyInferenceTask
from .tasks.task_cohort_denoising_calling import CohortDenoisingAndCallingTask
from .tasks.task_case_denoising_calling import CaseDenoisingCallingTask
from .tasks.task_case_ploidy_determination import CasePloidyInferenceTask

# model configs and workspaces
from .models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig, DenoisingCallingWorkspace
from .models.model_denoising_calling import DefaultInitialModelParametersSupplier as DefaultDenoisingModelInitializer
from .models.model_ploidy import PloidyModelConfig, PloidyWorkspace

# metadata
from .structs.metadata import TargetsIntervalListMetadata, SampleMetadataCollection,\
    SampleCoverageMetadata, SamplePloidyMetadata

# pre-processing and io
from . import preprocess
from .io import io_consts, io_ploidy, io_denoising_calling, io_intervals_and_counts, io_metadata
from .utils import cli_commons

assert pymc3_version in ["3.1", "3.2"], "Only PyMC3 3.1 and 3.2 are supported"
