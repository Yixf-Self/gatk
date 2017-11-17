from pymc3 import __version__ as pymc3_version
from ._version import __version__

# inference tasks
from .tasks.inference_task_base import HybridInferenceParameters
from .tasks.task_ploidy import PloidyInferenceTask
from .tasks.task_cohort_denoising_calling import CohortDenoisingAndCallingTask
from .tasks.task_case_sample_calling import CaseSampleCallingTask

# model configs and workspaces
from .models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig, DenoisingCallingWorkspace
from .models.model_denoising_calling import DefaultInitialModelParametersSupplier as DefaultDenoisingModelInitializer
from .models.model_ploidy import PloidyModelConfig, PloidyWorkspace

# metadata
from .structs.metadata import TargetsIntervalListMetadata, SampleMetadataCollection, \
    SampleCoverageMetadata, SamplePloidyMetadata

# pre-processing and io
from . import preprocess
from .utils import io
from .structs.interval import Interval, inherit_interval_annotations

assert pymc3_version in ["3.1", "3.2"], "Only PyMC3 3.1 and 3.2 are supported"
