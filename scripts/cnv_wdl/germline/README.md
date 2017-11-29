## Running the Germline CNV WDL

### Which WDL should you use?

- Calling a cohort of samples and building a model for denoising further case samples: ``cnv_germline_cohort_workflow.wdl``
- Calling a case sample using a previously built model for denoising: ``cnv_germline_case_workflow.wdl``

#### Setting up parameter json file for a run

To get started, copy the relevant ``*.template.json`` for the workflow you wish to run and adjust parameters accordingly.  

*Please note that there are optional workflow-level and task-level parameters that do not appear in the template files.  These are set to reasonable values by default, but can also be adjusted if desired.*

#### Required parameters in the germline cohort workflow

The reference used must be the same between PoN and case samples.

- ``CNVGermlineCohortWorkflow.gatk_docker`` -- GATK Docker image (e.g., ``broadinstitute/gatk:x.beta.x``).
- ``CNVGermlineCohortWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the chromosomes of interest.
- ``CNVGermlineCohortWorkflow.normal_bams_list`` -- TSV file consisting of corresponding bam and corresponding index files as described in cnv_germline_cohort_workflow.wdl.
- ``CNVGermlineCohortWorkflow.pon_entity_id`` -- Name of the final PoN file.
- ``CNVGermlineCohortWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVGermlineCohortWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVGermlineCohortWorkflow.ref_fasta`` -- Path to reference fasta file.

In additional, there are optional workflow-level and task-level parameters that may be set by advanced users; for example:

- ``CNVGermlineCohortWorkflow.do_explicit_gc_correction`` -- (optional) If true, perform explicit GC-bias correction when creating PoN and in subsequent denoising of case samples.  If false, rely on PCA-based denoising to correct for GC bias.
- ``CNVGermlineCohortWorkflow.PreprocessIntervals.bin_length`` -- Size of bins (in bp) for coverage collection.  *This must be the same value used for all case samples.*
- ``CNVGermlineCohortWorkflow.PreprocessIntervals.padding`` -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection.  *This must be the same value used for all case samples.*

Further explanation of other task-level parameters may be found by invoking the ``--help`` documentation available in the gatk.jar for each tool.  

#### Required parameters in the germline case workflow

The reference (and bins, if specified) used must be the same between PoN and case samples.

- ``CNVGermlineCaseWorkflow.gatk_docker`` -- GATK Docker image (e.g., "broadinstitute/gatk:x.beta.x").
- ``CNVGermlineCaseWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the chromosomes of interest.
- ``CNVGermlineCaseWorkflow.normal_bam`` -- File path or storage location (depending on backend) of the normal BAM file.
- ``CNVGermlineCaseWorkflow.normal_bam_idx`` -- File path or storage location (depending on backend) of the normal BAM file index.
- ``CNVGermlineCaseWorkflow.read_count_pon`` -- Path to read-count PoN created by the cohort workflow. 
- ``CNVGermlineCaseWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVGermlineCaseWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVGermlineCaseWorkflow.ref_fasta`` -- Path to reference fasta file.
- ``CNVGermlineCaseWorkflow.tumor_bam`` -- File path or storage location (depending on backend) of the tumor BAM file.
- ``CNVGermlineCaseWorkflow.tumor_bam_idx`` -- File path or storage location (depending on backend) of the tumor BAM file index.

In additional, there are several task-level parameters that may be set by advanced users as above.

Further explanation of these task-level parameters may be found by invoking the ``--help`` documentation available in the gatk.jar for each tool.