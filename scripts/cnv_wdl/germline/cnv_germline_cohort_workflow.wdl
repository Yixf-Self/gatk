# Workflow for creating a GATK GermlineCNVCaller denoising model and generating calls given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - Input file (normal_bams_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#    normal_bam_1    bam_idx_1
#    normal_bam_2    bam_idx_2
#    ...
#
# - The interval-list file is required for both WGS and WES workflows and should be a Picard or GATK-style interval list.
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning).
#   For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_germline_cohort_workflow.wdl myParameters.json
#   See cnv_germline_cohort_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCohortWorkflow {
    File intervals
    File normal_bams_list
    Array[Array[String]]+ normal_bams = read_tsv(normal_bams_list)
    String pon_entity_id
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker
    File? gatk4_jar_override
    Int? mem_for_create_read_count_pon

    # If true, AnnotateIntervals will be run to create GC annotations and explicit GC correction
    # will be performed by the model generated by
    Boolean? do_explicit_gc_correction

    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker
    }

    if (select_first([do_explicit_gc_correction, false])) {
        call CNVTasks.AnnotateIntervals {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker
        }
    }

    scatter (normal_bam in normal_bams) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam[0],
                bam_idx = normal_bam[1],
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker
        }
    }

    call Determine {
        input:
            pon_entity_id = pon_entity_id,
            read_count_files = CollectCounts.counts,
            annotated_intervals = AnnotateIntervals.annotated_intervals,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem = mem_for_create_read_count_pon
    }

    output {
        File read_count_pon = CreateReadCountPanelOfNormals.read_count_pon
    }
}

task CreateReadCountPanelOfNormals {
    String pon_entity_id
    Array[File] read_count_files
    Float? minimum_interval_median_percentile
    Float? maximum_zeros_in_sample_percentage
    Float? maximum_zeros_in_interval_percentage
    Float? extreme_sample_median_percentile
    Boolean? do_impute_zeros
    Float? extreme_outlier_truncation_percentile
    Int? number_of_eigensamples
    File? annotated_intervals   #do not perform explicit GC correction by default
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    Int machine_mem = if defined(mem) then select_first([mem]) else 8
    Float command_mem = machine_mem - 0.5

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${machine_mem}g -jar $GATK_JAR CreateReadCountPanelOfNormals \
            --input ${sep=" --input " read_count_files} \
            --minimumIntervalMedianPercentile ${default="10.0" minimum_interval_median_percentile} \
            --maximumZerosInSamplePercentage ${default="5.0" maximum_zeros_in_sample_percentage} \
            --maximumZerosInIntervalPercentage ${default="5.0" maximum_zeros_in_interval_percentage} \
            --extremeSampleMedianPercentile ${default="2.5" extreme_sample_median_percentile} \
            --doImputeZeros ${default="true" do_impute_zeros} \
            --extremeOutlierTruncationPercentile ${default="0.1" extreme_outlier_truncation_percentile} \
            --numberOfEigensamples ${default="20" number_of_eigensamples} \
            ${"--annotatedIntervals " + annotated_intervals} \
            --output ${pon_entity_id}.pon.hdf5
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: command_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File read_count_pon = "${pon_entity_id}.pon.hdf5"
    }
}