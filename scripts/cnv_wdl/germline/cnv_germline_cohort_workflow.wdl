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
    String cohort_entity_id
    File contig_ploidy_priors
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    Int num_intervals_per_scatter
    String gatk_docker
    File? gatk4_jar_override
    Int? mem_for_determine_germline_contig_ploidy
    Int? mem_for_germline_cnv_caller

    # If true, AnnotateIntervals will be run to create GC annotations and explicit GC correction
    # will be performed by the model generated by
    Boolean? do_explicit_gc_correction

    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
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

    call DetermineGermlineContigPloidyCohortMode {
        input:
            cohort_entity_id = cohort_entity_id,
            read_count_files = CollectCounts.counts,
            contig_ploidy_priors = contig_ploidy_priors,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem = mem_for_determine_germline_contig_ploidy
    }

    call CNVTasks.ScatterIntervals {
        input:
            interval_list = PreprocessIntervals.preprocessed_intervals,
            num_intervals_per_scatter = num_intervals_per_scatter,
            gatk_docker = gatk_docker
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCohortMode {
            input:
                scatter_index = scatter_index,
                cohort_entity_id = cohort_entity_id,
                read_count_files = CollectCounts.counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
                ref_fasta_dict = ref_fasta_dict,
                annotated_intervals = AnnotateIntervals.annotated_intervals,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem = mem_for_germline_cnv_caller
        }
    }
}

task DetermineGermlineContigPloidyCohortMode {
    String cohort_entity_id
    Array[File] read_count_files
    File contig_ploidy_priors
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    Int machine_mem = if defined(mem) then select_first([mem]) else 8
    Float command_mem = machine_mem - 0.5

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${machine_mem}g -jar $GATK_JAR DetermineGermlineContigPloidy \
            --input ${sep=" --input " read_count_files} \
            --contigPloidyPriors ${contig_ploidy_priors} \
            --output ${output_dir_} \
            --outputPrefix ${cohort_entity_id} \
            --verbosity DEBUG

        tar czf ${cohort_entity_id}-contig-ploidy-model.tar.gz -C ${output_dir_}/${cohort_entity_id}-model .
        tar czf ${cohort_entity_id}-contig-ploidy-calls.tar.gz -C ${output_dir_}/${cohort_entity_id}-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: command_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File contig_ploidy_model_tar = "${cohort_entity_id}-contig-ploidy-model.tar.gz"
        File contig_ploidy_calls_tar = "${cohort_entity_id}-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCohortMode {
    Int scatter_index
    String cohort_entity_id
    Array[File] read_count_files
    File contig_ploidy_calls_tar
    File intervals
    File ref_fasta_dict
    File? annotated_intervals
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    Int machine_mem = if defined(mem) then select_first([mem]) else 8
    Float command_mem = machine_mem - 0.5

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        mkdir contig-ploidy-calls
        tar xzf ${contig_ploidy_calls_tar} -C contig-ploidy-calls

        java -Xmx${machine_mem}g -jar $GATK_JAR GermlineCNVCaller \
            -L ${intervals} \
            --sequenceDictionary ${ref_fasta_dict} \
            --input ${sep=" --input " read_count_files} \
            --contigPloidyCalls contig-ploidy-calls \
            ${"--annotatedIntervals " + annotated_intervals} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output ${output_dir_} \
            --outputPrefix ${cohort_entity_id} \
            --verbosity DEBUG

        tar czf ${cohort_entity_id}-gcnv-model-${scatter_index}.tar.gz -C ${output_dir_}/${cohort_entity_id}-model .
        tar czf ${cohort_entity_id}-gcnv-calls-${scatter_index}.tar.gz -C ${output_dir_}/${cohort_entity_id}-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: command_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File gcnv_model_tar = "${cohort_entity_id}-gcnv-model-${scatter_index}.tar.gz"
        File gcnv_calls_tar = "${cohort_entity_id}-gcnv-calls-${scatter_index}.tar.gz"
    }
}