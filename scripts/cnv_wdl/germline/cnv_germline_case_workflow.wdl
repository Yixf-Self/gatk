# Workflow for running GATK GermlineCNVCaller on a single case sample. Supports both WGS and WES.
#
# Notes:
#
# - The interval-list file is required for both WGS and WES workflows and should be a Picard or GATK-style interval list.
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning).
#   For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_germline_case_workflow.wdl myParameters.json
#   See cnv_germline_case_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCaseWorkflow {
    File intervals
    File bam
    File bam_idx
    File contig_ploidy_model_tar
    Array[File]+ gcnv_model_tars
    Int num_intervals_per_scatter
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker
    File? gatk4_jar_override
    Int? mem_for_determine_germline_contig_ploidy
    Int? mem_for_germline_cnv_caller

    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker
    }

    call CNVTasks.CollectCounts {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = bam,
            bam_idx = bam_idx,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_count_files = CollectCounts.counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
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
        call GermlineCNVCallerCaseMode {
            input:
                scatter_index = scatter_index,
                read_count_files = CollectCounts.counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                gcnv_model_tar = gcnv_model_tars[scatter_index],
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
                ref_fasta_dict = ref_fasta_dict,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem = mem_for_germline_cnv_caller
        }
    }
}

task DetermineGermlineContigPloidyCaseMode {
    Array[File] read_count_files
    File contig_ploidy_model_tar
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

        mkdir input-contig-ploidy-model
        tar xzf ${contig_ploidy_model_tar} -C input-contig-ploidy-model

        java -Xmx${machine_mem}g -jar $GATK_JAR DetermineGermlineContigPloidy \
            --input ${sep=" --input " read_count_files} \
            --model input-contig-ploidy-model \
            --output ${output_dir_} \
            --outputPrefix case \
            --verbosity DEBUG

        tar czf case-contig-ploidy-calls.tar.gz -C ${output_dir_}/case-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: command_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    Int scatter_index
    Array[File] read_count_files
    File contig_ploidy_calls_tar
    File gcnv_model_tar
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

        mkdir gcnv-model
        tar xzf ${gcnv_model_tar} -C gcnv-model

        java -Xmx${machine_mem}g -jar $GATK_JAR GermlineCNVCaller \
            --input ${sep=" --input " read_count_files} \
            --contigPloidyCalls contig-ploidy-calls \
            --model gcnv-model \
            --output ${output_dir_} \
            --outputPrefix case \
            --verbosity DEBUG

        tar czf case-gcnv-calls-${scatter_index}.tar.gz -C ${output_dir_}/case-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: command_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File gcnv_calls_tar = "case-gcnv-calls-${scatter_index}.tar.gz"
    }
}