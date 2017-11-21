import os
import shutil

# set theano flags
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore"

import logging
import argparse
import gcnvkernel

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="gCNV cohort denoising and calling tool",
                                 formatter_class=gcnvkernel.cli_commons.GCNVHelpFormatter)

# logging args
gcnvkernel.cli_commons.add_logging_args_to_argparse(parser)

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--modeling_interval_list",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Full interval list, possibly including extra annotations (in .tsv format)")

group.add_argument("--read_count_tsv_files",
                   type=str,
                   required=True,
                   nargs='+',  # one or more
                   default=argparse.SUPPRESS,
                   help="List of read count files in the cohort (in .tsv format; must include sample name header)")

group.add_argument("--ploidy_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="The path to the results of ploidy determination tool")

group.add_argument("--output_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write model parameters")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write CNV calls")

# add denoising config args
gcnvkernel.DenoisingModelConfig.expose_args(parser)

# add calling config args
gcnvkernel.CopyNumberCallingConfig.expose_args(parser)

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(
    parser,
    hide={
        "--disable_sampler",
        "--disable_caller"
    })

if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()
    gcnvkernel.cli_commons.set_logging_config_from_args(args)

    # load modeling interval list
    modeling_interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(args.modeling_interval_list)

    # load sample names, truncated counts, and interval list from the sample read counts table
    logger.info("Loading {0} read counts file(s)...".format(len(args.read_count_tsv_files)))
    sample_names, n_st = gcnvkernel.io_intervals_and_counts.load_counts_in_the_modeling_zone(
        args.read_count_tsv_files, modeling_interval_list)

    # load read depth and ploidy metadata
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    gcnvkernel.io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection, args.ploidy_calls_path)

    # setup sample contig ploidy array
    contigs_set = {target.contig for target in modeling_interval_list}
    baseline_copy_number_s = None
    for contig in contigs_set:
        if baseline_copy_number_s is None:
            baseline_copy_number_s = sample_metadata_collection.get_sample_contig_ploidy_array(
                contig, sample_names)
        else:  # the target interval list has more than one contig
            other_baseline_copy_number_s = sample_metadata_collection.get_sample_contig_ploidy_array(
                contig, sample_names)
            assert all(baseline_copy_number_s == other_baseline_copy_number_s), \
                "Contig ploidy of one of more samples is variable across targets; cannot continue."

    # read depth array
    read_depth_s = sample_metadata_collection.get_sample_read_depth_array(sample_names)

    # setup the inference task
    args_dict = args.__dict__

    denoising_config = gcnvkernel.DenoisingModelConfig.from_args_dict(args_dict)
    calling_config = gcnvkernel.CopyNumberCallingConfig.from_args_dict(args_dict)
    inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)
    shared_workspace = gcnvkernel.DenoisingCallingWorkspace(
        denoising_config, calling_config, modeling_interval_list,
        n_st, baseline_copy_number_s, read_depth_s)
    initial_params_supplier = gcnvkernel.DefaultDenoisingModelInitializer(
        denoising_config, calling_config, shared_workspace)
    task = gcnvkernel.CohortDenoisingAndCallingTask(
        denoising_config, calling_config, inference_params,
        shared_workspace, initial_params_supplier)

    # go!
    task.engage()
    task.disengage()

    # save model
    gcnvkernel.io_denoising_calling.DenoisingModelExporter(
        denoising_config, calling_config,
        shared_workspace, task.continuous_model, task.continuous_model_approx,
        args.output_model_path)()

    # save a copy of targets in the model path
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_interval_list_filename))

    # save calls
    gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsExporter(
        shared_workspace, task.continuous_model, task.continuous_model_approx, sample_names,
        args.output_calls_path)()

    # save a copy of targets in the calls path
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_calls_path, gcnvkernel.io_consts.default_interval_list_filename))
