import os

# set theano flags
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore"

import logging
import argparse
import gcnvkernel


class GCNVHelpFormatter(argparse.HelpFormatter):

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

    def _get_default_metavar_for_optional(self, action):
        return action.type.__name__

    def _get_default_metavar_for_positional(self, action):
        return action.type.__name__


parser = argparse.ArgumentParser(description="gCNV case-sample calling tool",
                                 formatter_class=GCNVHelpFormatter)

# logging
parser.add_argument("--console_log_level",
                    type=str,
                    choices=["INFO", "WARNING", "DEBUG"],
                    default="INFO",
                    help="Console logging verbosity level")

parser.add_argument("--logfile_log_level",
                    type=str,
                    choices=["INFO", "WARNING", "DEBUG"],
                    default="DEBUG",
                    help="Logfile logging verbosity level")

parser.add_argument("--logfile",
                    type=str,
                    required=False,
                    default=argparse.SUPPRESS,
                    help="If provided, the output log will be written to file as well")

log_level_map = {
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "DEBUG": logging.DEBUG
}

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--targets_interval_list",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="Full interval list, possibly including extra annotations (in .tsv format)")

group.add_argument("--sample_read_counts_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Sample name and read counts table (in .tsv format)")

group.add_argument("--sample_ploidy_metadata_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Contig ploidy metadata of all samples (in .tsv format)")

group.add_argument("--sample_read_depth_metadata_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Read depth metadata of all samples (in .tsv format)")

group.add_argument("--first_target_index",
                   type=int,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="First target index to include in the analysis (inclusive)")

group.add_argument("--last_target_index",
                   type=int,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Last target index to include in the analysis (exclusive)")

group.add_argument("--input_model_path",
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
gcnvkernel.DenoisingModelConfig.expose_args(
    parser,
    hide={
        "--max_bias_factors",
        "--psi_scale",
        "--log_mean_bias_std",
        "--init_ard_rel_unexplained_variance"})

# add calling config args
gcnvkernel.CopyNumberCallingConfig.expose_args(
    parser,
    hide={
        '--p_flat',
        '--class_coherence_length',
        '--initialize_to_flat_class'
    })

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

    # file logger
    logging.basicConfig(level=log_level_map[args.logfile_log_level],
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=args.logfile if hasattr(args, 'logfile') else '/dev/null',
                        filemode='w')

    # console logger
    console = logging.StreamHandler()
    console.setLevel(log_level_map[args.console_log_level])
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    # load targets interval list if provided
    if 'targets_interval_list' in args:
        loaded_targets_interval_list = gcnvkernel.io.load_interval_list_tsv_file(args.targets_interval_list)
    else:
        loaded_targets_interval_list = None

    # load sample names, truncated counts, and interval list from the sample read counts table
    sample_names, n_st, targets_interval_list = gcnvkernel.io.load_truncated_counts(
        args.sample_read_counts_table, args.first_target_index, args.last_target_index,
        assert_targets_are_equal=False)

    # if an external interval list is provided, copy annotations to targets_interval_list
    if loaded_targets_interval_list is not None:
        gcnvkernel.inherit_interval_annotations(loaded_targets_interval_list, targets_interval_list)

    # load read depth and ploidy metadata
    sample_metadata_collection = gcnvkernel.SampleMetadataCollection()
    sample_metadata_collection.read_sample_read_depth_metadata(args.sample_read_depth_metadata_table)
    sample_metadata_collection.read_sample_ploidy_metadata(args.sample_ploidy_metadata_table)

    # setup sample contig ploidy array
    contigs_set = {target.contig for target in targets_interval_list}
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
        denoising_config, calling_config, targets_interval_list,
        n_st, baseline_copy_number_s, read_depth_s)

    initial_params_supplier = gcnvkernel.DefaultDenoisingModelInitializer(
        denoising_config, calling_config, shared_workspace)

    task = gcnvkernel.CaseSampleCallingTask(
        denoising_config, calling_config, inference_params,
        shared_workspace, initial_params_supplier, args.input_model_path)

    # go!
    task.engage()
    task.disengage()

    # save calls
    gcnvkernel.io.SamplePosteriorsExporter(
        shared_workspace, task.continuous_model, task.continuous_model_approx, sample_names,
        args.output_calls_path)()
