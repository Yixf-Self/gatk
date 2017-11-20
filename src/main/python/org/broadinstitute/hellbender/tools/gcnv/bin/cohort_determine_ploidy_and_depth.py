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


parser = argparse.ArgumentParser(description="gCNV contig ploidy and read depth determination tool",
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

group.add_argument("--interval_list",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Interval list of included genomic regions in the analysis (in .tsv format)")

group.add_argument("--sample_coverage_metadata",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Coverage metadata of all samples (in .tsv format)")

group.add_argument("--contig_ploidy_prior_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Contig ploidy prior probabilities (in .tsv format)")

group.add_argument("--output_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write the ploidy model for future single-sample ploidy determination use")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write posteriors")

# optional arguments
gcnvkernel.PloidyModelConfig.expose_args(parser)

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(
    parser,
    override_default={
        "--learning_rate": 0.1,
        "--adamax_beta2": 0.999,
        "--log_emission_samples_per_round": 1000,
        "--log_emission_sampling_rounds": 50,
        "--log_emission_sampling_median_rel_error": 1e-3,
        "--max_advi_iter_first_epoch": 1000,
        "--max_advi_iter_subsequent_epochs": 1000,
        "--convergence_snr_averaging_window": 5000,
        "--convergence_snr_countdown_window": 100,
        "--num_thermal_epochs": 10,
        "--max_calling_iters": 1,
        "--caller_update_convergence_threshold": 1e-3
    },
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

    # read contig ploidy prior map from file
    contig_ploidy_prior_map = gcnvkernel.PloidyModelConfig.get_contig_ploidy_prior_map_from_tsv_file(
        args.contig_ploidy_prior_table)

    # load targets interval list
    interval_list = gcnvkernel.io.load_interval_list_tsv_file(args.interval_list)

    # load sample coverage metadata
    sample_metadata_collection = gcnvkernel.SampleMetadataCollection()
    sample_names = sample_metadata_collection.read_sample_coverage_metadata(args.sample_coverage_metadata)

    # generate targets metadata
    intervals_metadata = gcnvkernel.TargetsIntervalListMetadata(interval_list)

    # inject ploidy prior map to the dictionary of parsed args
    args_dict = args.__dict__
    args_dict['contig_ploidy_prior_map'] = contig_ploidy_prior_map

    ploidy_config = gcnvkernel.PloidyModelConfig.from_args_dict(args_dict)
    ploidy_inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)
    ploidy_workspace = gcnvkernel.PloidyWorkspace(ploidy_config, intervals_metadata, sample_names,
                                                  sample_metadata_collection)
    ploidy_task = gcnvkernel.PloidyInferenceTask(ploidy_inference_params, ploidy_config, ploidy_workspace)

    # go!
    ploidy_task.engage()
    ploidy_task.disengage()

    # save model parameters
    gcnvkernel.io.PloidyModelExporter(ploidy_config, ploidy_workspace,
                                      ploidy_task.continuous_model, ploidy_task.continuous_model_approx,
                                      args.output_model_path)()

    # sample sample-specific posteriors
    gcnvkernel.io.SamplePloidyExporter(ploidy_config, ploidy_workspace,
                                       ploidy_task.continuous_model, ploidy_task.continuous_model_approx,
                                       args.output_calls_path)()
