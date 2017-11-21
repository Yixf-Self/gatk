import argparse
import logging
import sys


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


def add_logging_args_to_argparse(parser: argparse.ArgumentParser):
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


def set_logging_config_from_args(args):
    # file logger
    logging.basicConfig(level=log_level_map[args.logfile_log_level],
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=args.logfile if hasattr(args, 'logfile') else '/dev/null',
                        filemode='w')

    # console logger
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(log_level_map[args.console_log_level])
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
