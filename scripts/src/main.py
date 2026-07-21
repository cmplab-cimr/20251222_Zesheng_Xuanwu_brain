import argparse


class ComputationFunction(object):
    experimental_mfa = 'experimental_mfa'
    metabolomics = 'metabolomics'
    machine_learning = 'machine_learning'


def arg_setting(subparsers):
    def print_computation_help(args):
        computation_parser.print_help()

    computation_parser = subparsers.add_parser('computation', help='Run computation functions')
    computation_subparsers = computation_parser.add_subparsers(
        title='Commands',
        description='Different content for computation',
        help='Decide to run different analysis functions')
    from scripts.src.mfa_analysis.argument_and_entrance import arg_setting as \
        experiments_arg_setting
    experiments_arg_setting(computation_subparsers, ComputationFunction.experimental_mfa)

    from scripts.src.metabolomics_analysis.argument_and_entrance import arg_setting as \
        metabolomics_arg_setting
    metabolomics_arg_setting(computation_subparsers, ComputationFunction.metabolomics)

    # from scripts.src.machine_learning.machine_learning_main import arg_setting as \
    #     ml_arg_setting
    # ml_arg_setting(computation_subparsers, ComputationFunction.machine_learning)

    computation_parser.set_defaults(func=print_computation_help)
