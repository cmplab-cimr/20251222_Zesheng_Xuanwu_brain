from ..inventory import ExperimentName, MFARunningMode, data_model_comment
from ..common.packages import argparse


def arg_setting(computation_subparsers, experimental_operation_label):
    def flux_analysis(args):
        main(flux_analysis_parser, args)

    flux_analysis_parser = computation_subparsers.add_parser(
        experimental_operation_label, prog='Analysis of experimental data',
        description='Run MFA for several experimental data analyses',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Definition of data_model_name:\n\n{}'.format(
            '\n'.join([
                f'{enum_item.name:<50}{data_model_comment[enum_item]}'
                for enum_item in ExperimentName
            ])
        )
    )
    flux_analysis_parser.add_argument(
        '-s', '--suffix', default=None, help='Suffix of stored folders')
    flux_analysis_parser.add_argument(
        '-i', '--dataset_id', default=None, help='Specific dataset id')
    flux_analysis_parser.add_argument(
        '-t', '--test_mode', action='store_true', default=False,
        help='Whether the code is executed in test mode, which means less sample number and shorter time.')
    flux_analysis_parser.add_argument(
        '-d', '--docker_mode', default=None,
        help='Whether the code is executed in docker mode, which means no parallel settings.')
    flux_analysis_parser.add_argument(
        '-p', '--parallel_num', type=int, default=None,
        help='Number of parallel processes. If not provided, it will be selected according to CPU cores.'
    )
    running_mode_display = '{}'.format(',  '.join([running_mode.value for running_mode in MFARunningMode]))
    flux_analysis_parser.add_argument(
        'running_mode', nargs='?', type=MFARunningMode, choices=list(MFARunningMode),
        help='Running mode of flux analysis', default=None, metavar=running_mode_display)
    data_source_display = '{}'.format(',  '.join([data_source.value for data_source in ExperimentName]))
    flux_analysis_parser.add_argument(
        'experiment_name', nargs='?', type=ExperimentName, choices=list(ExperimentName),
        help='The name of experiment you want to compute.', default=None,
        metavar=data_source_display)
    flux_analysis_parser.set_defaults(func=flux_analysis)


def main(experimental_data_analysis_parser, args):
    running_mode = args.running_mode
    if running_mode is None:
        experimental_data_analysis_parser.print_help()
    else:
        from .mfa_main_functions import multiple_tissue_model_analysis_main
        multiple_tissue_model_analysis_main(
            running_mode=running_mode, experiment_name=args.experiment_name, dataset_id=args.dataset_id,
            suffix=args.suffix, test_mode=args.test_mode, docker_mode=args.docker_mode,
            parallel_num=args.parallel_num)
