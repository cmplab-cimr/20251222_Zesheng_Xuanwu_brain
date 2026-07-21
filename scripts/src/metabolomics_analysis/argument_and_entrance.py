from ..inventory import MetabolomicsExperimentName, MetabolomicsRunningMode, metabolomics_data_model_comment
from ..common.packages import argparse


def arg_setting(computation_subparsers, experimental_operation_label):
    def run_analysis(args):
        main(metabolomics_parser, args)

    metabolomics_parser = computation_subparsers.add_parser(
        experimental_operation_label, prog='Metabolomics analysis',
        description='Run metabolomics data analysis',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Definition of data_model_name:\n\n{}'.format(
            '\n'.join([
                f'{enum_item.name:<50}{metabolomics_data_model_comment[enum_item]}'
                for enum_item in MetabolomicsExperimentName
            ])
        )
    )
    metabolomics_parser.add_argument(
        '-s', '--suffix', default=None, help='Suffix of stored folders')
    metabolomics_parser.add_argument(
        '-t', '--tissue', default=None,
        help='Specific tissue to analyze (e.g. serum, cortex, hippocampus). If not provided, all tissues are analyzed.')
    metabolomics_parser.add_argument(
        '--method', default='ttest', choices=['ttest', 'mannwhitney'],
        help='Statistical test method (default: ttest)')
    metabolomics_parser.add_argument(
        '--p_threshold', type=float, default=0.1,
        help='Raw p-value threshold (default: 0.1)')
    metabolomics_parser.add_argument(
        '--fc_threshold', type=float, default=1.5,
        help='Fold change threshold (default: 1.5)')
    running_mode_display = '{}'.format(
        ',  '.join([mode.value for mode in MetabolomicsRunningMode]))
    metabolomics_parser.add_argument(
        'running_mode', nargs='?', type=MetabolomicsRunningMode,
        choices=list(MetabolomicsRunningMode),
        help='Running mode of metabolomics analysis', default=None,
        metavar=running_mode_display)
    data_source_display = '{}'.format(
        ',  '.join([ds.value for ds in MetabolomicsExperimentName]))
    metabolomics_parser.add_argument(
        'experiment_name', nargs='?', type=MetabolomicsExperimentName,
        choices=list(MetabolomicsExperimentName),
        help='The name of experiment dataset.', default=None,
        metavar=data_source_display)
    metabolomics_parser.add_argument(
        '--n_permutations', type=int, default=200,
        help='Number of permutations for OPLS-DA (default: 200, 0 to skip)')
    metabolomics_parser.add_argument(
        '--vip_threshold', type=float, default=1.0,
        help='VIP score threshold for clustered heatmap (default: 1.0)')
    metabolomics_parser.add_argument(
        '--msi_level', default='3',
        choices=['1a', '1b', '2', '3'],
        help='Min MSI identification level: keep metabolites <= this level (default: 3 = keep all)')
    metabolomics_parser.set_defaults(func=run_analysis)


def main(metabolomics_parser, args):
    running_mode = args.running_mode
    if running_mode is None:
        metabolomics_parser.print_help()
    else:
        from .metabolomics_main_functions import metabolomics_analysis_main
        metabolomics_analysis_main(
            running_mode=running_mode, experiment_name=args.experiment_name,
            suffix=args.suffix, tissue=args.tissue,
            method=args.method, p_threshold=args.p_threshold,
            fc_threshold=args.fc_threshold,
            n_permutations=args.n_permutations,
            vip_threshold=args.vip_threshold,
            msi_level=args.msi_level)
