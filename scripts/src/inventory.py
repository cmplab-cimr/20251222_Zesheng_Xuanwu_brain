from .common.packages import ValueEnum


class ExperimentName(ValueEnum):
    example_data_20250815 = 'data_20251115_ad_sema'

default_data_model_name = ExperimentName.example_data_20250815


class MLRunningMode(ValueEnum):
    met_gat = 'met_gat'


class MFARunningMode(ValueEnum):
    flux_analysis = 'flux_analysis'
    raw_experimental_data_plotting = 'raw_experimental_data_plotting'
    result_process = 'result_process'
    solver_output = 'solver_output'
    flux_analysis_test = 'flux_analysis_test'


data_model_comment = {
    ExperimentName.example_data_20250815:
        'Data from mice',
}
