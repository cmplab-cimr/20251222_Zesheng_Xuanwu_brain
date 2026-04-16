from ..core.common.classes import OptionDict, MFAConfig, TransformDict
from ..core.common.functions import (
    tissue_specific_name_constructor, check_if_mix_flux, mid_name_process, tissue_name_breakdown)
from ..core.common.config import ParamName, ModelKeyword, CoreConstants
from ..common_parallel_solver.config import Keywords as ParallelSolverKeywords
from ..common_parallel_solver.result_class import FinalResult as GeneralFinalResult
from ..io_functions.config import check_and_mkdir_of_direct
from ..io_functions.result_processing_functions import (
    multiple_repeat_result_processing, loss_data_distribution_plotting,
    time_distribution_plotting, figure_plotting,
    experimental_data_plotting_func_template, experimental_mid_prediction_plotting_func_template,
    flux_comparison_plotting_func_template
)
from ..io_functions.result_output_functions import (
    output_raw_flux_data, output_predicted_mid_data, output_flux_comparison_data, output_raw_experimental_mid_data)

from ..common.config import (
    Direct as CommonDirect, Keywords as GeneralKeywords, Color, FigureData, FigureDataKeywords,
)
from ..common.packages import mp, np, os
from ..common.functions import (
    update_parameter_object, reversible_flux_pair_dict_generator,
    tissue_specific_reversible_flux_title_constructor_generator,)
from scripts.model.model_loader import ModelList, model_loader
from scripts.data.data_loader import DataSource, common_data_loader

from ..inventory import MFARunningMode, ExperimentName, default_data_model_name
from ..model_data_combination.common_data_model_loader import (
    common_data_model_function_loader)


class Direct(object):
    name = 'mfa_analysis'
    output_direct = f'{CommonDirect.output_direct}/{name}'
    common_data_direct = f'{CommonDirect.common_submitted_raw_data_direct}/{name}'
    experimental_result_display = 'experimental_result_display'


class Keywords(GeneralKeywords):
    pass


# test_mode = True
# running_mode = MFARunningMode.flux_analysis
# running_mode = MFARunningMode.result_process
running_mode = MFARunningMode.raw_experimental_data_plotting
