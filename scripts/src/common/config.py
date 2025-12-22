from scripts.config import (
    Direct as GeneralDirect, np, RawBasicFigureData, RawFigureData, FigureDataKeywords,)


class Keywords(object):
    flux_raw_data = 'flux_raw_data'
    mid_raw_data = 'mid_raw_data'
    solver_descriptions = 'solver_descriptions'
    model_metabolites_reactions_standard_name = 'model_metabolites_reactions_standard_name'

    fasted = 'fast'
    fed = 'fed'

    optimized = 'optimized'
    unoptimized = 'unoptimized'
    experimental = 'experimental'

    mid_prediction = 'mid_prediction'
    flux_comparison = 'flux_comparison'
    glucose__ad__ctrl__1 = 'glucose__ad__ctrl__1'
    glucose__ad__semaglutide__1 = 'glucose__ad__semaglutide__1'
    glucose__ad__semaglutide__2 = 'glucose__ad__semaglutide__2'
    glucose__ad__semaglutide__3 = 'glucose__ad__semaglutide__3'



class DataType(object):
    test = 'test'


def rgba_to_rgb(raw_rgb, alpha, background=None):
    if background is None:
        background = np.array([1, 1, 1])
    return raw_rgb * alpha + background * (1 - alpha)


class Color(object):
    white = np.array([1, 1, 1])
    blue = np.array([21, 113, 177]) / 255
    orange = np.array([251, 138, 68]) / 255
    purple = np.array([112, 48, 160]) / 255
    light_blue = np.array([221, 241, 255]) / 255
    green = np.array([26, 150, 82]) / 255
    pink = np.array([213, 43, 132]) / 255

    alpha_value = 0.3
    alpha_for_bar_plot = alpha_value + 0.1
    alpha_for_heatmap = alpha_value + 0.2

    color_list = [
        rgba_to_rgb(blue, alpha_for_heatmap, white), white,
        rgba_to_rgb(orange, alpha_for_heatmap, white)]


class Direct(GeneralDirect):
    common_submitted_raw_data_direct = GeneralDirect.common_submitted_raw_data_direct
    data_direct = 'scripts/data'
    output_direct = 'scripts/output'
    raw_flux_analysis = 'raw_flux_analysis'
    flux_comparison = 'flux_comparison'
    predicted_experimental_mid_comparison_direct = 'predicted_experimental_mid_comparison'
    raw_and_mid_experimental_data_display_direct = 'raw_and_mid_experimental_data_display'
    metabolic_network_visualization_direct = 'metabolic_network_visualization'

    solution_array = 'solution_array'
    time_array = 'time_array'
    loss_array = 'loss_array'
    result_information = 'result_information'
    predicted_dict = 'predicted_dict'
    flux_name_index_dict = 'flux_name_index_dict'
    experimental_data = 'experimental_data'
    solution_id_array = 'solution_id_array'

    simulated_input_file_name = 'simulated_flux_vector_and_mid_data.py'
    simulated_input_file_path = f'scripts/src/simulated_data/{simulated_input_file_name}'
    simulated_output_py_file_direct = 'scripts/data/simulated_data'
    simulated_data_direct_name = 'simulated_data'
    simulated_output_xlsx_file_direct = f'{common_submitted_raw_data_direct}/{simulated_data_direct_name}'
    simulated_output_pickle_direct = simulated_output_py_file_direct


class FigureData(RawFigureData):
    def __init__(self, data_prefix, data_name):
        super().__init__(GeneralDirect.figure_raw_data_direct, data_prefix, data_name)


def check_and_mkdir_of_direct(direct_str, file_path=False):
    import pathlib
    direct_obj = pathlib.Path(direct_str)
    if file_path:
        direct_obj = direct_obj.parent
    dir_stack = []
    while not direct_obj.exists():
        dir_stack.append(direct_obj)
        direct_obj = direct_obj.parent
    while len(dir_stack) != 0:
        missed_direct = dir_stack.pop()
        missed_direct.mkdir()

