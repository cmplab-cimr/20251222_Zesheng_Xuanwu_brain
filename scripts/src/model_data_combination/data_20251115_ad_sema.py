from scripts.src.common.config import Color, Keywords, Direct
from scripts.src.common.functions import data_param_list_generator_func_template, \
    collect_results_func_template
from scripts.src.common.packages import it
from scripts.src.io_functions.result_processing_functions import common_flux_comparison_func, \
    experimental_data_plotting_func_template, mid_name_list_generator_for_multiple_labeling_substrate

from scripts.data.data_loader import DataSource, common_data_loader

from scripts.model.model_loader import model_loader, ModelList
from scripts.src.core.model.model_constructor import common_model_constructor
from scripts.src.core.common.config import ParamName, ModelKeyword
from scripts.src.core.common.classes import OptionDict, DefaultDict
from scripts.src.core.data.data_class import average_multiple_mfa_data

from ..inventory import MFARunningMode

from scripts.src.mfa_analysis.mfa_config_generator import MFASettings

average_str = 'average'

data_name = DataSource.data_20251115_ad_sema
model_name = ModelList.liver_cortex_compartmentalized_model
data_content, keyword = common_data_loader(data_name, test_mode=False, natural_anti_correction=False)
user_defined_model = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)
def return_data_content(running_mode, dataset_id=None, **kwargs):
    if running_mode == MFARunningMode.raw_experimental_data_plotting:
        return data_content
    else:
        (mfa_data_dict, data_content_1, data_content_2, mfa_information_dict, data_content_4) = data_content
        each_data_target_mfa_data_name_list_dict = {
            'ctrl_vs_sema': [
                'glucose__ad__ctrl__1',
                'glucose__ad__ctrl__2',
                'glucose__ad__ctrl__3',
                # 'glucose__ad__ctrl__4',
                # 'glucose__ad__ctrl__5',
                # 'glucose__ad__ctrl__6',
                'glucose__ad__semaglutide__1',
                'glucose__ad__semaglutide__2',
                'glucose__ad__semaglutide__3',
                # 'glucose__ad__semaglutide__4',
                # 'glucose__ad__semaglutide__5',
                # 'glucose__ad__semaglutide__6',
            ],
        }
        if dataset_id is None:
            target_mfa_data_name_list = list(it.chain(*each_data_target_mfa_data_name_list_dict.values()))
        else:
            target_mfa_data_name_list = each_data_target_mfa_data_name_list_dict[dataset_id]
        if len(target_mfa_data_name_list) == 0:
            new_mfa_data_dict = mfa_data_dict
            new_mfa_information_dict = mfa_information_dict
        else:
            new_mfa_data_dict = {
                mfa_data_name: mfa_data_dict[mfa_data_name]
                for mfa_data_name in target_mfa_data_name_list
            }
            new_mfa_information_dict = {
                mfa_data_name: mfa_information_dict[mfa_data_name]
                for mfa_data_name in target_mfa_data_name_list
            }
        new_data_content = (new_mfa_data_dict, data_content_1, data_content_2, new_mfa_information_dict, data_content_4)
        return new_data_content


class UserDefinedMFASettings(MFASettings):
    loss_percentile = 0.005
    user_defined_constant_flux_dict = {
        **MFASettings.user_defined_constant_flux_dict,
    }
    user_defined_specific_flux_range_dict = {
        **MFASettings.user_defined_specific_flux_range_dict,
        'Salvage_c': (None, (1, 10)),
    }
    optimization_num_for_each_suffix_dict = {
        None: 10000,
    }
    average_optimized_results = 50

user_defined_mfa_settings = UserDefinedMFASettings()

def result_process_information_dict_analysis(complete_result_information_dict, target=Keywords.mid_prediction):
    from scripts.figures.common.config import ColorConfig

    color_dict_without_index = {
        (keyword.no_exec, keyword.ctrl): ColorConfig.normal_blue,
        (keyword.no_exec, keyword.mr): ColorConfig.orange,
        (keyword.no_exec, keyword.tau): ColorConfig.light_bright_mixed_mid_red,
        (keyword.no_exec, keyword.mr_tau): ColorConfig.light_green,
        (keyword.with_exec, keyword.ctrl): ColorConfig.dark_blue,
        (keyword.with_exec, keyword.mr): ColorConfig.dark_orange,
        (keyword.with_exec, keyword.tau): ColorConfig.mid_red,
        (keyword.with_exec, keyword.mr_tau): ColorConfig.green,
    }
    color_order_dict = {}
    total_index = 0
    each_group_replicate_num = 3
    for (exercise_name, condition_name), color in color_dict_without_index.items():
        for current_index in range(1, 1 + each_group_replicate_num):
            color_order_dict[(exercise_name, condition_name, current_index)] = (total_index + current_index, color)
        total_index += each_group_replicate_num

    return_parameter_list = []
    if target == Keywords.mid_prediction:
        for result_label, current_result_information_dict in complete_result_information_dict.items():
            major_key = data_name
            # exercise = current_result_information_dict[keyword.exercise]
            # exer_str = '+' if exercise == keyword.with_exec else '-'
            exer_str = result_label[-1]
            exercise = keyword.with_exec if exer_str == '+' else keyword.no_exec
            condition = current_result_information_dict[keyword.condition]
            index = current_result_information_dict[keyword.index]
            minor_key_list = []
            minor_key_str = f'{condition}_{index}{exer_str}'
            if index == average_str:
                color_index = 1
            else:
                color_index = index
            order_index, current_color = color_order_dict[(exercise, condition, color_index)]
            return_parameter_list.append(
                (result_label, major_key, minor_key_list, minor_key_str, current_color, order_index))
    elif target == Keywords.flux_comparison:
        for result_label, current_result_information_dict in complete_result_information_dict.items():
            major_key = data_name
            # exercise = current_result_information_dict[keyword.exercise]
            # exer_str = '+' if exercise == keyword.with_exec else '-'
            exer_str = result_label[-1]
            exercise = keyword.with_exec if exer_str == '+' else keyword.no_exec
            condition = current_result_information_dict[keyword.condition]
            index = current_result_information_dict[keyword.index]
            minor_key_list = []
            minor_key_str = f'{condition}_{index}{exer_str}'
            if index == average_str:
                color_index = 1
            else:
                color_index = index
            order_index, current_color = color_order_dict[(exercise, condition, color_index)]
            return_parameter_list.append(
                (result_label, major_key, minor_key_list, minor_key_str, current_color, order_index))
    else:
        raise ValueError()
    return return_parameter_list


def result_process_figure_config_dict_generator():
    figure_config_dict = {
        'figure_size': (8, 18),
        # 'size': [0.8, 1.1],
        # 'scale': 0.7,
        'legend': True,
    }
    return {
        'predicted_mid_figure_config_dict': figure_config_dict,
        'flux_comparison_figure_config_dict': figure_config_dict,
    }

flux_name_list_in_one_tissue = [
    ['GLC_input', ('FBA_c', 'FBA_c__R'), 'PYK_c', 'LAC_output'],
    ['CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), 'ALA_input'],
]


display_flux_name_dict = {
    ('LDH_e', 'LDH_e__R'): 'LDH_e net',
    ('FBA_c', 'FBA_c__R'): 'FBA net',
    ('LAC_exchange', 'LAC_exchange__R'): 'LAC output net',
    ('MDH_c', 'MDH_c__R'): 'MDH net',
    ('FA_exchange', 'FA_exchange__R'): 'FA net',
    ('GLC_input', 'GLC_output'): 'GLC input net',
    ('GLU_input', 'GLU_output'): 'GLU input net',
    ('GLN_input', 'GLN_output'): 'GLN input net',
    ('ASP_input', 'ASP_output'): 'ASP input net',
    ('PYR_input', 'PYR_output'): 'PYR input net',
    ('ALA_input', 'ALA_output'): 'ALA input net',
}


def tissue_mid_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    mid_name_list_in_serum = [
        ['glucose', 'pyruvate', 'lactate', 'alanine', 'serine', ],
        ['glutamate', 'glutamine', 'aspartate',  ],
    ]
    mid_name_list_in_one_tissue = [
        [
            'glucose', 'pyruvate', 'lactate', 'alanine', 'acetyl-coa', ],
        [
            'citrate/isocitrate', 'glutamate', 'glutamine', 'succinate', 'fumarate',
            'malate', 'aspartate', ],'gut--citrate/isocitrate'
    ]

    if current_model_name == ModelList.complete_model_whole_body:
        tissue = keyword.whole_body
    elif current_model_name == ModelList.complete_model_gut:
        tissue = keyword.gut
    else:
        raise ValueError()
    base_tissue_mid_name_list_dict = {
        tissue: mid_name_list_in_one_tissue,
    }
    tissue_mid_name_list_dict = base_tissue_mid_name_list_dict
    return tissue_mid_name_list_dict, keyword.display_mid_name_dict, keyword.display_tissue_name_dict


def tissue_flux_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    if current_model_name == ModelList.complete_model_whole_body:
        tissue = keyword.whole_body
    elif current_model_name == ModelList.complete_model_gut:
        tissue = keyword.gut
    else:
        raise ValueError()
    all_tissue_flux_name_list_dict = {
        tissue: flux_name_list_in_one_tissue,
    }
    specific_tissue_reversible_flux_pair_dict = {
        tissue: [
            # ('GLC_input', 'GLC_output'),
            # ('GLU_input', 'GLU_output'),
            # ('GLN_input', 'GLN_output'),
            # ('ASP_input', 'ASP_output'),
            # ('PYR_input', 'PYR_output'),
            # ('SER_input', 'SER_output'),
            # ('ALA_input', 'ALA_output'),
        ],
    }
    return all_tissue_flux_name_list_dict, display_flux_name_dict, specific_tissue_reversible_flux_pair_dict


def metabolic_network_parameter_generator():
    experimental_mid_metabolite_set = {
        'AKG_c', 'AKG_m',
        'ALA_c',
        'ASP_c', 'ASP_m',
        'ATP_c',
        'CIT_c', 'CIT_m',
        'FBP_c',
        'GLC_c',
        'GLU_c', 'GLU_m',
        'GLN_c', 'GLN_m',
        'GLY_c',
        'LAC_c',
        'MAL_c', 'MAL_m',
        'PYR_c', 'PYR_m',
        'SAM_c',
        'SER_c',
        'SUC_m',
    }
    experimental_mixed_mid_metabolite_set = {
        'AKG_m', 'AKG_c',
        'ASP_m', 'ASP_c',
        'CIT_m', 'CIT_c',
        'GLU_m', 'GLU_c',
        'GLN_m', 'GLN_c',
        'MAL_m', 'MAL_c',
        'PYR_c', 'PYR_m',
    }
    biomass_metabolite_set = {
        'ALA_c', 'RIB5P_c', 'GLY_c', 'SER_c', 'ASP_c',
        'GLU_c', 'GLN_c', 'ATP_c', 'MET_c'
    }
    input_metabolite_set = {
        'GLC_e', 'GLC_unlabelled_e', 'GLN_e', 'GLU_e', 'ASP_e', 'SER_e', 'GLY_e', 'ALA_e', 'LAC_e', 'MET_e', 'PRP_e'
    }
    c13_labeling_metabolite_set = {
        # 'GLC_e',
        'PRP_e',
    }
    boundary_flux_set = {
        'GLC_input'
    }
    return experimental_mid_metabolite_set, experimental_mixed_mid_metabolite_set, biomass_metabolite_set, \
        input_metabolite_set, c13_labeling_metabolite_set, boundary_flux_set, False
