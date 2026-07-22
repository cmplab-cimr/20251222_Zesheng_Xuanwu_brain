from .data_20251115_ad_sema import (
    DataSource, common_data_loader, model_loader, ModelList, common_model_constructor,
    Keywords, MFARunningMode, it, ModelKeyword,
    UserDefinedMFASettings as OldUserDefinedMFASettings, MFASettings,
    result_process_information_dict_analysis, flux_process_map_func,
    result_process_figure_config_dict_generator,
    display_flux_name_dict as old_display_flux_name_dict,
    metabolic_network_parameter_generator, dataset_extractor,
)
from ..mfa_analysis.mfa_config_generator import ModelDataSpecificSettings


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    loss_percentile = 0.005
    optimization_num_for_each_suffix_dict = {
        None: 5000,
    }
    average_optimized_results = None
    user_defined_constant_flux_dict = {
        **MFASettings.user_defined_constant_flux_dict,
    }
    user_defined_specific_flux_range_dict = {
        **MFASettings.user_defined_specific_flux_range_dict,
        'Salvage_c': (None, (1, 10)),
    }

user_defined_mfa_settings = UserDefinedMFASettings()


class UserDefinedModelDataSpecificSettings(ModelDataSpecificSettings):
    def __init__(self, model_tissue, data_content_parameter_dict):
        user_defined_constant_flux_dict = {
            **user_defined_mfa_settings.user_defined_constant_flux_dict,
        }
        user_defined_specific_flux_range_dict = {
            **user_defined_mfa_settings.user_defined_specific_flux_range_dict,
        }
        super().__init__(
            user_defined_specific_flux_range_dict=user_defined_specific_flux_range_dict,
            user_defined_constant_flux_dict=user_defined_constant_flux_dict,
            data_content_parameter_dict=data_content_parameter_dict)


data_name = DataSource.data_20251115_ad_sema
model_name = ModelList.tissue_specific_model
_, keyword = common_data_loader(data_name, keyword_only=True)
model_tissue_list = [
    ModelKeyword.liver,
    # ModelKeyword.kidney,
    ModelKeyword.cortex,
    ModelKeyword.hippocampus,
    ModelKeyword.hypothalamus,
    # ModelKeyword.q_muscle,
    # ModelKeyword.soleus_muscle,
    # ModelKeyword.brown_adipose,
    # ModelKeyword.colon,
]
def return_result_label_model_data_config_dict(
        running_mode, experiment_name=None, dataset_id=None, suffix=None):
    (
        mfa_data_dict, experimental_result_information_dict_analysis,
        experimental_tissue_mid_name_list_dict_constructor,
        result_information_dict, experimental_figure_config_dict,
    ), _ = common_data_loader(data_name, test_mode=False, natural_anti_correction=False)
    if running_mode == MFARunningMode.raw_experimental_data_plotting:
        default_model = model_tissue_list[0]
        user_defined_model, labeled_buffer_input_metabolite_obj_dict = model_loader(
            model_name, {Keywords.tissue: default_model})
        mfa_model_obj = common_model_constructor(user_defined_model)

        # 构建 config
        model_data_specific_setting_obj = UserDefinedModelDataSpecificSettings(
            default_model, user_defined_mfa_settings.data_content_parameter_dict)
        mfa_config_parameter_dict = model_data_specific_setting_obj.generate_mfa_config_parameter_dict(
            user_defined_model, mfa_model_obj)
        slsqp_mfa_config = user_defined_mfa_settings.construct_mfa_config(mfa_config_parameter_dict)

        result_label_model_data_config_dict = {data_result_label:
            (
                None, mfa_model_obj,
                data_result_label, mfa_data_obj,
                None, slsqp_mfa_config,
            ) for data_result_label, mfa_data_obj in mfa_data_dict.items()
        }
        returned_content = (
            result_label_model_data_config_dict, experimental_result_information_dict_analysis,
            experimental_tissue_mid_name_list_dict_constructor,
            result_information_dict, experimental_figure_config_dict,
        )
    else:
        complete_dataset_dict = {
            'part1': [
                'glucose__ad__ctrl__1',
                'glucose__ad__ctrl__2',
                'glucose__ad__ctrl__3',
                'glucose__ad__semaglutide__1',
                'glucose__ad__semaglutide__2',
                'glucose__ad__semaglutide__3',
            ],
            'part2': [
                'glucose__ad__ctrl__4',
                'glucose__ad__ctrl__5',
                'glucose__ad__ctrl__6',
                'glucose__ad__semaglutide__4',
                'glucose__ad__semaglutide__5',
                'glucose__ad__semaglutide__6',
            ]
        }
        target_mfa_data_dict, target_mfa_information_dict = dataset_extractor(
            mfa_data_dict, result_information_dict, dataset_id, complete_dataset_dict)

        result_label_model_data_config_dict = {}
        final_result_information_dict = {}
        # 构建 model
        for model_tissue in model_tissue_list:
            user_defined_model, labeled_buffer_input_metabolite_obj_dict = model_loader(
                model_name, {Keywords.tissue: model_tissue})
            mfa_model_obj = common_model_constructor(user_defined_model)

            # 构建 config
            model_data_specific_setting_obj = UserDefinedModelDataSpecificSettings(
                model_tissue, user_defined_mfa_settings.data_content_parameter_dict)
            mfa_config_parameter_dict = model_data_specific_setting_obj.generate_mfa_config_parameter_dict(
                user_defined_model, mfa_model_obj)
            slsqp_mfa_config = user_defined_mfa_settings.construct_mfa_config(mfa_config_parameter_dict)

            # 构建 result_label → (model_label, user_defined_model, mfa_model_obj, data_label, mfa_data_obj, config_label, slsqp_mfa_config)
            for data_result_label, mfa_data_obj in target_mfa_data_dict.items():
                result_label = f'{data_result_label}__{model_tissue}'
                mfa_data_obj.add_input_metabolite_obj_dict(labeled_buffer_input_metabolite_obj_dict)
                result_label_model_data_config_dict[result_label] = (
                    model_tissue, mfa_model_obj,
                    data_result_label, mfa_data_obj,
                    model_tissue, slsqp_mfa_config,
                )
                final_result_information_dict[result_label] = target_mfa_information_dict[data_result_label]
        returned_content = (
            result_label_model_data_config_dict, experimental_result_information_dict_analysis,
            experimental_tissue_mid_name_list_dict_constructor,
            final_result_information_dict, experimental_figure_config_dict,
        )

    return returned_content


display_flux_name_dict = {
    **old_display_flux_name_dict,
    ('LAC_c_isf_exchange', 'LAC_c_isf_exchange__R'): 'LAC isf output net',
    ('PYR_c_isf_exchange', 'PYR_c_isf_exchange__R'): 'PYR isf output net',
    ('ALA_c_isf_exchange', 'ALA_c_isf_exchange__R'): 'ALA isf output net',
    'GLU_c_isf_exchange': 'GLU isf output net',
    'GLN_c_isf_exchange__R': 'GLN isf input net',
}


def tissue_mapper(current_experiment_name):
    for specific_tissue in model_tissue_list:
        if current_experiment_name.endswith(specific_tissue):
            return specific_tissue
    raise ValueError()


def tissue_mid_name_list_dict_constructor(current_experiment_name=model_name, **kwargs):
    """单组织模型 MID 展示: serum + liver + 脑组织(cortex/hippocampus/hypothalamus)。"""
    mid_name_list_in_serum = [
        ['glucose', 'pyruvate', 'lactate', 'alanine', 'serine'],
        ['glutamate', 'glutamine', 'aspartate', '3-hydroxybutyrate'],
    ]
    mid_name_list_in_one_tissue = [
        ['glucose', 'pyruvate', 'lactate', 'alanine', 'acetyl-coa', '3-hydroxybutyrate'],
        ['citrate', 'glutamate', 'glutamine', 'succinate', 'fumarate', 'malate', 'aspartate'],
    ]
    tissue_specific_mid_name_list_dict = {
        each_model_name: mid_name_list_in_one_tissue
        for each_model_name in model_tissue_list
    }
    tissue_mid_name_list_dict = {
        ModelKeyword.serum: mid_name_list_in_serum,
        **tissue_specific_mid_name_list_dict,
    }
    return tissue_mid_name_list_dict, keyword.display_mid_name_dict, keyword.display_tissue_name_dict


def tissue_flux_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    """单组织模型重要 flux 展示: serum/liver/脑组织(糖酵解 + TCA + 神经递质交换)。"""
    flux_name_list_in_serum = [
        ['GLYCOLYSIS_c', ('LDH_e', 'LDH_e__R')],
    ]
    flux_name_list_in_liver = [
        [('GLC_exchange', 'GLC_exchange__R'), ('FBA_c', 'FBA_c__R'), 'PYK_c', ('LDH_c', 'LDH_c__R'),
         ('LAC_exchange', 'LAC_exchange__R'), 'BHB_output'],
        ['CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), ('ALA_exchange', 'ALA_exchange__R')],
    ]
    # 脑组织 (cortex/hippocampus/hypothalamus): 糖酵解 + TCA + 神经递质/乳酸交换
    # single_tissue 单区室模型用 *_exchange (非 neuron_astrocyte 的 *_c_isf_exchange)
    flux_name_list_in_brain = [
        [('GLC_exchange', 'GLC_exchange__R'), ('FBA_c', 'FBA_c__R'), 'PYK_c', ('LDH_c', 'LDH_c__R'),
         ('LAC_exchange', 'LAC_exchange__R'), 'PDH_m'],
        ['CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), ('ALA_exchange', 'ALA_exchange__R'),
         'GLU_exchange', 'GLN_exchange__R'],
    ]

    # current_tissue_name = tissue_mapper(current_model_name)
    tissue_specific_flux_name_list_dict = {
        ModelKeyword.liver: flux_name_list_in_liver,
        ModelKeyword.cortex: flux_name_list_in_brain,
        ModelKeyword.hippocampus: flux_name_list_in_brain,
        ModelKeyword.hypothalamus: flux_name_list_in_brain,
    }
    all_tissue_flux_name_list_dict = {
        ModelKeyword.serum: flux_name_list_in_serum,
        **tissue_specific_flux_name_list_dict,
    }
    specific_tissue_reversible_flux_pair_dict = {}
    return (
        all_tissue_flux_name_list_dict, display_flux_name_dict, keyword.display_tissue_name_dict,
        specific_tissue_reversible_flux_pair_dict, flux_process_map_func)

