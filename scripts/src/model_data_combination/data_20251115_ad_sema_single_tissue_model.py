from .data_20251115_ad_sema import (
    DataSource, common_data_loader, model_loader, ModelList, common_model_constructor,
    Keywords, MFARunningMode, it, ModelKeyword,
    UserDefinedMFASettings as OldUserDefinedMFASettings,
    result_process_information_dict_analysis, flux_process_map_func,
    result_process_figure_config_dict_generator, tissue_mid_name_list_dict_constructor,
    tissue_flux_name_list_dict_constructor, metabolic_network_parameter_generator,
)
from ..mfa_analysis.mfa_config_generator import ModelDataSpecificSettings


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    loss_percentile = 0.005
    optimization_num_for_each_suffix_dict = {
        None: 5000,
    }
    average_optimized_results = None

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
def return_result_label_model_data_config_dict(
        running_mode, experiment_name=None, dataset_id=None, suffix=None):
    model_tissue_list = [
        ModelKeyword.liver,
        # ModelKeyword.kidney,
        # ModelKeyword.cortex,
        ModelKeyword.hippocampus,
        # ModelKeyword.hypothalamus,
        # ModelKeyword.q_muscle,
        # ModelKeyword.soleus_muscle,
        # ModelKeyword.brown_adipose,
        # ModelKeyword.colon,
    ]
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
        each_data_target_mfa_data_name_list_dict = {
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
        if dataset_id is None:
            target_mfa_data_name_list = list(it.chain(*each_data_target_mfa_data_name_list_dict.values()))
        else:
            try:
                target_mfa_data_name_list = each_data_target_mfa_data_name_list_dict[dataset_id]
            except KeyError:
                raise KeyError(
                    f'The dataset_id {dataset_id} is unrecognized. Feasible dataset_id '
                    f'are {list(each_data_target_mfa_data_name_list_dict.keys())}')

        if len(target_mfa_data_name_list) == 0:
            target_mfa_data_dict = mfa_data_dict
            target_result_information_dict = result_information_dict
        else:
            target_mfa_data_dict = {
                mfa_data_name: mfa_data_dict[mfa_data_name]
                for mfa_data_name in target_mfa_data_name_list
            }
            target_result_information_dict = {
                mfa_data_name: result_information_dict[mfa_data_name]
                for mfa_data_name in target_mfa_data_name_list
            }

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
                final_result_information_dict[result_label] = target_result_information_dict[data_result_label]
        returned_content = (
            result_label_model_data_config_dict, experimental_result_information_dict_analysis,
            experimental_tissue_mid_name_list_dict_constructor,
            final_result_information_dict, experimental_figure_config_dict,
        )

    return returned_content


def return_data_content(running_mode, dataset_id=None, **kwargs):
    data_content, _ = common_data_loader(data_name, test_mode=False, natural_anti_correction=False)
    if running_mode == MFARunningMode.raw_experimental_data_plotting:
        return data_content
    else:
        (mfa_data_dict, data_content_1, data_content_2, mfa_information_dict, data_content_4) = data_content
        each_data_target_mfa_data_name_list_dict = {
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
        if dataset_id is None:
            target_mfa_data_name_list = list(it.chain(*each_data_target_mfa_data_name_list_dict.values()))
        else:
            try:
                target_mfa_data_name_list = each_data_target_mfa_data_name_list_dict[dataset_id]
            except KeyError:
                raise KeyError(
                    f'The dataset_id {dataset_id} is unrecognized. Feasible dataset_id '
                    f'are {list(each_data_target_mfa_data_name_list_dict.keys())}')
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
        for mfa_data_name, mfa_data_obj in mfa_data_dict.items():
            mfa_data_obj.add_input_metabolite_obj_dict(labeled_buffer_input_metabolite_obj_dict)
        new_data_content = (new_mfa_data_dict, data_content_1, data_content_2, new_mfa_information_dict, data_content_4)
        return new_data_content

