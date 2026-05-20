from .data_20251115_ad_sema import (
    DataSource, common_data_loader, model_loader, ModelList, common_model_constructor,
    Keywords, MFARunningMode, it, ModelKeyword,
    UserDefinedMFASettings as OldUserDefinedMFASettings,
    result_process_information_dict_analysis, flux_process_map_func,
    result_process_figure_config_dict_generator, tissue_mid_name_list_dict_constructor,
    tissue_flux_name_list_dict_constructor, metabolic_network_parameter_generator,
)

data_name = DataSource.data_20251115_ad_sema
model_name = ModelList.liver_cortex_compartmentalized_model_with_serum_supplement
_, keyword = common_data_loader(
    data_name, keyword_only=True, test_mode=False, natural_anti_correction=False)
user_defined_model, labeled_buffer_input_metabolite_obj_dict = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)
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
            ],
            'part3': [
                'glucose__ad__semaglutide__4',
                'glucose__ad__semaglutide__5',
                'glucose__ad__semaglutide__6',
            ]
        }
        if dataset_id is None:
            target_mfa_data_name_list = []
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
        for mfa_data_obj in new_mfa_data_dict.values():
            mfa_data_obj.add_input_metabolite_obj_dict(labeled_buffer_input_metabolite_obj_dict)
        new_data_content = (new_mfa_data_dict, data_content_1, data_content_2, new_mfa_information_dict, data_content_4)
        return new_data_content


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    optimization_num_for_each_suffix_dict = {
        None: 10000,
    }
    average_optimized_results = None

user_defined_mfa_settings = UserDefinedMFASettings()
