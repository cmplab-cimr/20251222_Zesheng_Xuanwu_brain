from .data_20251115_ad_sema import (
    DataSource, common_data_loader, model_loader, ModelList, common_model_constructor,
    Keywords, MFARunningMode, it, ModelKeyword,
    UserDefinedMFASettings as OldUserDefinedMFASettings, MFASettings,
    result_process_information_dict_analysis, flux_process_map_func,
    result_process_figure_config_dict_generator, tissue_mid_name_list_dict_constructor,
    tissue_flux_name_list_dict_constructor, metabolic_network_parameter_generator,
    dataset_extractor,
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
            ],
            'part3': [
                'glucose__ad__semaglutide__4',
                'glucose__ad__semaglutide__5',
                'glucose__ad__semaglutide__6',
            ]
        }
        new_mfa_data_dict, new_mfa_information_dict = dataset_extractor(
            mfa_data_dict, mfa_information_dict, dataset_id, complete_dataset_dict)
        for mfa_data_obj in new_mfa_data_dict.values():
            mfa_data_obj.add_input_metabolite_obj_dict(labeled_buffer_input_metabolite_obj_dict)
        new_data_content = (new_mfa_data_dict, data_content_1, data_content_2, new_mfa_information_dict, data_content_4)
        return new_data_content


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    optimization_num_for_each_suffix_dict = {
        None: 10000,
    }
    user_defined_constant_flux_dict = {
        **MFASettings.user_defined_constant_flux_dict,
        'GLC_infuse': (ModelKeyword.serum, 200)
    }
    user_defined_specific_flux_range_dict = {
        **MFASettings.user_defined_specific_flux_range_dict,
        'Salvage_c': (None, (1, 10)),
    }
    average_optimized_results = None

user_defined_mfa_settings = UserDefinedMFASettings()
