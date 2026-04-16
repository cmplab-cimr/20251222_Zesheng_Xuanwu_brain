class DataSource(object):
    data_20251115_ad_sema = 'data_20251115_ad_sema'


def specific_data_dict_loader(specific_parameter, reload=False, **kwargs):
    from .config import (
        pickle_save, pickle_load, check_file_existence,)

    tissue_labeling_pickle_file_path = specific_parameter.pickle_file_path

    data_dict_str = 'data_dict'
    info_dict_str = 'info_dict'

    if not check_file_existence(tissue_labeling_pickle_file_path) or reload:
        for sheet_name, current_data_parameter_dict in specific_parameter.return_data_parameter_dict().items():
            current_data_dict = specific_parameter.target_metabolite_data_loader(**current_data_parameter_dict)
            specific_parameter.add_data_sheet(sheet_name, current_data_dict)
        (
            final_target_mfa_data_dict, final_result_information_dict
        ) = specific_parameter.return_multi_tissue_mfa_data_dict()
        output_dict = {
            data_dict_str: final_target_mfa_data_dict,
            info_dict_str: final_result_information_dict,
        }
        pickle_save(output_dict, tissue_labeling_pickle_file_path)
    else:
        output_dict = pickle_load(tissue_labeling_pickle_file_path)
        final_target_mfa_data_dict = output_dict[data_dict_str]
        final_result_information_dict = output_dict[info_dict_str]
    result_information_dict_analysis = specific_parameter.result_information_dict_analysis
    tissue_mid_name_list_dict_constructor = specific_parameter.tissue_mid_name_list_dict_constructor
    experimental_figure_config_dict = specific_parameter.experimental_figure_config_dict_generator()
    return (
        final_target_mfa_data_dict, result_information_dict_analysis,
        tissue_mid_name_list_dict_constructor, final_result_information_dict, experimental_figure_config_dict)


def common_data_loader(data_type, keyword_only=False, **kwargs):
    if data_type == DataSource.data_20251115_ad_sema:
        from .data_20251115_ad_sema.specific_data_parameters import SpecificParameters, Keyword
    else:
        raise ValueError()
    if not keyword_only:
        specific_parameters = SpecificParameters()
        content = specific_data_dict_loader(specific_parameters, **kwargs)
    else:
        content = None
    keyword = Keyword()
    return content, keyword
