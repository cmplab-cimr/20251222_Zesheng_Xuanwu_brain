class DataSource(object):
    data_20251115_ad_sema = 'data_20251115_ad_sema'
    data_20260412_metabolomics = 'data_20260412_metabolomics'
    metabolomics_behavioral = 'metabolomics_behavioral'
    flux_behavioral = 'flux_behavioral'


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
    from .data_20251115_ad_sema.specific_data_parameters import Keyword
    if data_type == DataSource.data_20251115_ad_sema:
        from .data_20251115_ad_sema.specific_data_parameters import SpecificParameters
        content = (specific_data_dict_loader(SpecificParameters(), **kwargs)
                   if not keyword_only else None)
        return content, Keyword()
    if data_type == DataSource.data_20260412_metabolomics:
        from .data_20260412_metabolomics.metabolomic_data_loader import MetabolomicDataset
        content = (None if keyword_only
                   else MetabolomicDataset(min_msi_level=kwargs.get('min_msi_level', '3')))
        return content, Keyword()
    if data_type == DataSource.metabolomics_behavioral:
        from .metabolomics_flux_behaviroal_data_loader.behavioral_data_loader import load_behavioral_data
        content = None if keyword_only else load_behavioral_data('metabolomics')
        return content, Keyword()
    if data_type == DataSource.flux_behavioral:
        from .metabolomics_flux_behaviroal_data_loader.behavioral_data_loader import load_behavioral_data
        content = None if keyword_only else load_behavioral_data('flux')
        return content, Keyword()
    raise ValueError(f'未知 data_type: {data_type}')
