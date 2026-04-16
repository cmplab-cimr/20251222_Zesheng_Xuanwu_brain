from ..core.common.functions import tissue_specific_name_constructor

def tissue_specific_metabolite_mid_list_constructor(
        formatted_mean_data_dict, formatted_std_data_dict, data_len_dict, tissue_mid_name_list_dict):
    mid_name_list = []
    filtered_mean_data_dict = {}
    filtered_std_data_dict = {}
    for tissue_name, each_tissue_mid_name_list in tissue_mid_name_list_dict.items():
        for mid_name_row in each_tissue_mid_name_list:
            new_mid_name_row = []
            for metabolite_name in mid_name_row:
                new_tissue_specific_mid_name = tissue_specific_name_constructor(metabolite_name, tissue_name)
                try:
                    this_mid_data_dict = formatted_mean_data_dict[new_tissue_specific_mid_name]
                except KeyError:
                    this_mid_data_dict = None
                    try:
                        mid_data_len = data_len_dict[metabolite_name]
                    except KeyError:
                        mid_data_len = None
                    new_mid_name_row.append(mid_data_len)
                else:
                    new_mid_name_row.append(new_tissue_specific_mid_name)
                filtered_mean_data_dict[new_tissue_specific_mid_name] = this_mid_data_dict
                try:
                    this_mid_std_data_dict = formatted_std_data_dict[new_tissue_specific_mid_name]
                except KeyError:
                    this_mid_std_data_dict = None
                filtered_std_data_dict[new_tissue_specific_mid_name] = this_mid_std_data_dict
            mid_name_list.append(new_mid_name_row)
    return mid_name_list, filtered_mean_data_dict, filtered_std_data_dict


def tissue_specific_flux_list_constructor(all_tissue_flux_name_list_dict, display_name_dict):
    flux_name_list = []
    tissue_specific_display_flux_dict = {}
    for tissue_name, each_tissue_flux_name_list in all_tissue_flux_name_list_dict.items():
        for flux_name_row in each_tissue_flux_name_list:
            new_mid_name_row = []
            for flux_name in flux_name_row:
                if isinstance(flux_name, tuple):
                    flux_name_str = f'{flux_name[0]} - {flux_name[1]}'
                elif isinstance(flux_name, str):
                    flux_name_str = flux_name
                else:
                    raise ValueError()
                new_tissue_specific_flux_name = tissue_specific_name_constructor(flux_name_str, tissue_name)
                new_mid_name_row.append(new_tissue_specific_flux_name)
                if flux_name in display_name_dict:
                    display_flux_name_str = display_name_dict[flux_name]
                else:
                    display_flux_name_str = flux_name_str
                tissue_specific_display_flux_name = f'{tissue_name}: {display_flux_name_str}'
                tissue_specific_display_flux_dict[new_tissue_specific_flux_name] = tissue_specific_display_flux_name
            flux_name_list.append(new_mid_name_row)
    return flux_name_list, tissue_specific_display_flux_dict
