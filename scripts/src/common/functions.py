from scripts.src.core.common.config import CoreConstants
from scripts.src.core.common.functions import reverse_reaction_name, tissue_specific_name_constructor

from .packages import it


def mid_name_process(raw_mid_name):
    emu_sep = CoreConstants.emu_carbon_list_str_sep
    modified_str_list = []
    str_start = 0
    while str_start < len(raw_mid_name):
        emu_location = raw_mid_name.find(emu_sep, str_start)
        if emu_location == -1:
            modified_str_list.append(raw_mid_name[str_start:])
            break
        modified_str_list.append(raw_mid_name[str_start:emu_location])
        new_start = emu_location + len(emu_sep)
        while new_start != len(raw_mid_name) and raw_mid_name[new_start] == '1':
            new_start += 1
        str_start = new_start
    modified_str = ''.join(modified_str_list)
    return modified_str


def default_parameter_extract(
        option_dict: dict, key, default_value=None, force=False, pop=False, repeat_default_value=False):
    def single_extract(_option_dict, _key, _default_value):
        if force or _key in _option_dict:
            if pop:
                _value = _option_dict.pop(_key)
            else:
                _value = _option_dict[_key]
            return _value
        else:
            return _default_value

    if isinstance(key, str):
        return single_extract(option_dict, key, default_value)
    elif isinstance(key, (list, tuple)):
        result_list = []
        if isinstance(default_value, (list, tuple)):
            default_value_iter = default_value
        elif force or repeat_default_value:
            default_value_iter = it.repeat(default_value)
        else:
            raise ValueError()
        for each_key, each_default_value in zip(key, default_value_iter):
            result_list.append(single_extract(option_dict, each_key, each_default_value))
        return result_list
    else:
        raise ValueError()


def update_parameter_object(original_parameter_object, new_parameter_object):
    for item_key, item_value in new_parameter_object.__class__.__dict__.items():
        if not item_key.startswith('__'):
            if hasattr(original_parameter_object, item_key) and isinstance(
                    getattr(original_parameter_object, item_key), dict):
                getattr(original_parameter_object, item_key).update(item_value)
            else:
                original_parameter_object.__setattr__(item_key, item_value)
    return original_parameter_object


def data_param_list_generator_func_template(keyword_list, extra_key_default_value_dict=None):
    empty_str = ''
    miscellaneous_str = 'miscellaneous'

    def data_param_list_generator_iter_func(current_dict_list, current_index, current_complete_param_dict, final_list):
        current_keyword = keyword_list[current_index]
        for current_simplified_param_dict in current_dict_list:
            current_key = current_simplified_param_dict[current_keyword]
            new_complete_param_dict = {
                **current_complete_param_dict,
                current_keyword: current_key
            }
            if empty_str in current_simplified_param_dict:
                next_layer_dict_list = current_simplified_param_dict[empty_str]
                data_param_list_generator_iter_func(
                    next_layer_dict_list, current_index + 1, new_complete_param_dict, final_list)
            else:
                if extra_key_default_value_dict is not None:
                    for extra_key, default_value in extra_key_default_value_dict.items():
                        new_complete_param_dict[extra_key] = default_parameter_extract(
                            current_simplified_param_dict, extra_key, default_value)
                new_complete_param_dict[miscellaneous_str] = {}
                for other_key, value in current_simplified_param_dict.items():
                    if (other_key not in new_complete_param_dict) or (
                            extra_key_default_value_dict is not None and other_key in extra_key_default_value_dict):
                        new_complete_param_dict[miscellaneous_str][other_key] = value
                final_list.append(new_complete_param_dict)

    def data_param_list_generator(param_raw_list):
        data_param_list = []
        data_param_list_generator_iter_func(
            param_raw_list, 0, {}, data_param_list)
        return data_param_list

    return data_param_list_generator


def collect_results_func_template(
        project_name_generator, total_param_list, project_parameter_key_list, obj_threshold=False,
        different_final_data_obj=None):
    obj_threshold_key = 'obj_threshold'

    def collect_results(final_data_obj):
        final_mapping_dict = {}
        for param_dict in total_param_list:
            project_parameter_tuple = tuple(param_dict[key] for key in project_parameter_key_list)
            project_name = project_name_generator(*project_parameter_tuple)
            if project_name not in final_mapping_dict:
                final_mapping_dict[project_name] = project_parameter_tuple
            if different_final_data_obj is not None:
                different_final_data_obj.load_current_result_label(project_name)
                final_data_obj.share_data(different_final_data_obj)
            else:
                final_data_obj.load_current_result_label(project_name)
            if obj_threshold:
                final_data_obj.final_information_dict[project_name][obj_threshold_key] = default_parameter_extract(
                    param_dict, obj_threshold_key, None)
        return final_mapping_dict

    return collect_results


def tissue_name_breakdown(raw_name):
    tissue_sep = CoreConstants.specific_tissue_sep
    if tissue_sep not in raw_name:
        return None, raw_name
    plus = '+'
    prefix = None
    search_start = 0
    main_body_list = []
    while True:
        sep_location = raw_name.find(tissue_sep, search_start)
        if sep_location == -1:
            break
        else:
            this_prefix = raw_name[search_start:sep_location]
            if prefix is None:
                prefix = this_prefix
            else:
                assert this_prefix == prefix
            body_start = sep_location + len(tissue_sep)
            next_plus_location = raw_name.find(plus, sep_location)
            if next_plus_location == -1:
                main_body_list.append(raw_name[body_start:])
                break
            else:
                next_search_start = next_plus_location + 1
                main_body_list.append(raw_name[body_start:next_search_start])
                search_start = next_search_start
    main_body = ''.join(main_body_list)
    return prefix, main_body


def excel_column_letter_to_0_index(raw_column_str):
    final_index = -1
    for loc_index, letter in enumerate(raw_column_str[::-1]):
        letter_index = ord(letter) - ord('A') + 1
        if not 0 <= letter_index <= 26:
            raise ValueError('Should pass letter between A to Z: {}'.format(raw_column_str))
        final_index += letter_index * (26 ** loc_index)
    return final_index


def reversible_flux_pair_dict_generator(
        flux_name_index_dict, flux_name_tissue_name_dict, specific_tissue_reversible_flux_pair_dict):
    def is_reverse_flux(_flux_name):
        return _flux_name.endswith(CoreConstants.reverse_reaction_name_suffix)

    def remove_reverse_suffix(_reverse_flux_name):
        return _reverse_flux_name[:-len(CoreConstants.reverse_reaction_name_suffix)]

    reversible_flux_name_dict = {}
    all_flux_name_list = []
    for tissue_name, each_tissue_specific_reversible_flux_pair_list in specific_tissue_reversible_flux_pair_dict.items():
        for bare_forward_flux_name, bare_reverse_flux_name in each_tissue_specific_reversible_flux_pair_list:
            forward_flux_name = tissue_specific_name_constructor(bare_forward_flux_name, (tissue_name,))
            reverse_flux_name = tissue_specific_name_constructor(bare_reverse_flux_name, (tissue_name,))
            reversible_flux_pair = (forward_flux_name, reverse_flux_name)
            all_flux_name_list.append(reversible_flux_pair)
            reversible_flux_name_dict[forward_flux_name] = reversible_flux_pair
            reversible_flux_name_dict[reverse_flux_name] = reversible_flux_pair
    for complete_flux_name in flux_name_index_dict.keys():
        if complete_flux_name in reversible_flux_name_dict:
            continue
        if is_reverse_flux(complete_flux_name):
            forward_flux_name = remove_reverse_suffix(complete_flux_name)
            reverse_flux_name = complete_flux_name
        else:
            forward_flux_name = complete_flux_name
            reverse_flux_name = reverse_reaction_name(complete_flux_name)
        if forward_flux_name in flux_name_index_dict and reverse_flux_name in flux_name_index_dict:
            reversible_flux_pair = (forward_flux_name, reverse_flux_name)
            reversible_flux_name_dict[forward_flux_name] = reversible_flux_pair
            reversible_flux_name_dict[reverse_flux_name] = reversible_flux_pair
            all_flux_name_list.append(reversible_flux_pair)
        else:
            all_flux_name_list.append(complete_flux_name)
        tissue_name, raw_flux_name = flux_name_tissue_name_dict[complete_flux_name]
    return reversible_flux_name_dict, all_flux_name_list


def tissue_specific_reversible_flux_title_constructor_generator(flux_name_tissue_name_dict):
    def reversible_flux_title_constructor(flux_name_0, flux_name_1):
        tissue_name_0, raw_flux_name_0 = flux_name_tissue_name_dict[flux_name_0]
        tissue_name_1, raw_flux_name_1 = flux_name_tissue_name_dict[flux_name_1]
        if tissue_name_0 == tissue_name_1:
            new_raw_flux_name = f'{raw_flux_name_0} - {raw_flux_name_1}'
            new_flux_name = tissue_specific_name_constructor(new_raw_flux_name, tissue_name_0)
        else:
            new_flux_name = f'{flux_name_0} - {flux_name_1}'
        return new_flux_name
    return reversible_flux_title_constructor
