from .packages import np, mp, scipy_comb, tqdm
from .config import CoreConstants


def remove_numerical_error(array, eps):
    array[np.abs(array) < eps] = 0


def natural_dist(carbon_num, c13_ratio=CoreConstants.natural_c13_ratio):
    c12_ratio = 1 - c13_ratio
    total_num = carbon_num + 1
    output = []
    for index in range(total_num):
        output.append(
            scipy_comb(carbon_num, index) * c13_ratio ** index * c12_ratio ** (carbon_num - index))
    return np.array(output)


def natural_dist_array(atom_num, isotope_abundance_vector):
    final_array = isotope_abundance_vector
    for _ in range(1, atom_num):
        final_array = np.convolve(final_array, isotope_abundance_vector, mode='full')
    return final_array


def np_list_conv(array_list):
    result_array = array_list[0]
    for input_array in array_list[1:]:
        result_array = np.convolve(result_array, input_array, mode='full')
    return result_array


def np_log_eps(experimental_vector, predicted_vector, eps):
    return -np.sum((experimental_vector + eps) * np.log(predicted_vector + eps))


def full_emu_name_constructor(metabolite_name, carbon_num=None):
    if carbon_num is None:
        return metabolite_name
    else:
        return '{}{}{}'.format(
            metabolite_name, CoreConstants.emu_carbon_list_str_sep, '1' * carbon_num)


def mix_flux_name_constructor(unique_id_name, index, case_name=None):
    if case_name is None:
        case_name = ''
    else:
        case_name = f'{case_name}_'
    return f'{CoreConstants.mix_ratio_prefix}{CoreConstants.mix_flux_sep}{case_name}{unique_id_name}_{index}'


def compartmental_mid_name_constructor(metabolite_name, compartment_name_list=None):
    if compartment_name_list is None:
        return metabolite_name
    if isinstance(compartment_name_list, str):
        compartment_name_str = compartment_name_list
    else:
        compartment_name_str = '_'.join(compartment_name_list)
    return '{}{}{}'.format(
        metabolite_name, CoreConstants.compartmental_mid_name_sep, compartment_name_str)


def tissue_specific_name_constructor(metabolite_or_reaction_name, tissue_name_list=None):
    if tissue_name_list is None:
        return metabolite_or_reaction_name
    if isinstance(tissue_name_list, str):
        tissue_name_str = tissue_name_list
    elif len(tissue_name_list) == 1 and tissue_name_list[0] is None:
        return metabolite_or_reaction_name
    else:
        tissue_name_str = '_'.join(tissue_name_list)
    return '{}{}{}'.format(
        tissue_name_str, CoreConstants.specific_tissue_sep, metabolite_or_reaction_name)


def group_emu_name_constructor(metabolite_name, case_name=None):
    if case_name is None:
        return metabolite_name
    else:
        return f'{case_name}_{metabolite_name}'


def isdigit(number):
    return isinstance(number, int) or isinstance(number, float)


def default_parameter_extract(option_dict, key, default_value):
    if key in option_dict:
        return option_dict[key]
    else:
        return default_value


def reverse_reaction_name(raw_reaction_id):
    return raw_reaction_id + CoreConstants.reverse_reaction_name_suffix


def search_sorted(item, seq, start, end):
    seq_len = end - start
    if seq_len <= 3:
        for index in range(start, end):
            if item == seq[index]:
                return index
        return -1
    else:
        mid_index = int((start + end) / 2)
        if item == seq[mid_index]:
            while mid_index != 0 and item == seq[mid_index - 1]:
                mid_index -= 1
            return mid_index
        else:
            index = search_sorted(item, seq, start, mid_index)
            if index != -1:
                return index
            index = search_sorted(item, seq, mid_index + 1, end)
            if index != -1:
                return index
            return -1


def check_if_subsequence(subseq, whole_seq, reverse=False):
    """
    Both subseq and whole_seq should be sorted in ascending order
    """
    subsequence = True
    if len(subseq) > len(whole_seq):
        return check_if_subsequence(whole_seq, subseq, True)
    elif len(subseq) == len(whole_seq):
        for subseq_item, whole_seq_item in zip(subseq, whole_seq):
            if subseq_item != whole_seq_item:
                subsequence = False
                break
    else:
        start_index = 0
        end_index = len(whole_seq)
        for element in subseq:
            current_index = search_sorted(element, whole_seq, start_index, end_index)
            if current_index == -1:
                subsequence = False
                break
            start_index = current_index + 1
    return subsequence, reverse


def check_if_mix_flux(flux_name: str):
    if flux_name.startswith(f'{CoreConstants.mix_ratio_prefix}{CoreConstants.mix_flux_sep}'):
        return True
    else:
        return False


def check_if_biomass_flux(flux_name: str):
    # return flux_name == CoreConstants.biomass_flux_id
    return flux_name.startswith(CoreConstants.biomass_flux_id)


def biomass_reaction_dict_constructor(biomass_metabolite_name_list):
    biomass_reaction_dict_list = []
    for biomass_reaction_index, biomass_metabolite_name in enumerate(biomass_metabolite_name_list):
        biomass_reaction_dict_list.append(
            {
                'id': f'{CoreConstants.biomass_flux_id}_{biomass_reaction_index + 1}',
                'sub': [
                    (biomass_metabolite_name, '',),
                    ],
                'pro': [(CoreConstants.biomass_metabolite_id, '')],
        },)
    return biomass_reaction_dict_list


def exchange_reaction_constructor(
        metabolite_carbon_num_list, original_compartment, target_compartment,
        reaction_specific_suffix=None, reversible=True):
    underline_original_compartment = f'_{original_compartment}'
    underline_target_compartment = f'_{target_compartment}'
    exchange_reaction_dict_list = []
    new_metabolite_name_list = []
    ord_a = ord('a')
    for (metabolite_name, carbon_num) in metabolite_carbon_num_list:
        assert metabolite_name.endswith(underline_original_compartment)
        bare_metabolite_name = metabolite_name[:-len(underline_original_compartment)]
        reaction_id = f'EXCHANGE_{bare_metabolite_name}_{original_compartment}_{target_compartment}'
        if reaction_specific_suffix is not None:
            reaction_id = f'{reaction_id}_{reaction_specific_suffix}'
        carbon_mapping_str = ''.join([chr(ord_a + i) for i in range(carbon_num)])
        new_metabolite_name = f'{bare_metabolite_name}{underline_target_compartment}'
        new_reaction_dict = {
            'id': reaction_id,
            'sub': [(metabolite_name, carbon_mapping_str,),],
            'pro': [(new_metabolite_name, carbon_mapping_str,)],
        }
        if reversible:
            new_reaction_dict['reverse'] = True
        exchange_reaction_dict_list.append(new_reaction_dict)
        new_metabolite_name_list.append(new_metabolite_name)
    return exchange_reaction_dict_list, new_metabolite_name_list


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


def split_total_num_to_process(total_num, processes_num):
    common_size = total_num // processes_num
    rest_size = total_num % processes_num
    each_process_num_list = (
            [common_size + 1] * rest_size +
            [common_size] * (processes_num - rest_size))
    return each_process_num_list


def tqdm_daemon_process_func(receive_pipe, total_num, display_title=None):
    pbar = tqdm.tqdm(
        total=total_num, smoothing=0, maxinterval=5,
        desc=display_title)
    while True:
        update_num = receive_pipe.recv()
        if update_num < 0:
            break
        pbar.update(update_num)


class ProgressBarExtraProcess(object):
    def __init__(self, display_progress_bar=False, target_size=0, display_title='', manual_start=False):
        self.display_progress_bar = display_progress_bar
        self.target_size = target_size
        self.display_title = display_title
        self.receive_pipe = None
        self.send_pipe = None
        self.new_daemon_process = None
        self.manual_start = manual_start
        self.started = False
        if self.display_progress_bar:
            (self.receive_pipe, self.send_pipe) = mp.Pipe()

    def __enter__(self):
        if self.display_progress_bar:
            self.new_daemon_process = mp.Process(
                target=tqdm_daemon_process_func, args=(self.receive_pipe, self.target_size, self.display_title))
            if not self.manual_start:
                self.started = True
                self.new_daemon_process.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.display_progress_bar:
            self.send_pipe.send(-1)
            if self.started:
                self.new_daemon_process.join()
            self.new_daemon_process.close()

    def start(self):
        if self.display_progress_bar and self.manual_start:
            self.started = True
            self.new_daemon_process.start()
            self.send_pipe.send(0)
