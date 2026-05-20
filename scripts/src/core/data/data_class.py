from ..common.packages import np, optimize, scipy_comb
from ..common.classes import default_transform_dict, DefaultDict
from ..common.config import CoreConstants
from ..common.functions import (
    natural_dist, np_list_conv, compartmental_mid_name_constructor, tissue_specific_name_constructor)

from .chemical_formula_dict import (metabolite_name_formula_dict, specific_dimension_metabolite_dict)


class InputMetaboliteData(object):
    def __init__(self, metabolite_name, c13_ratio_abundance_list=None, carbon_num=None):
        """
        If c13_ratio_abundance_list is not set, then carbon_num is required and this input metabolite will be regarded
        as unlabeled.
        :param metabolite_name: Metabolite name.
        :param c13_ratio_abundance_list: Ratio list of 13C for each carbon in this metabolite.
        :param carbon_num: Number of carbon atom in this metabolite.
        """
        self.name = metabolite_name
        if c13_ratio_abundance_list is None:
            if carbon_num is None:
                raise ValueError(
                    'Carbon number is not set for unlabeled input metabolite! {}'.format(metabolite_name))
            c13_ratio_list = [np.array([CoreConstants.natural_c13_ratio] * carbon_num)]
            normalized_abundance_list = np.array([1])
        else:
            self.distribution_num = len(c13_ratio_abundance_list)
            c13_ratio_list = []
            abundance_list = []
            carbon_num = -1
            for data_dict in c13_ratio_abundance_list:
                current_ratio_list = data_dict['ratio_list']
                if carbon_num != -1:
                    if carbon_num != len(current_ratio_list):
                        raise ValueError(
                            'InputMetabolite_init: Not valid 13C ratio list: {}'.format(current_ratio_list))
                else:
                    carbon_num = len(current_ratio_list)
                c13_ratio_list.append(np.array(current_ratio_list))
                abundance_list.append(data_dict['abundance'])
            if abs(sum(abundance_list) - 1) > 1e-3:
                normalized_abundance_list = None
                raise ValueError(
                    'InputMetabolite_init: Not valid abundance distribution: {}'.format(abundance_list))
            else:
                abundance_list = np.array(abundance_list)
                normalized_abundance_list = abundance_list / np.sum(abundance_list)
        self.c13_ratio_list = c13_ratio_list
        self.abundance_list = normalized_abundance_list
        self.carbon_num = carbon_num
        self.full_carbon_mid = None
        self.generate_mid_vector()

    def __hash__(self):
        return self.name.__hash__()

    def generate_mid_vector(self, emu_selected_carbon_list=None):
        if emu_selected_carbon_list is None and self.full_carbon_mid is not None:
            final_mid_vector = self.full_carbon_mid
        else:
            tmp_sum_mid_vector = None
            for c13_ratio_vector, abundance in zip(self.c13_ratio_list, self.abundance_list):
                if emu_selected_carbon_list is None:
                    convolution_array = np.vstack([1 - c13_ratio_vector, c13_ratio_vector]).T
                    convolution_list = list(convolution_array)
                else:
                    convolution_list = []
                    for carbon_index, selection_indicator in enumerate(emu_selected_carbon_list):
                        if selection_indicator == 1:
                            c13_ratio = c13_ratio_vector[carbon_index]
                            convolution_list.append(np.array([1 - c13_ratio, c13_ratio]))
                new_mid_array = np_list_conv(convolution_list)
                if tmp_sum_mid_vector is None:
                    tmp_sum_mid_vector = new_mid_array * abundance
                else:
                    tmp_sum_mid_vector += new_mid_array * abundance
            final_mid_vector = tmp_sum_mid_vector / np.sum(tmp_sum_mid_vector)
            if self.full_carbon_mid is None:
                self.full_carbon_mid = final_mid_vector
        return final_mid_vector


def vector_normalize(raw_data_vector, eps_for_mid, metabolite_name=None, total_sum_threshold=None, warn_func=None):
    raw_data_sum = np.sum(raw_data_vector)
    # data_length = len(raw_data_vector)
    correct_sum = raw_data_sum * (1 + eps_for_mid)
    if correct_sum < 1000 * eps_for_mid or (total_sum_threshold is not None and correct_sum < total_sum_threshold):
        if warn_func is not None:
            warn_func('[Normalization] Data from {} is left out'.format(metabolite_name))
        return None
    return (raw_data_vector + raw_data_sum * eps_for_mid) / correct_sum
    # return (raw_data_vector + raw_data_sum / data_length * eps_for_mid) / correct_sum


def abundance_dict_mass_dict_generator(isotope_abundance_mass_dict):
    isotope_abundance_dict = {}
    isotope_mass_dict = {}

    for atom_name, atom_mass_ratio_pair_list in isotope_abundance_mass_dict.items():
        mass_list, abundance_list = zip(*atom_mass_ratio_pair_list)
        isotope_abundance_dict[atom_name] = np.array(abundance_list)
        isotope_mass_dict[atom_name] = np.array(mass_list)
    return isotope_abundance_dict, isotope_mass_dict


class Atom(object):
    carbon = 'c'
    oxygen = 'o'
    nitrogen = 'n'
    sulfur = 's'
    hydrogen = 'h'
    phosphorus = 'p'

    m_plus_h = 'm+h'
    m_minus_h = 'm-h'
    proton_mass = 1.007_276

    natural_c13_ratio = CoreConstants.natural_c13_ratio
    natural_o17_ratio = 0.00038
    natural_o18_ratio = 0.002
    natural_n15_ratio = 0.0038
    natural_h2_ratio = 0.000145
    natural_s33_ratio = 0.00763
    natural_s34_ratio = 0.04365

    isotope_offset_abundance_dict = {
        carbon: [(1, natural_c13_ratio)],
        oxygen: [(2, natural_o18_ratio)],
        nitrogen: [(1, natural_n15_ratio)],
        sulfur: [(1, natural_s33_ratio), (2, natural_s34_ratio)],
        hydrogen: [(1, natural_h2_ratio)],
    }

    isotope_precise_mass_abundance_dict = {
        carbon: [
            (12.0, 1 - natural_c13_ratio),
            (13.003_355, natural_c13_ratio),
        ],
        oxygen: [
            (15.994_914, 1 - natural_o17_ratio - natural_o18_ratio),
            (16.999_131, natural_o17_ratio),
            (17.999_159, natural_o18_ratio),
        ],
        nitrogen: [
            (14.003_074, 1 - natural_n15_ratio),
            (15.000_108, natural_n15_ratio),
        ],
        hydrogen: [
            (1.007_825, 1 - natural_h2_ratio),
            (2.014_101, natural_h2_ratio),
        ],
        sulfur: [
            (31.972_071, 1 - natural_s33_ratio - natural_s34_ratio),
            (32.971_458, natural_s33_ratio),
            (33.967_867, natural_s34_ratio),
        ],
        phosphorus: [
            (30.973_761, 1),
        ],
    }

    isotope_abundance_dict, isotope_mass_dict = abundance_dict_mass_dict_generator(
        isotope_precise_mass_abundance_dict)

    complete_atom_dict = {carbon: 0, hydrogen: 0, oxygen: 0, nitrogen: 0, sulfur: 0, phosphorus: 0}

    isotope_num_mid_dict = {
        atom_name: [single_mid_array]
        for atom_name, single_mid_array in isotope_abundance_dict.items()
    }
    computed_atom_dict = {
        carbon: None, oxygen: None, nitrogen: None, sulfur: None, hydrogen: None
    }
    computed_chemical_str_formula_dict = {}
    computed_chemical_formula_mid_dict = {}

    isotope_num_precise_mass_list_dict = {
        atom_name: []
        for atom_name in isotope_mass_dict.keys()
    }

    @classmethod
    def chemical_formula_dict_generator(cls, chemical_formula_str):
        if chemical_formula_str is None or chemical_formula_str == '':
            return None
        elif not isinstance(chemical_formula_str, str):
            raise TypeError('')
        else:
            chemical_formula_str = chemical_formula_str.lower()
            if chemical_formula_str in cls.computed_chemical_str_formula_dict:
                return cls.computed_chemical_str_formula_dict[chemical_formula_str]
            else:
                chemical_formula_dict = {}
                str_index = 0
                formula_str_len = len(chemical_formula_str)
                while str_index < formula_str_len:
                    current_str = chemical_formula_str[str_index]
                    if current_str not in cls.complete_atom_dict or current_str in chemical_formula_dict:
                        raise ValueError()

                    next_str_index = str_index + 1
                    while next_str_index < formula_str_len and chemical_formula_str[next_str_index].isdigit():
                        next_str_index += 1
                    this_chemical_num_str = chemical_formula_str[str_index + 1:next_str_index]
                    if this_chemical_num_str == '':
                        this_chemical_num = 1
                    else:
                        this_chemical_num = int(this_chemical_num_str)
                    chemical_formula_dict[current_str] = this_chemical_num
                    str_index = next_str_index
                ordered_chemical_formula_dict = {
                    atom_label: chemical_formula_dict[atom_label]
                    for atom_label in cls.complete_atom_dict.keys()
                    if atom_label in chemical_formula_dict
                }
                cls.computed_chemical_str_formula_dict[chemical_formula_str] = ordered_chemical_formula_dict
                return ordered_chemical_formula_dict


    @classmethod
    def rare_isotope_vector_generator_old(cls, chemical_formula_dict):
        final_mid_conv_list = []
        for atom_label, atom_num in chemical_formula_dict.items():
            if atom_label not in cls.isotope_offset_abundance_dict:
                continue
            for current_isotope_offset, current_isotope_abundance in cls.isotope_offset_abundance_dict[atom_label]:
                if current_isotope_offset > 0:
                    raw_mid = natural_dist(atom_num, current_isotope_abundance)
                    if current_isotope_offset > 1:
                        inserted_zero_vector = np.zeros(current_isotope_offset - 1, dtype=float)
                        mid_vector_list = []
                        for element_index, each_element in enumerate(raw_mid):
                            if element_index != 0:
                                mid_vector_list.append(inserted_zero_vector)
                            mid_vector_list.append([each_element])
                        current_mid_vector = np.concatenate(mid_vector_list)
                    else:
                        current_mid_vector = raw_mid
                else:
                    raise ValueError()
                final_mid_conv_list.append(current_mid_vector)
        final_mid_array = np_list_conv(final_mid_conv_list)
        return final_mid_array

    @classmethod
    def rare_isotope_vector_generator(cls, chemical_formula_dict):
        def compute_num_repeat_mid(target_num, complete_mid_list):
            while len(complete_mid_list) < target_num:
                complete_mid_list.append(None)
            current_mid = complete_mid_list[target_num - 1]
            if current_mid is None:
                input_num_1 = target_num // 2
                input_num_2 = target_num - input_num_1
                input_mid_1 = compute_num_repeat_mid(input_num_1, complete_mid_list)
                input_mid_2 = compute_num_repeat_mid(input_num_2, complete_mid_list)
                current_mid = np.convolve(input_mid_1, input_mid_2, mode='full')
                complete_mid_list[target_num - 1] = current_mid
            return current_mid

        chemical_formula_str = ''.join([
            f'{atom_label}{str(atom_num)}'
            for atom_label, atom_num in chemical_formula_dict.items()])
        if chemical_formula_str in cls.computed_chemical_formula_mid_dict:
            final_mid_array = cls.computed_chemical_formula_mid_dict[chemical_formula_str]
        else:
            final_mid_conv_list = []
            for atom_label, atom_num in chemical_formula_dict.items():
                if atom_label not in cls.isotope_num_mid_dict:
                    raise ValueError()
                if atom_label not in cls.computed_atom_dict:
                    continue
                current_mid_vector = compute_num_repeat_mid(atom_num, cls.isotope_num_mid_dict[atom_label])
                final_mid_conv_list.append(current_mid_vector)
            final_mid_array = np_list_conv(final_mid_conv_list)
            cls.computed_chemical_formula_mid_dict[chemical_formula_str] = final_mid_array
        return final_mid_array

    @staticmethod
    def calculate_multinomial_probability(total_atom_num, counter_array, atom_abundance_array):
        if len(atom_abundance_array) == 1:
            prob = 1.0
        elif len(atom_abundance_array) == 2:
            prob = scipy_comb(total_atom_num, counter_array[1]) * np.prod(atom_abundance_array ** counter_array)
        else:
            inverse_exp_array = (atom_abundance_array ** counter_array)[::-1]
            comb_upper_list = []
            rest_atom = total_atom_num
            for current_atom in counter_array:
                comb_upper_list.append(rest_atom)
                rest_atom -= current_atom
            comb_upper_array = np.array(comb_upper_list)
            comb_array = scipy_comb(comb_upper_array, counter_array)
            prob = np.prod(comb_array * inverse_exp_array)
        return prob

    @classmethod
    def precise_mass_calculater(
            cls, chemical_formula_dict, adduct, mainly_traced_atom_tuple=('c',), max_len=None,
            probability_threshold=-1):
        def calculate_element_total_mass_value_abundance_list(
                current_atom_num, atom_mass_vector, atom_abundance_vector,
                current_element_isotope_precise_list, max_len):
            isotope_num = len(atom_mass_vector)
            if len(current_element_isotope_precise_list) > current_atom_num:
                current_list = current_element_isotope_precise_list[current_atom_num]
            else:
                current_list = current_element_isotope_precise_list[-1]
                while len(current_element_isotope_precise_list) <= current_atom_num:
                    last_max_atom_total_list = current_list
                    max_atom_num = len(current_element_isotope_precise_list)
                    max_offset = min(max_atom_num * (isotope_num - 1) + 1, max_len)
                    current_list = [{} for _ in range(max_offset)]
                    for atom_index, (atom_mass, atom_abundance) in enumerate(
                            zip(atom_abundance_vector, atom_mass_vector)):
                        for last_offset_index, last_each_offset_dict in enumerate(last_max_atom_total_list):
                            new_index = atom_index + last_offset_index
                            if new_index >= max_len:
                                continue
                            current_dict = current_list[new_index]
                            for last_counter_array, (last_total_mass, last_total_abundance) in last_each_offset_dict.items():
                                new_counter_array = np.array(last_counter_array)
                                new_counter_array[atom_index] += 1
                                new_counter_tuple = tuple(new_counter_array)
                                if new_counter_tuple not in current_dict:
                                    new_total_mass = np.sum(new_counter_array * atom_mass_vector)
                                    new_total_abundance = cls.calculate_multinomial_probability(
                                        max_atom_num, new_counter_array, atom_abundance_vector)
                                    current_dict[new_counter_tuple] = [new_total_mass, new_total_abundance]
                    current_element_isotope_precise_list.append(current_list)
            return current_list

        def compute_one_element_mass(
                element_isotope_precise_list, atom_mass_vector, atom_abundance_vector, repeat_num, max_len):
            if len(element_isotope_precise_list) == 0:
                one_list = []
                for atom_index, (atom_mass, atom_abundance) in enumerate(
                        zip(atom_mass_vector, atom_abundance_vector)):
                    current_one_hot_array = np.zeros_like(atom_mass_vector, dtype=int)
                    current_one_hot_array[atom_index] = 1
                    one_list.append({tuple(current_one_hot_array): (atom_mass, atom_abundance)})
                element_isotope_precise_list.extend([[], one_list])
            else:
                assert len(element_isotope_precise_list) >= 2
            final_mass_abundance_list = calculate_element_total_mass_value_abundance_list(
                repeat_num, atom_mass_vector, atom_abundance_vector,
                element_isotope_precise_list, max_len)
            return final_mass_abundance_list

        def calculate_rest_element_mass(
                result_element_list, molecule_max_len, total_mass_num_offset, total_mass_value, total_probability,
                total_untraced_probability,
                all_atom_category_combination_tuple, final_mass_list, atom_name_mass_abundance_vector_dict,
                molecule_mass_offset, total_probability_threshold=-1):
            if total_mass_num_offset > molecule_max_len or total_untraced_probability < total_probability_threshold:
                return
            elif len(result_element_list) == 0:
                while len(final_mass_list) <= total_mass_num_offset:
                    final_mass_list.append([])
                real_mass = total_mass_value + molecule_mass_offset
                final_mass_list[total_mass_num_offset].append(
                    (real_mass, (total_probability, total_untraced_probability), all_atom_category_combination_tuple))
            else:
                current_element = result_element_list[0]
                current_element_mass_list = atom_name_mass_abundance_vector_dict[current_element]
                rest_result_atom_list = result_element_list[1:]
                for mass_num_offset, atom_mass_dict in enumerate(current_element_mass_list):
                    new_mass_offset = total_mass_num_offset + mass_num_offset
                    for atom_category_counter_tuple, (each_element_mass, each_element_abundance) in atom_mass_dict.items():
                        if current_element in traced_atom_num_dict:
                            new_untraced_probability = total_untraced_probability
                        else:
                            new_untraced_probability = each_element_abundance * total_untraced_probability
                        new_probability = each_element_abundance * total_probability
                        new_mass_value = total_mass_value + each_element_mass
                        new_all_atom_category_combination_tuple = (
                            *all_atom_category_combination_tuple, (current_element, atom_category_counter_tuple))
                        calculate_rest_element_mass(
                            rest_result_atom_list, molecule_max_len, new_mass_offset, new_mass_value, new_probability,
                            new_untraced_probability, new_all_atom_category_combination_tuple, final_mass_list,
                            atom_name_mass_abundance_vector_dict, molecule_mass_offset, total_probability_threshold)

        if adduct == cls.m_plus_h:
            molecule_mass_offset = cls.proton_mass
        elif adduct == cls.m_minus_h:
            molecule_mass_offset = -cls.proton_mass
        else:
            raise ValueError()
        traced_atom_num_dict = {
            atom: chemical_formula_dict[atom] for atom in mainly_traced_atom_tuple}
        untraced_atom_num_dict = {
            atom_name: atom_num for atom_name, atom_num in chemical_formula_dict.items()
            if atom_name not in mainly_traced_atom_tuple}
        if max_len is None:
            traced_atom_offset_list = []
            for atom_label, atom_num in traced_atom_num_dict.items():
                current_isotope_dim = len(cls.isotope_abundance_dict[atom_label]) - 1
                traced_atom_offset_list.append(current_isotope_dim * atom_num)
            max_len = sum(traced_atom_offset_list) + 1 + 3
        all_atom_name_mass_abundance_vector_dict = {}
        for atom_label, atom_num in chemical_formula_dict.items():
            mass_vector = cls.isotope_mass_dict[atom_label]
            abundance_vector = cls.isotope_abundance_dict[atom_label]
            element_isotope_precise_list = cls.isotope_num_precise_mass_list_dict[atom_label]
            all_atom_name_mass_abundance_vector_dict[atom_label] = compute_one_element_mass(
                element_isotope_precise_list, mass_vector, abundance_vector, atom_num, max_len)
            # all_atom_name_mass_abundance_vector_dict[atom_label] = compute_one_atom_mass(
            #     mass_vector, abundance_vector, atom_num, max_len)
        raw_all_possible_mass_list = []
        initial_atom_list = list(chemical_formula_dict.keys())
        calculate_rest_element_mass(
            initial_atom_list, max_len, 0, 0, 1, 1,
            tuple(), raw_all_possible_mass_list, all_atom_name_mass_abundance_vector_dict,
            molecule_mass_offset, total_probability_threshold=probability_threshold)
        all_possible_mass_atom_category_list = [
            sorted(each_mass_offset_mass_list, key=lambda x: x[0])
            for each_mass_offset_mass_list in raw_all_possible_mass_list]
        return all_possible_mass_atom_category_list


complete_metabolite_dim_dict = {
    **{
        metabolite_name: Atom.chemical_formula_dict_generator(chemical_formula_str)[Atom.carbon]
        for metabolite_name, chemical_formula_str in metabolite_name_formula_dict.items()
        if metabolite_name not in specific_dimension_metabolite_dict
    },
    **specific_dimension_metabolite_dict
}


class MIDData(object):
    def __init__(
            self, raw_data_list=None, raw_data_vector: np.ndarray = None, raw_name: str = None,
            data_vector: np.ndarray = None, combined_raw_name_list=(), compartment_list=None, tissue_list=None,
            to_standard_name_dict=None, ms_total_sum_threshold=None, excluded_from_mfa=False, invalid_index_list=None,
            chemical_formula_str=None):
        if to_standard_name_dict is None:
            to_standard_name_dict = default_transform_dict
        combined = False
        if len(combined_raw_name_list) != 0:
            sorted_combined_name_list = sorted([to_standard_name_dict[name] for name in combined_raw_name_list])
            metabolite_name = CoreConstants.standard_metabolite_sep.join(sorted_combined_name_list)
            combined = True
            combined_standard_name_list = sorted_combined_name_list
        else:
            if raw_name is None:
                raise ValueError('Name must be given for non-combined MID model!')
            metabolite_name = to_standard_name_dict[raw_name]
            combined_standard_name_list = ()
        if raw_data_list is not None:
            true_raw_data_vector = np.array(raw_data_list)
            carbon_num = len(true_raw_data_vector) - 1
        elif raw_data_vector is not None:
            true_raw_data_vector = raw_data_vector
            carbon_num = len(true_raw_data_vector) - 1
        elif data_vector is not None:
            true_raw_data_vector = None
            carbon_num = len(data_vector) - 1
        else:
            raise ValueError('One of raw_data_list and raw_data_vector should be fed')
        if compartment_list is None:
            raise ValueError('Compartment is not specified')
        if tissue_list is None:
            tissue_list = (None,)

        self.combined = combined
        self.combined_standard_name_list = combined_standard_name_list
        self.raw_data_vector = true_raw_data_vector
        self.carbon_num = carbon_num
        self.data_vector = data_vector
        self.compartment = compartment_list
        self.tissue = tissue_list
        self.name = metabolite_name
        self.to_standard_name_dict = to_standard_name_dict
        self.ms_total_sum_threshold = ms_total_sum_threshold
        if len(compartment_list) > 1 or compartment_list[0] is None:
            compartment_list = None
        self.full_name = self.generate_full_name(metabolite_name, tissue_list, compartment_list)
        self.excluded_from_mfa = excluded_from_mfa
        self.invalid_index_list = invalid_index_list
        self.chemical_formula_dict = Atom.chemical_formula_dict_generator(chemical_formula_str)

    def generate_full_name(self, metabolite_name, tissue_list, compartment_list):
        return compartmental_mid_name_constructor(
            tissue_specific_name_constructor(metabolite_name, tissue_list),
            compartment_list)

    def exclude_this_mid(self):
        self.excluded_from_mfa = True

    def normalize(self, eps_for_mid=CoreConstants.eps_for_mid, deconvolve_other_atoms=False):
        if self.data_vector is None and self.raw_data_vector is not None:
            if deconvolve_other_atoms:
                real_eps_for_mid = 0
            else:
                real_eps_for_mid = eps_for_mid
            normalized_data_array = vector_normalize(
                self.raw_data_vector, real_eps_for_mid, self.name, total_sum_threshold=self.ms_total_sum_threshold,
                warn_func=print)
            if normalized_data_array is not None:
                if deconvolve_other_atoms and self.chemical_formula_dict is not None:
                    normalized_data_array = deconvolve_real_mid(
                        normalized_data_array, all_chemical_formula_dict=self.chemical_formula_dict,
                        minimal_mid=eps_for_mid)
                self.data_vector = normalized_data_array
                return True
        return False

    def copy(self):
        new_invalid_index_list = None if self.invalid_index_list is None else list(self.invalid_index_list)
        new_mid_data_obj = MIDData(
            raw_data_vector=np.copy(self.raw_data_vector),
            raw_name=self.name, combined_raw_name_list=list(self.combined_standard_name_list),
            compartment_list=self.compartment, tissue_list=self.tissue,
            to_standard_name_dict=self.to_standard_name_dict, ms_total_sum_threshold=self.ms_total_sum_threshold,
            invalid_index_list=new_invalid_index_list)
        if self.data_vector is not None:
            new_mid_data_obj.data_vector = np.copy(self.data_vector)
        new_mid_data_obj.excluded_from_mfa = self.excluded_from_mfa
        new_mid_data_obj.chemical_formula_dict = dict(self.chemical_formula_dict)
        return new_mid_data_obj

    def change_tissue(self, tissue_list):
        if tissue_list is None:
            tissue_list = (None,)
        else:
            assert isinstance(tissue_list, (list, tuple))
        self.tissue = tissue_list
        self.full_name = self.generate_full_name(self.name, self.tissue, self.compartment)

    def __repr__(self):
        if self.compartment is not None:
            name_str = '{}__{}'.format(self.name, self.compartment)
        else:
            name_str = self.name
        if self.data_vector is not None:
            data_str = self.data_vector.__repr__()
        else:
            data_str = self.raw_data_vector.__repr__()
        return '{}: {}'.format(name_str, data_str)


def deconvolve_real_mid(raw_mid, all_chemical_formula_dict, minimal_mid):
    counted_carbon_num = len(raw_mid) - 1
    deconvolve_chemical_formula_dict = {}
    for chemical_formula_key, atom_num in all_chemical_formula_dict.items():
        if chemical_formula_key == Atom.carbon:
            deconvolve_atom_num = atom_num - counted_carbon_num
            if deconvolve_atom_num > 0:
                deconvolve_chemical_formula_dict[chemical_formula_key] = deconvolve_atom_num
        else:
            deconvolve_chemical_formula_dict[chemical_formula_key] = atom_num
    rare_isotope_vector = Atom.rare_isotope_vector_generator(deconvolve_chemical_formula_dict)
    effective_isotope_vector = rare_isotope_vector[:counted_carbon_num + 1]
    original_real_mid = slsqp_deconvolution_solving(effective_isotope_vector, raw_mid, len(raw_mid))
    real_mid = np.maximum(original_real_mid, minimal_mid)
    real_mid /= np.sum(real_mid)
    return real_mid


def real_mid_correct_c13(real_mid_vector_with_natural_c13, carbon_num=None, minimal_mid=CoreConstants.eps_for_mid):
    if carbon_num is None:
        carbon_num = len(real_mid_vector_with_natural_c13) - 1
    assert carbon_num < len(real_mid_vector_with_natural_c13)
    natural_c13_isotope_vector = Atom.rare_isotope_vector_generator({Atom.carbon: carbon_num})
    original_real_mid = slsqp_deconvolution_solving(
        natural_c13_isotope_vector, real_mid_vector_with_natural_c13, len(real_mid_vector_with_natural_c13))
    real_mid = np.maximum(original_real_mid, minimal_mid)
    real_mid /= np.sum(real_mid)
    return real_mid


def slsqp_deconvolution_solving(effective_isotope_vector, measured_mid, real_mid_dim):
    def calculate_residue(x):
        return obj_matrix @ x - measured_mid

    def obj_func(x):
        residue = calculate_residue(x)
        return np.sum(residue ** 2)

    def jac_func(x):
        residue = calculate_residue(x)
        return 2 * obj_matrix.T @ residue

    def eq_cons(x):
        return sum_one_matrix @ x - 1

    def eq_jac(x):
        return sum_one_matrix

    x0 = np.zeros(real_mid_dim)
    x0[0] = 1
    effective_dim = len(effective_isotope_vector)
    before_pad_size = max(real_mid_dim - effective_dim, 0)
    padded_effective_vector = np.pad(effective_isotope_vector[::-1], (before_pad_size, real_mid_dim - 1))
    total_effective_dim = len(padded_effective_vector)
    obj_matrix_list = [
        padded_effective_vector[total_effective_dim - row_index - real_mid_dim : total_effective_dim - row_index]
        for row_index in range(real_mid_dim)
    ]
    obj_matrix = np.array(obj_matrix_list)
    sum_one_matrix = np.ones((1, real_mid_dim))
    cons = ({'type': 'eq', 'fun': eq_cons, 'jac': eq_jac},)

    res = optimize.minimize(
        obj_func, x0=x0, method='SLSQP', jac=jac_func, constraints=cons, bounds=optimize.Bounds(0, 1))
    if not res.success:
        raise ValueError(res.message)
    return res.x


def average_multiple_mid_data(mid_data_list, output_std=False):
    first_mid_obj = mid_data_list[0]
    num_mid_obj = len(mid_data_list)
    assert first_mid_obj.data_vector is not None
    data_vector_list = [first_mid_obj.data_vector]
    complete_data_vector = np.copy(first_mid_obj.data_vector)
    common_name = first_mid_obj.name
    common_combined_standard_name_list = first_mid_obj.combined_standard_name_list
    common_compartment = first_mid_obj.compartment
    common_tissue = first_mid_obj.tissue
    for new_mid_obj in mid_data_list[1:]:
        assert new_mid_obj.name == common_name
        assert new_mid_obj.combined_standard_name_list == common_combined_standard_name_list
        assert new_mid_obj.compartment == common_compartment
        assert new_mid_obj.tissue == common_tissue
        assert new_mid_obj.data_vector is not None
        complete_data_vector += new_mid_obj.data_vector
        data_vector_list.append(new_mid_obj.data_vector)
    final_data_vector = complete_data_vector / num_mid_obj
    new_mid_data_obj = MIDData(
        raw_name=common_name,
        data_vector=final_data_vector, combined_raw_name_list=common_combined_standard_name_list,
        compartment_list=common_compartment,
        tissue_list=common_tissue, to_standard_name_dict=first_mid_obj.to_standard_name_dict,
        ms_total_sum_threshold=first_mid_obj.ms_total_sum_threshold,
        invalid_index_list=first_mid_obj.invalid_index_list)
    if not output_std:
        return new_mid_data_obj
    else:
        std_data_vector = np.std(data_vector_list, axis=0)
        return new_mid_data_obj, std_data_vector


class MFAData(object):
    def __init__(
            self, data_name, experimental_mid_data_obj_dict, input_metabolite_obj_data_dict,
            combined_data=False, ratio_dict_to_objective_func=None):
        self.data_name = data_name
        self.experimental_mid_data_obj_dict = experimental_mid_data_obj_dict
        self.input_metabolite_obj_data_dict = input_metabolite_obj_data_dict
        self.combined_data = combined_data
        if combined_data:
            list_of_case_name = list(experimental_mid_data_obj_dict.keys())
            example_mid_data_obj_key = list_of_case_name[0]
            mid_data_obj_dict_for_construction = experimental_mid_data_obj_dict[example_mid_data_obj_key]
            if ratio_dict_to_objective_func is None:
                ratio_dict_to_objective_func = DefaultDict(default_value=1)
        else:
            list_of_case_name = [data_name]
            mid_data_obj_dict_for_construction = experimental_mid_data_obj_dict
        self.ratio_dict_to_objective_func = ratio_dict_to_objective_func
        self.list_of_case_name = list_of_case_name
        self.mid_data_obj_dict_for_construction = mid_data_obj_dict_for_construction

    def copy(self):
        return MFAData(
            data_name=self.data_name,
            experimental_mid_data_obj_dict=self.experimental_mid_data_obj_dict.copy(),
            input_metabolite_obj_data_dict=self.input_metabolite_obj_data_dict.copy(),
            combined_data=self.combined_data,
            ratio_dict_to_objective_func=self.ratio_dict_to_objective_func,
        )

    def add_input_metabolite_obj_dict(
            self, new_labeled_buffer_input_metabolite_obj_dict):
        self.input_metabolite_obj_data_dict = {
            **self.input_metabolite_obj_data_dict,
            **new_labeled_buffer_input_metabolite_obj_dict,
        }


def average_multiple_mfa_data(mfa_data_list, output_std=False, common_data_name=None, filter_excluded_from_mfa=False):
    def average_experimental_mid_data_obj_dict(mid_data_obj_dict_list):
        mid_data_obj_list_dict = {}
        for mid_data_obj_dict in mid_data_obj_dict_list:
            for mid_name, mid_data_obj in mid_data_obj_dict.items():
                if filter_excluded_from_mfa and mid_data_obj.excluded_from_mfa:
                    continue
                if mid_name not in mid_data_obj_list_dict:
                    mid_data_obj_list_dict[mid_name] = []
                mid_data_obj_list_dict[mid_name].append(mid_data_obj)
        mean_mid_data_obj_dict = {}
        std_mid_data_vector_dict = {}
        for mid_name, mid_data_obj_list in mid_data_obj_list_dict.items():
            new_mid_data_obj, std_data_vector = average_multiple_mid_data(mid_data_obj_list, output_std=True)
            mean_mid_data_obj_dict[mid_name] = new_mid_data_obj
            std_mid_data_vector_dict[mid_name] = std_data_vector
        return mean_mid_data_obj_dict, std_mid_data_vector_dict

    first_mfa_obj = mfa_data_list[0]
    assert len(first_mfa_obj.experimental_mid_data_obj_dict) != 0
    mid_data_obj_dict_list = [first_mfa_obj.experimental_mid_data_obj_dict]
    if common_data_name is None:
        common_data_name = 'Average'
    common_input_metabolite_obj_data_dict = first_mfa_obj.input_metabolite_obj_data_dict
    common_combined_data = first_mfa_obj.combined_data
    common_ratio_dict_to_objective_func = first_mfa_obj.ratio_dict_to_objective_func
    common_list_of_case_name = first_mfa_obj.list_of_case_name
    for new_mfa_obj in mfa_data_list[1:]:
        assert new_mfa_obj.input_metabolite_obj_data_dict == common_input_metabolite_obj_data_dict
        assert new_mfa_obj.combined_data == common_combined_data
        if new_mfa_obj.combined_data:
            assert new_mfa_obj.list_of_case_name == common_list_of_case_name
        assert len(new_mfa_obj.experimental_mid_data_obj_dict) != 0
        assert new_mfa_obj.ratio_dict_to_objective_func == common_ratio_dict_to_objective_func
        mid_data_obj_dict_list.append(new_mfa_obj.experimental_mid_data_obj_dict)
    mean_mid_data_obj_dict, std_mid_data_vector_dict = average_experimental_mid_data_obj_dict(mid_data_obj_dict_list)
    new_mfa_data_obj = MFAData(
        common_data_name, mean_mid_data_obj_dict,
        common_input_metabolite_obj_data_dict, combined_data=common_combined_data,
        ratio_dict_to_objective_func=common_ratio_dict_to_objective_func)
    if not output_std:
        return new_mfa_data_obj
    else:
        return new_mfa_data_obj, std_mid_data_vector_dict


def input_metabolite_obj_dict_batch_generator(
        metabolite_name_carbon_num_list, isotope_label_ratio=0.99):
    labeled_input_metabolite_dict = {
        metabolite_name: InputMetaboliteData(
            metabolite_name, [
                {
                    "ratio_list": [isotope_label_ratio] * carbon_num,
                    "abundance": 1,
                },
            ])
        for metabolite_name, carbon_num in metabolite_name_carbon_num_list
    }
    return labeled_input_metabolite_dict
