from scripts.src.common.config import Direct, Keywords as CommonKeywords
from scripts.src.core.common.classes import TransformDict
from scripts.src.core.data.data_class import MIDData, MFAData, InputMetaboliteData
from ..common_input_metabolites import input_mid_data_processor
from ..config import np, pd, CoreConstants, ModelKeyword, CompleteDataset, standardize_metabolite_name


data_metabolite_to_standard_name_dict = TransformDict({
    'aspartate(1-)': 'aspartate',
    'd-aspartate': 'aspartate',
    '(r,s)-lactate': 'lactate',
    'd-ribulose-5-phosphate': 'ribulose 5-phosphate',
    '3pg': '3-phosphoglycerate',
    '2pg': '2-phosphoglycerate',
    '3-hydroxybutanoic acid': '3-hydroxybutyrate',
    'glutathione disulfide_pos': 'gssg',
    'glutathione': 'gsh',
    '-atp': 'atp',
    '-ump': 'ump',
    'c5h11o8p(5)': 'ribose 5-phosphate',
    'c6h13o9p(13)': 'g6p or f6p',
    '(s)-malate': 'malate',
    'd-glucose': 'glucose',
    'g6p': 'glucose 6-phosphate',
    'f6p': 'fructose 6-phosphate',
    'd-erythrose 4-phosphate': 'erythrose 4-phosphate',
    'd-ribose-5-phosphate': 'ribose 5-phosphate',
    'c5h9no3(5)': 'l-hydroxyproline',
    'c6h14o12p2(7)': 'fructose 1,6-bisphosphate',
    'fructose-1,6-bisphosphate': 'fructose 1,6-bisphosphate',
    'c20h32n6o12s2': 'l-palmitoylcarnitine',
})

glucose_6_labeled_input_metabolite_dict = {
    "GLC_e": [
        {
            "ratio_list":
                [
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                ],
            "abundance": 1,
        },
    ],
}

class Keyword(object):
    labeling = 'labeling'
    tissue = 'tissue'
    exercise = 'exercise'
    condition = 'condition'
    index = 'index'

    sucrose = 'sucrose'
    sucrose_2 = 'sucrose_2'

    whole_body = 'whole_body'
    gut = 'gut'

    ctrl = 'ctrl'
    mr = 'mr'
    tau = 'tau'
    mr_tau = 'mr_tau'

    no_exec = 'no_exec'
    with_exec = 'with_exec'

    prefix = 'prefix'

    labeling_mapping_dict = TransformDict(**{})
    strain_mapping_dict = TransformDict(**{

    })
    two_list = [1, 2]
    three_list = [1, 2, 3]
    four_list = [1, 2, 3, 4]
    index_average_dict = {

    }

    display_mid_name_dict = {
        '2-phosphoglycerate/3-phosphoglycerate': '2pg/3pg',
        'fructose 6-phosphate/glucose 6-phosphate': 'g6p/f6p',
        'glycerol 3-phosphate': 'glycerol 3-p',
        'glucose 6-phosphate': 'glucose 6-p',
        'fructose 6-phosphate': 'fructose 6-p',
        'fructose 1,6-bisphosphate': 'fructose 16-bp',
    }

    display_tissue_name_dict = {
        whole_body: 'whole_body',
        gut: 'gut',
    }


class SpecificParameters(CompleteDataset):
    data_name = 'example_data_20250815'

    def __init__(self):
        super(SpecificParameters, self).__init__()
        self.mixed_compartment_list = ('c', 'm')
        self.special_compartment_tissue_dict = {}

        self.current_direct = f'{Direct.data_direct}/{self.data_name}'
        self.pickle_file_path = f'{self.current_direct}/complete_data.pickle'

        self.input_metabolite_data_dict = input_mid_data_processor(glucose_6_labeled_input_metabolite_dict)
        self.combined_data = False
        self.deconvolve_other_atoms = False

        first_sheet_name_list = ['WB - + full', 'gut- full', 'gut (ctrl, MR-Tau) - + full']
        file_name_sheet_name_dict = {
            '13C6-Sucrose WB Gut full.xlsx': first_sheet_name_list,
        }
        self.first_sheet_name_list = first_sheet_name_list
        self.file_name_list = list(file_name_sheet_name_dict.keys())
        self.file_path_list = [f'{self.current_direct}/{file_name}' for file_name in file_name_sheet_name_dict.keys()]
        index_col_name = 'Name'
        gut_mfa_excluded_metabolite_dict = {
            'fructose 1,6-bisphosphate': None,
            'glucose 6-phosphate': None,
            'fructose 6-phosphate': None,
            'a-ketoglutarate': None,
            '2-phosphoglycerate/3-phosphoglycerate': None,
        }
        self.mfa_exclude_metabolites_dict = {
            Keyword.gut: gut_mfa_excluded_metabolite_dict,
        }
        self.specific_dimension_metabolite_dict = {
            'acetyl-coa': 2,
            'alanine': 3,
            'aspartate': 4,
            'serine': 3,
            'propionyl-coa': 3,
            'succinyl-coa': 4,
            'pyruvate': 3,
            'a-ketoglutarate': 5,
            'oxaloacetate': 4,
            'propionate': 3,
            'mma': 4,
            '3-hp': 3,
            'glucose': 6,
            'fructose 6-phosphate/glucose 6-phosphate': 6,
            '2-phosphoglycerate/3-phosphoglycerate': 3,
            'lactate': 3,
            'glycerol 3-phosphate': 3,
            'citrate/isocitrate': 6,
            'glutamate': 5,
            'glutamine': 5,
            'succinate': 4,
            'malate': 4,
        }
        self.compound_replace_dict = {
            'C4N7NO4': 'C4H7NO4'
        }
        self.condition_standard_name_dict = {
            'Ctrl': Keyword.ctrl,
            'MR': Keyword.mr,
            'Tau': Keyword.tau,
            'MR-Tau': Keyword.mr_tau,
            'MR+Tau': Keyword.mr_tau,
        }
        self.sheet_name_tissue_dict = {
            'WB - + full': Keyword.whole_body,
            'gut- full': Keyword.gut,
            'gut (ctrl, MR-Tau) - + full': Keyword.gut,
        }
        self.excluded_column_name_dict = {
            'ID': None,
            'MZ': None,
        }
        complete_data_parameter_dict_dict = {}
        for (file_name, each_file_sheet_name_list), file_path in zip(
                file_name_sheet_name_dict.items(), self.file_path_list):
            for sheet_index, each_sheet_name in enumerate(each_file_sheet_name_list):
                if sheet_index == 0:
                    skip_row_num = 0
                else:
                    skip_row_num = 1
                complete_data_parameter_dict_dict[f'{file_name}_{each_sheet_name}'] = {
                    'xlsx_file_path': file_path,
                    'xlsx_sheet_name': each_sheet_name,
                    'skip_row_num': skip_row_num,
                    'index_col_name': index_col_name,
                    'mixed_compartment_list': self.mixed_compartment_list,
                    'to_standard_name_dict': data_metabolite_to_standard_name_dict,
                    'specific_dimension_metabolite_dict': self.specific_dimension_metabolite_dict,
                    'compound_replace_dict': self.compound_replace_dict,
                    'excluded_column_name_dict': self.excluded_column_name_dict,
                    'decipher_column_name': self._decipher_column_name,
                    'mfa_exclude_metabolites_dict': self.mfa_exclude_metabolites_dict,
                    'special_compartment_tissue_dict': self.special_compartment_tissue_dict,
                }
        self._complete_data_parameter_dict_dict = complete_data_parameter_dict_dict

    @staticmethod
    def project_name_generator(labeling, exercise, condition_name, index):
        return f'{labeling}__{condition_name}{'+' if exercise == Keyword.with_exec else '-'}__{index}'

    def add_data_sheet(self, sheet_name, current_data_dict):
        if sheet_name.endswith(self.first_sheet_name_list[-1]):
            labeling_name = Keyword.sucrose_2
        else:
            labeling_name = Keyword.sucrose
        if labeling_name not in self.complete_dataset:
            self.complete_dataset[labeling_name] = {}
        for exercise, each_exercise_dict in current_data_dict.items():
            if exercise not in self.complete_dataset[labeling_name]:
                self.complete_dataset[labeling_name][exercise] = {}
            for condition_name, each_condition_data in each_exercise_dict.items():
                if condition_name not in self.complete_dataset[labeling_name][exercise]:
                    self.complete_dataset[labeling_name][exercise][condition_name] = {}
                for index_num, each_index_data in each_condition_data.items():
                    if index_num not in self.complete_dataset[labeling_name][exercise][condition_name]:
                        self.complete_dataset[labeling_name][exercise][condition_name][index_num] = each_index_data
                    else:
                        self.complete_dataset[labeling_name][exercise][condition_name][index_num].update(each_index_data)

    def return_multi_tissue_mfa_data_dict(self):
        current_input_metabolite_data_dict = input_mid_data_processor(glucose_6_labeled_input_metabolite_dict)
        final_mfa_data_dict = {}
        final_result_information_dict = {}
        for labeling, each_labeling_data_dict in self.complete_dataset.items():
            for exercise, each_exercise_dict in each_labeling_data_dict.items():
                for condition, each_condition_data_dict in each_exercise_dict.items():
                    for index_num, current_mid_data_dict in each_condition_data_dict.items():
                        project_name = self.project_name_generator(labeling, exercise, condition, index_num)
                        mfa_data_obj = MFAData(
                            project_name, current_mid_data_dict, current_input_metabolite_data_dict,
                            combined_data=self.combined_data, ratio_dict_to_objective_func=None)
                        final_mfa_data_dict[project_name] = mfa_data_obj
                        final_result_information_dict[project_name] = {
                            Keyword.labeling: labeling,
                            Keyword.exercise: exercise,
                            Keyword.condition: condition,
                            Keyword.index: index_num,
                        }
        return final_mfa_data_dict, final_result_information_dict

    def _decipher_column_name(self, raw_column_name):
        if raw_column_name.endswith('+'):
            exercise = Keyword.with_exec
        elif raw_column_name.endswith('-'):
            exercise = Keyword.no_exec
        else:
            raise ValueError()
        condition_str, index_str = raw_column_name[:-1].split(' ')
        condition_name = self.condition_standard_name_dict[condition_str]
        assert index_str.isdigit()
        index_num = int(index_str)
        return exercise, condition_name, index_num

    def target_metabolite_data_loader(
            self, xlsx_file_path, xlsx_sheet_name, mixed_compartment_list=('c', 'm'),
            index_col_name='Name', formula_col_name=None, skip_row_num=None, row_num=None, col_range=None, to_standard_name_dict=None,
            excluded_metabolite_name_set=None, specific_dimension_metabolite_dict=None, compound_replace_dict=None,
            excluded_column_name_dict=None, decipher_column_name=None,  mfa_exclude_metabolites_dict=None,
            special_compartment_tissue_dict=None):
        if col_range is not None:
            col_range = list(range(*col_range))
        if excluded_metabolite_name_set is None:
            excluded_metabolite_name_set = set()
        if to_standard_name_dict is None:
            to_standard_name_dict = TransformDict()
        tissue_name = self.sheet_name_tissue_dict[xlsx_sheet_name]
        raw_data_frame = pd.read_excel(
            xlsx_file_path, sheet_name=xlsx_sheet_name, skiprows=skip_row_num,
            index_col=index_col_name, nrows=row_num, usecols=col_range)
        group_metabolite_row_dict = {}
        metabolite_name_formula_dict = {}
        current_metabolite_name = ''
        current_metabolite_group_list = []
        c13_label = '[13C]'
        for row_index, complete_metabolite_name in enumerate(raw_data_frame.index):
            if complete_metabolite_name is np.nan:
                continue
            if formula_col_name is not None:
                chemical_formula = raw_data_frame.loc[complete_metabolite_name][formula_col_name].strip()
                if chemical_formula in compound_replace_dict:
                    chemical_formula = compound_replace_dict[chemical_formula]
            else:
                chemical_formula = None
            if c13_label in complete_metabolite_name:
                new_chemical = False
                c13_loc = complete_metabolite_name.index(c13_label)
                metabolite_name = complete_metabolite_name[c13_loc + len(c13_label):]
                num_str = complete_metabolite_name[:c13_loc]
            elif complete_metabolite_name.endswith('-0'):
                new_chemical = True
                metabolite_name = complete_metabolite_name[:-len('-0')]
                num_str = ''
            else:
                if '-' in complete_metabolite_name:
                    right_minus_loc = complete_metabolite_name.rindex('-')
                    tail_str = complete_metabolite_name[right_minus_loc + 1:]
                    if tail_str.isdigit():
                        new_chemical = False
                        metabolite_name = complete_metabolite_name[:right_minus_loc]
                        num_str = tail_str
                    else:
                        new_chemical = True
                        metabolite_name = complete_metabolite_name
                        num_str = ''
                else:
                    new_chemical = True
                    metabolite_name = complete_metabolite_name
                    num_str = ''
            metabolite_name = to_standard_name_dict[metabolite_name.strip().lower()]
            if not new_chemical:
                assert current_metabolite_name == metabolite_name
                if not num_str.isdigit():
                    raise ValueError()
                else:
                    isotopomer_index = int(num_str)
                    while len(current_metabolite_group_list) <= isotopomer_index:
                        current_metabolite_group_list.append([])
                    current_metabolite_group_list[isotopomer_index].append(row_index)
            else:
                group_metabolite_row_dict[current_metabolite_name] = current_metabolite_group_list
                current_metabolite_name = metabolite_name
                current_metabolite_group_list = [[row_index]]
                metabolite_name_formula_dict[current_metabolite_name] = chemical_formula
        group_metabolite_row_dict[current_metabolite_name] = current_metabolite_group_list
        final_raw_data_dict = {}
        for column_name in raw_data_frame.columns:
            if column_name in excluded_column_name_dict:
                continue
            exercise, condition_name, index_num = decipher_column_name(column_name)
            if exercise not in final_raw_data_dict:
                final_raw_data_dict[exercise] = {}
            if condition_name not in final_raw_data_dict[exercise]:
                final_raw_data_dict[exercise][condition_name] = {}
            if index_num not in final_raw_data_dict[exercise][condition_name]:
                final_raw_data_dict[exercise][condition_name][index_num] = {}
            current_raw_data_dict = final_raw_data_dict[exercise][condition_name][index_num]
            for complete_metabolite_name, this_metabolite_all_row_list in group_metabolite_row_dict.items():
                if complete_metabolite_name == '':
                    continue
                metabolite_formula_str = metabolite_name_formula_dict[complete_metabolite_name]
                standard_metabolite_name, combined_metabolite_list = standardize_metabolite_name(complete_metabolite_name)
                if standard_metabolite_name in excluded_metabolite_name_set:
                    continue
                raw_data_list = []
                invalid_data_index_list = []
                for isotopomer_index, each_row_list in enumerate(this_metabolite_all_row_list):
                    total_raw_data_element = 0
                    if len(each_row_list) == 0:
                        invalid_data_index_list.append(isotopomer_index)
                    else:
                        for metabolite_row in each_row_list:
                            current_raw_data_element = raw_data_frame[column_name].iloc[metabolite_row]
                            if (
                                    current_raw_data_element is None or
                                    current_raw_data_element == 'N/F' or
                                    np.isnan(current_raw_data_element)):
                                current_raw_data_element = 0
                            total_raw_data_element += current_raw_data_element
                    raw_data_list.append(total_raw_data_element)
                if len(raw_data_list) == 1 or sum(raw_data_list) == 0:
                    print(
                        f'[Data loader] Zero sum: Data from metabolite {complete_metabolite_name} '
                        f'in condition {column_name} in sheet {xlsx_sheet_name} is left out')
                    continue
                if standard_metabolite_name in specific_dimension_metabolite_dict:
                    specific_dim_num = specific_dimension_metabolite_dict[standard_metabolite_name]
                    while len(raw_data_list) <= specific_dim_num:
                        invalid_data_index_list.append(len(raw_data_list))
                        raw_data_list.append(0)
                if len(invalid_data_index_list) == 0:
                    invalid_data_index_list = None
                else:
                    invalid_data_index_list = tuple(invalid_data_index_list)
                if tissue_name in special_compartment_tissue_dict:
                    current_compartment_list = special_compartment_tissue_dict[tissue_name]
                else:
                    current_compartment_list = mixed_compartment_list
                new_mid_data_obj = MIDData(
                    raw_data_list=raw_data_list, raw_name=standard_metabolite_name,
                    combined_raw_name_list=combined_metabolite_list, to_standard_name_dict=to_standard_name_dict,
                    compartment_list=current_compartment_list, tissue_list=(tissue_name,),
                    invalid_index_list=invalid_data_index_list, chemical_formula_str=metabolite_formula_str)
                if invalid_data_index_list is not None:
                    new_mid_data_obj.normalize(0, deconvolve_other_atoms=self.deconvolve_other_atoms)
                else:
                    new_mid_data_obj.normalize(CoreConstants.eps_for_mid, deconvolve_other_atoms=self.deconvolve_other_atoms)
                current_raw_data_dict[new_mid_data_obj.full_name] = new_mid_data_obj
                if tissue_name in mfa_exclude_metabolites_dict:
                    if new_mid_data_obj.name in mfa_exclude_metabolites_dict[tissue_name]:
                        new_mid_data_obj.exclude_this_mid()
        return final_raw_data_dict

    @property
    def _color_dict(self):
        from scripts.figures.common.config import ColorConfig
        color_dict_without_index = {
            (Keyword.no_exec, Keyword.ctrl): ColorConfig.normal_blue,
            (Keyword.no_exec, Keyword.mr): ColorConfig.orange,
            (Keyword.no_exec, Keyword.tau): ColorConfig.light_green,
            (Keyword.no_exec, Keyword.mr_tau): ColorConfig.light_bright_mixed_mid_red,
            (Keyword.with_exec, Keyword.ctrl): ColorConfig.dark_blue,
            (Keyword.with_exec, Keyword.mr): ColorConfig.dark_orange,
            (Keyword.with_exec, Keyword.tau): ColorConfig.green,
            (Keyword.with_exec, Keyword.mr_tau): ColorConfig.mid_red,
        }
        common_legend_color_dict = {
            f'{condition}{'+' if exercise == Keyword.with_exec else '-'}': color
            for (exercise, condition), color in color_dict_without_index.items()
        }
        legend_color_dict = {
            Keyword.sucrose: common_legend_color_dict,
            Keyword.sucrose_2: common_legend_color_dict,
        }
        return color_dict_without_index, legend_color_dict

    def result_information_dict_analysis(self, complete_result_information_dict):
        keyword = Keyword
        data_name = SpecificParameters.data_name

        color_dict_without_index, _ = self._color_dict
        color_order_dict = {}
        total_index = 0
        each_group_replicate_num = 3
        for (exercise_name, condition_name), color in color_dict_without_index.items():
            for current_index in range(1, 1 + each_group_replicate_num):
                color_order_dict[(exercise_name, condition_name, current_index)] = (total_index + current_index, color)
            total_index += each_group_replicate_num

        return_parameter_list = []
        for result_label, current_result_information_dict in complete_result_information_dict.items():
            major_key = current_result_information_dict[keyword.labeling]
            exercise = current_result_information_dict[keyword.exercise]
            exer_str = '+' if exercise == keyword.with_exec else '-'
            condition = current_result_information_dict[keyword.condition]
            index = current_result_information_dict[keyword.index]
            minor_key_list = []
            minor_key_str = f'{condition}_{index}{exer_str}'
            order_index, current_color = color_order_dict[(exercise, condition, index)]
            return_parameter_list.append(
                (result_label, major_key, minor_key_list, minor_key_str, current_color, order_index))
        return return_parameter_list

    def tissue_mid_name_list_dict_constructor(self, **kwargs):
        mid_name_list_in_one_tissue = [
            [
                'glucose', 'fructose 6-phosphate/glucose 6-phosphate', '2-phosphoglycerate/3-phosphoglycerate',
                'pyruvate', 'lactate', 'glycerol 3-phosphate', 'alanine', ],
            [
                'citrate/isocitrate', 'glutamate', 'glutamine', 'succinate',
                'malate', 'aspartate', ],
        ]
        complete_tissue_mid_name_list_dict = {
            Keyword.whole_body: mid_name_list_in_one_tissue,
            Keyword.gut: mid_name_list_in_one_tissue,
        }
        tissue_mid_name_list_dict = complete_tissue_mid_name_list_dict
        return (
            tissue_mid_name_list_dict, self.specific_dimension_metabolite_dict,
            Keyword.display_mid_name_dict, Keyword.display_tissue_name_dict)

    def experimental_figure_config_dict_generator(self):
        _, legend_color_dict = self._color_dict
        return {
            'figure_size': (8, 18),
            # 'size': [0.8, 1.1],
            # 'scale': 0.7,
            'legend': True,
            'legend_color_dict': legend_color_dict,
        }
