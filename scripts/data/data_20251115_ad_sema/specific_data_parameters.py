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

u_labeled_glucose_ratio_list = [0.99] * 6
u_glucose_labeled_input_metabolite_dict = {
    "GLC_inf": [
        {
            "ratio_list": u_labeled_glucose_ratio_list,
            "abundance": 1,
        },
    ],
}

class Keyword(object):
    labeling = 'labeling'
    tissue = 'tissue'
    strain = 'strain'
    condition = 'condition'
    index = 'index'

    glucose = 'glucose'
    pyruvate = 'pyruvate'

    wt = 'wt'
    ad = 'ad'
    ctrl = 'ctrl'
    sema = 'semaglutide'
    ctrl_vs_sema = 'ctrl_vs_semaglutide'

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
        'glycerol 3-phosphate': 'glycerol 3-p',
    }

    display_tissue_name_dict = {
        ModelKeyword.hippocampus: 'hippo',
        ModelKeyword.hypothalamus: 'hypotha',
        ModelKeyword.q_muscle: 'quads',
        ModelKeyword.soleus_muscle: 'soleus',
        ModelKeyword.brown_adipose: 'BAT',
        ModelKeyword.colon: 'ileum',
    }


class SpecificParameters(CompleteDataset):
    data_name = 'data_20251115_ad_sema'

    def __init__(self):
        super(SpecificParameters, self).__init__()
        self.deconvolve_other_atoms = True
        self.mixed_compartment_list = ('c', 'm')
        self.special_compartment_tissue_dict = {
            ModelKeyword.serum: ('e', 'c', 'm')
        }
        self.current_direct = f'{Direct.data_direct}/{self.data_name}'
        file_name_list = [
            'TSQ2511-288-ST-GLC-AA-FLUX-raw.xlsx',
            'TSQ2511-288-ST-GLC-FLUX-raw.xlsx',
        ]
        self.file_name_list = file_name_list
        self.file_path_list = [f'{self.current_direct}/{file_name}' for file_name in file_name_list]
        self.pickle_file_path = f'{self.current_direct}/complete_data.pickle'
        self.input_metabolite_data_dict = input_mid_data_processor(u_glucose_labeled_input_metabolite_dict)
        self.combined_data = False

        # self.sheet_name_dict = {
        #     'Sheet-1': None,
        # }
        common_sheet_name = 'Sheet-1'
        index_col_name = 'Compound'
        serum_mfa_excluded_metabolite_dict = {}
        common_mfa_excluded_metabolite_dict = {
            # 'fructose 1,6-bisphosphate': None,
            # 'glucose 6-phosphate': None,
            # 'fructose 6-phosphate': None,
            # 'a-ketoglutarate': None,
            # '2-phosphoglycerate/3-phosphoglycerate': None,
            # 'glucose',
            # 'pyruvate': None,
        }
        self.mfa_exclude_metabolites_dict = {
            ModelKeyword.serum: serum_mfa_excluded_metabolite_dict,
            ModelKeyword.cortex: common_mfa_excluded_metabolite_dict,
            ModelKeyword.hippocampus: common_mfa_excluded_metabolite_dict,
            ModelKeyword.hypothalamus: common_mfa_excluded_metabolite_dict,
            ModelKeyword.liver: common_mfa_excluded_metabolite_dict,
            ModelKeyword.kidney: common_mfa_excluded_metabolite_dict,
            ModelKeyword.q_muscle: common_mfa_excluded_metabolite_dict,
            ModelKeyword.soleus_muscle: common_mfa_excluded_metabolite_dict,
            ModelKeyword.brown_adipose: common_mfa_excluded_metabolite_dict,
            ModelKeyword.colon: common_mfa_excluded_metabolite_dict,
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
        self.tissue_label_name_dict = {
            '0': ModelKeyword.serum,
            '1': ModelKeyword.cortex,
            '2': ModelKeyword.hippocampus,
            '3': ModelKeyword.hypothalamus,
            '4': ModelKeyword.liver,
            '5': ModelKeyword.kidney,
            '6': ModelKeyword.q_muscle,
            '7': ModelKeyword.soleus_muscle,
            '8': ModelKeyword.brown_adipose,
            '9': ModelKeyword.colon,
        }
        ctrl_condition_list = ['A', 'B', 'C', 'D', 'E', 'F']
        sema_condition_list = ['G', 'H', 'I', 'J', 'K', 'L']
        group_condition_name_dict = {
            **{key: (Keyword.ad, Keyword.ctrl, index + 1) for index, key in enumerate(ctrl_condition_list)},
            **{key: (Keyword.ad, Keyword.sema, index + 1) for index, key in enumerate(sema_condition_list)},
        }
        self.group_condition_name_dict = {
            **group_condition_name_dict
        }
        self.common_prefix = 'M-1ul-'
        excluded_column_name_dict = {}
        for pf1 in range(10):
            for pf2 in range(1, 3):
                column_name = f'{pf1}QC{pf2}-1ul'
                excluded_column_name_dict[column_name] = None
        self.excluded_column_name_dict = {
            'Expected RT': None,
            'm/z (Expected)': None,
            **excluded_column_name_dict,
        }
        self.compound_replace_dict = {
            'Aspartic acid*C4N7NO4': 'Aspartic acid*C4H7NO4'
        }
        self._complete_data_parameter_dict_dict = {
            file_name: {
                'xlsx_file_path': file_path,
                'xlsx_sheet_name': common_sheet_name,
                'index_col_name': index_col_name,
                'mixed_compartment_list': self.mixed_compartment_list,
                'to_standard_name_dict': data_metabolite_to_standard_name_dict,
                'specific_dimension_metabolite_dict': self.specific_dimension_metabolite_dict,
                'compound_replace_dict': self.compound_replace_dict,
                'excluded_column_name_dict': self.excluded_column_name_dict,
                'decipher_column_name': self._decipher_column_name,
                'mfa_exclude_metabolites_dict': self.mfa_exclude_metabolites_dict,
                'special_compartment_tissue_dict': self.special_compartment_tissue_dict
            }
            for file_name, file_path in zip(file_name_list, self.file_path_list)
        }

    @staticmethod
    def project_name_generator(labeling, strain_name, condition_name, index):
        return f'{labeling}__{strain_name}__{condition_name}__{index}'

    def add_data_sheet(self, sheet_name, current_data_dict):
        labeling_name = Keyword.glucose
        if labeling_name not in self.complete_dataset:
            self.complete_dataset[labeling_name] = {}
        for strain_name, each_strain_data in current_data_dict.items():
            if strain_name not in self.complete_dataset[labeling_name]:
                self.complete_dataset[labeling_name][strain_name] = {}
            for condition_name, each_condition_data in each_strain_data.items():
                if condition_name not in self.complete_dataset[labeling_name][strain_name]:
                    self.complete_dataset[labeling_name][strain_name][condition_name] = {}
                for index_num, each_index_data in each_condition_data.items():
                    if index_num not in self.complete_dataset[labeling_name][strain_name][condition_name]:
                        self.complete_dataset[labeling_name][strain_name][condition_name][index_num] = each_index_data
                    else:
                        self.complete_dataset[labeling_name][strain_name][condition_name][index_num].update(
                            each_index_data)

    def return_multi_tissue_mfa_data_dict(self):
        current_input_metabolite_data_dict = input_mid_data_processor(u_glucose_labeled_input_metabolite_dict)
        final_mfa_data_dict = {}
        final_result_information_dict = {}
        for labeling, each_labeling_data_dict in self.complete_dataset.items():
            for strain, each_strain_data_dict in each_labeling_data_dict.items():
                for condition, each_condition_data_dict in each_strain_data_dict.items():
                    for index_num, current_mid_data_dict in each_condition_data_dict.items():
                        project_name = self.project_name_generator(labeling, strain, condition, index_num)
                        mfa_data_obj = MFAData(
                            project_name, current_mid_data_dict, current_input_metabolite_data_dict,
                            combined_data=self.combined_data, ratio_dict_to_objective_func=None)
                        final_mfa_data_dict[project_name] = mfa_data_obj
                        final_result_information_dict[project_name] = {
                            Keyword.labeling: labeling,
                            Keyword.strain: strain,
                            Keyword.condition: condition,
                            Keyword.index: index_num,
                        }
        return final_mfa_data_dict, final_result_information_dict

    def _decipher_column_name(self, raw_column_name):
        assert raw_column_name.startswith(self.common_prefix)
        group_tissue_str = raw_column_name[len(self.common_prefix):]
        group_label = group_tissue_str[0]
        tissue_label = group_tissue_str[1:]
        strain_name, condition_name, index_num = self.group_condition_name_dict[group_label]
        tissue_name = self.tissue_label_name_dict[tissue_label]
        return strain_name, condition_name, index_num, tissue_name

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
        raw_data_frame = pd.read_excel(
            xlsx_file_path, sheet_name=xlsx_sheet_name, skiprows=skip_row_num,
            index_col=index_col_name, nrows=row_num, usecols=col_range)
        group_metabolite_row_dict = {}
        metabolite_name_formula_dict = {}
        current_metabolite_name = ''
        current_metabolite_group_list = []
        for row_index, complete_metabolite_name in enumerate(raw_data_frame.index):
            if complete_metabolite_name is np.nan:
                continue
            if complete_metabolite_name in compound_replace_dict:
                complete_metabolite_name = compound_replace_dict[complete_metabolite_name]
            metabolite_name, formula_or_mid_index, *possible_index = complete_metabolite_name.split('*')
            if formula_or_mid_index.startswith('m'):
                assert current_metabolite_name == metabolite_name
                num_str = formula_or_mid_index[1:]
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
                metabolite_name_formula_dict[current_metabolite_name] = formula_or_mid_index
        group_metabolite_row_dict[current_metabolite_name] = current_metabolite_group_list
        final_raw_data_dict = {}
        for column_name in raw_data_frame.columns:
            if column_name in excluded_column_name_dict:
                continue
            strain_name, condition_name, index_num, tissue_name = decipher_column_name(column_name)
            if strain_name not in final_raw_data_dict:
                final_raw_data_dict[strain_name] = {}
            if condition_name not in final_raw_data_dict[strain_name]:
                final_raw_data_dict[strain_name][condition_name] = {}
            if index_num not in final_raw_data_dict[strain_name][condition_name]:
                final_raw_data_dict[strain_name][condition_name][index_num] = {}
            current_raw_data_dict = final_raw_data_dict[strain_name][condition_name][index_num]
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
            (Keyword.ad, Keyword.ctrl): ColorConfig.normal_blue,
            (Keyword.ad, Keyword.sema): ColorConfig.orange,
        }
        color_name_without_index = {
            f'{strain}_{condition}': color
            for (strain, condition), color in color_dict_without_index.items()
        }
        legend_color_dict = {
            Keyword.ctrl_vs_sema: color_name_without_index,
        }
        return color_dict_without_index, legend_color_dict

    def result_information_dict_analysis(self, complete_result_information_dict):
        color_dict_without_index, _ = self._color_dict
        color_order_dict = {}
        total_index = 0
        each_group_replicate_num = 6
        for (exercise_name, condition_name), color in color_dict_without_index.items():
            for current_index in range(1, 1 + each_group_replicate_num):
                color_order_dict[(exercise_name, condition_name, current_index)] = (total_index + current_index, color)
            total_index += each_group_replicate_num

        return_parameter_list = []
        for result_label, current_result_information_dict in complete_result_information_dict.items():
            major_key = Keyword.ctrl_vs_sema
            strain = current_result_information_dict[Keyword.strain]
            condition = current_result_information_dict[Keyword.condition]
            index = current_result_information_dict[Keyword.index]
            minor_key_list = []
            minor_key_str = f'{strain}_{condition}_{index}'
            order_index, current_color = color_order_dict[(strain, condition, index)]
            return_parameter_list.append(
                (result_label, major_key, minor_key_list, minor_key_str, current_color, order_index))
        return return_parameter_list

    def tissue_mid_name_list_dict_constructor(self):
        mid_name_list_in_serum = [
            ['glucose', 'pyruvate', 'lactate', 'alanine', ],
            ['glutamate', 'glutamine', 'aspartate', 'glycerol 3-phosphate', ]
        ]

        mid_name_list_in_one_tissue = [
            ['glucose', '2-phosphoglycerate/3-phosphoglycerate', 'pyruvate', 'lactate', 'glycerol 3-phosphate',
             'citrate', ],
            ['glutamate', 'a-ketoglutarate', 'succinate', 'malate', 'aspartate', 'acetyl-coa'],
        ]

        brain_mid_name_list_in_one_tissue = [
            ['glucose', '2-phosphoglycerate/3-phosphoglycerate', 'pyruvate', 'lactate', 'glycerol 3-phosphate', ],
            ['citrate', 'a-ketoglutarate', 'glutamate', 'malate', 'aspartate', 'acetyl-coa'],
        ]

        tissue_mid_name_list_dict = {
            ModelKeyword.serum: mid_name_list_in_serum,
            ModelKeyword.cortex: mid_name_list_in_one_tissue,
            ModelKeyword.hippocampus: mid_name_list_in_one_tissue,
            ModelKeyword.hypothalamus: mid_name_list_in_one_tissue,
            ModelKeyword.liver: mid_name_list_in_one_tissue,
            ModelKeyword.kidney: mid_name_list_in_one_tissue,
            ModelKeyword.q_muscle: mid_name_list_in_one_tissue,
            ModelKeyword.soleus_muscle: mid_name_list_in_one_tissue,
            ModelKeyword.brown_adipose: mid_name_list_in_one_tissue,
            ModelKeyword.colon: mid_name_list_in_one_tissue,
        }
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
