from scripts.src.core.common.config import ModelKeyword, ParamName, CoreConstants
from scripts.src.core.common.functions import tissue_specific_name_constructor, natural_dist
from scripts.src.core.data.data_class import MIDData as BasicMIDData, MFAData, average_multiple_mfa_data
from scripts.src.core.common.classes import MFAConfig, OptionDict

from scripts.src.common_parallel_solver.config import pickle_load, pickle_save, check_and_mkdir_of_direct
from ..src.common.config import Direct
from ..src.common.packages import defaultdict, datetime, pd, openpyxl, np, linregress, pathlib
from ..src.common.functions import excel_column_letter_to_0_index as letter_to_index


class Constants(CoreConstants):
    ms_total_intensity_threshold = 5 * 10 ** 6
    eps_for_mid = 1e-10
    normalization_mid = 1e-5


def standardize_metabolite_name(raw_name: str):
    standard_name = raw_name
    standard_name = standard_name.strip()
    standard_name = standard_name.strip('"')
    standard_name = standard_name.lower()
    for common_metabolite_sep in Constants.common_metabolite_sep_list:
        standard_name = standard_name.replace(common_metabolite_sep, Constants.standard_metabolite_sep)
    if standard_name.endswith('-)'):
        standard_name = standard_name[:standard_name.rindex('(')]
    if standard_name.endswith('-pos'):
        standard_name = standard_name[:-len('-pos')]
    elif standard_name.endswith('_pos'):
        standard_name = standard_name[:-len('_pos')]
    combined_metabolite_list = []
    if Constants.standard_metabolite_sep in standard_name:
        combined_metabolite_list = standard_name.split(Constants.standard_metabolite_sep)
    return standard_name, combined_metabolite_list


class MIDData(BasicMIDData):
    def __init__(self, *args, compartment_list=None, **kwargs):
        if compartment_list is None:
            compartment_list = ('c', 'm')
        BasicMIDData.__init__(self, *args, compartment_list=compartment_list, **kwargs)

    def sum(self):
        sum_data = np.sum(self.raw_data_vector)
        return MetabolomicData(raw_data=sum_data, mid_data=self)


class MetabolomicData(MIDData):
    def __init__(
            self, raw_data=None, raw_name=None, combined_raw_name_list=(), compartment_list=None, tissue_list=None,
            mid_data=None, to_standard_name_dict=None, ms_total_sum_threshold=None):
        if raw_data is None:
            raise ValueError('Must contain raw data!')
        if mid_data is not None:
            raw_name = mid_data.name
            combined_raw_name_list = mid_data.combined_standard_name_list
            compartment_list = mid_data.compartment
            tissue_list = mid_data.tissue
            to_standard_name_dict = mid_data.to_standard_name_dict
            ms_total_sum_threshold = mid_data.ms_total_sum_threshold
        self.raw_data = raw_data
        super().__init__(
            raw_data_vector=np.array([raw_data]), raw_name=raw_name,
            combined_raw_name_list=combined_raw_name_list, compartment_list=compartment_list,
            tissue_list=tissue_list, to_standard_name_dict=to_standard_name_dict,
            ms_total_sum_threshold=ms_total_sum_threshold)

    def normalize(self, eps_for_mid=CoreConstants.eps_for_mid):
        raise ValueError('Metabolomic data cannot be normalized!')

    def copy(self):
        new_metabolomic_data_obj = super().copy()
        return MetabolomicData(raw_data=self.raw_data, mid_data=new_metabolomic_data_obj)

    def modify_data(self, new_raw_data):
        self.raw_data = new_raw_data

    def __repr__(self):
        if self.compartment is not None:
            name_str = '{}__{}'.format(self.name, self.compartment)
        else:
            name_str = self.name
        data_str = self.raw_data.__repr__()
        return '{}: {}'.format(name_str, data_str)


def check_file_existence(file_path):
    file_path_obj = pathlib.Path(file_path)
    return file_path_obj.exists()



class CompleteDataset(object):
    def __init__(self):
        self.complete_dataset = {}
        self.anti_correction = False
        self._complete_data_parameter_dict_dict = None
        self._test_data_parameter_dict_dict = None

    def set_anti_correction(self, anti_correction=False):
        self.anti_correction = anti_correction

    def add_data_sheet(self, sheet_name, current_data_dict):
        pass

    def return_multi_tissue_mfa_data_dict(self):
        pass

    def return_data_parameter_dict(self):
        return self._complete_data_parameter_dict_dict

    def target_metabolite_data_loader(self, **kwargs):
        pass


class NaturalDistDict(dict):
    def __init__(self, *args, **kwargs):
        super(NaturalDistDict, self).__init__(*args, **kwargs)

    def __getitem__(self, item):
        if item not in self:
            new_natural_dist = natural_dist(item)
            self.__setitem__(item, new_natural_dist)
        return super(NaturalDistDict, self).__getitem__(item)


natural_dist_dict = NaturalDistDict()


def natural_distribution_anti_correction(raw_data_dict):
    for metabolite, corrected_mid_data_obj in raw_data_dict.items():
        current_natural_dist = natural_dist_dict[corrected_mid_data_obj.carbon_num].copy()
        corrected_data_vector = corrected_mid_data_obj.data_vector
        zero_item = current_natural_dist[0]
        current_natural_dist[0] = 0
        raw_data_array = corrected_data_vector * zero_item + current_natural_dist
        corrected_mid_data_obj.data_vector = raw_data_array


def normalize_negative_data_array(raw_data_array):
    minimal_abs = -np.min(raw_data_array)
    new_data_array = raw_data_array + minimal_abs
    normalized_new_data_array = vector_normalize(new_data_array, CoreConstants.eps_for_mid)
    return normalized_new_data_array



metabolite_name_formula_dict = {
    'glucose': 'C6H12O6',
    'glucose 6-phosphate': 'C6H13O9P',
    'fructose 6-phosphate': 'C6H13O9P',
    'fructose 1,6-bisphosphate': 'C6H14O12P2',
    'glyceraldehyde 3-phosphate': 'C3H7O6P',
    '2-phosphoglycerate': 'C3H7O7P',
    '3-phosphoglycerate': 'C3H7O7P',
    'phosphoenolpyruvate': 'C3H5O6P',
    'pyruvate': 'C3H4O3',
    'lactate': 'C3H6O3',
    'glycerol 3-phosphate': 'C3H9O6P',
    'dihydroxyacetone phosphate': 'C3H7O6P',
    '6-phosphogluconolactone': 'C6H9O9P',
    '6-phospho-d-gluconate': 'C6H11O10P',

    # PPP
    'erythrose 4-phosphate': 'C4H9O7P',
    'sedoheptulose 7-phosphate': 'C7H15O10P',
    'ribose 5-phosphate': 'C5H11O8P',
    'xylulose 5-phosphate': 'C5H11O8P',
    'ribulose 5-phosphate': 'C5H11O8P',

    # TCA
    'citrate': 'C6H8O7',
    'isocitrate': 'C6H8O7',
    'aconitate': 'C6H6O6',
    'a-ketoglutarate': 'C5H6O5',
    'beta-hydroxybutyrate': 'C4H8O3',
    'succinate': 'C4H6O4',
    'fumarate': 'C4H4O4',
    'malate': 'C4H6O5',
    'oxaloacetate': 'C4H4O5',
    'acetoacetate': 'C4H6O3',
    'udp-d-glucose': 'C15H24N2O17P2',
    'glucose 1-phosphate': 'C6H13O9P',
    'glycerate': 'C3H6O4',

    # AAs
    'adp': 'C10H15N5O10P2',
    'gdp': 'C10H15N5O11P2',
    'atp': 'C10H16N5O13P3',
    'gtp': 'C10H16N5O14P3',
    'amp': 'C10H14N5O7P',
    'glutamate': 'C5H9NO4',
    'glutamine': 'C5H10N2O3',
    'serine': 'C3H7NO3',
    'cysteine': 'C3H7NO2S',
    'aspartate': 'C4H7NO4',
    'carnitine': 'C7H15NO3',
    'ornithine': 'C5H12N2O2',
    'citrulline': 'C6H13N3O3',
    'creatine': 'C4H9N3O2',
    'glycine': 'C2H5NO2',
    'alanine': 'C3H7NO2',
    'nadh': 'C21H29N7O14P2',
    'nadph': 'C21H30N7O17P3',
    'nad': 'C21H27N7O14P2',
    'nadp': 'C21H28N7O17P3',

    'acetyl-coa': 'C23H38N7O17P3S',
    'succinyl-coa': 'C25H40N7O19P3S',
    'coa': 'C21H36N7O16P3S',
}

specific_dimension_metabolite_dict = {
    'glucose': 6,
    'glucose 6-phosphate': 6,
    'fructose 6-phosphate': 6,
    'fructose 1,6-bisphosphate': 6,
    'glyceraldehyde 3-phosphate': 3,
    '2-phosphoglycerate': 3,
    '3-phosphoglycerate': 3,
    'phosphoenolpyruvate': 3,
    'pyruvate': 3,
    'lactate': 3,
    'glycerol 3-phosphate': 3,
    'dihydroxyacetone phosphate': 3,
    '6-phosphogluconolactone': 6,
    '6-phospho-d-gluconate': 6,
    'erythrose 4-phosphate': 4,
    'sedoheptulose 7-phosphate': 7,
    'ribose 5-phosphate': 5,
    'xylulose 5-phosphate': 5,
    'ribulose 5-phosphate': 5,
    'citrate': 6,
    'isocitrate': 6,
    'aconitate': 6,
    'a-ketoglutarate': 5,
    'beta-hydroxybutyrate': 4,
    'succinate': 4,
    'fumarate': 4,
    'malate': 4,
    'oxaloacetate': 4,
    'acetoacetate': 4,
    'udp-d-glucose': 15,
    'glucose 1-phosphate': 6,
    'glycerate': 3,
    'adp': 10,
    'gdp': 10,
    'atp': 10,
    'gtp': 10,
    'amp': 10,
    'glutamate': 5,
    'glutamine': 5,
    'serine': 3,
    'cysteine': 3,
    'aspartate': 4,
    'carnitine': 7,
    'ornithine': 5,
    'citrulline': 6,
    'creatine': 4,
    'glycine': 2,
    'alanine': 3,
    'nadh': 21,
    'nadph': 21,
    'nad': 21,
    'nadp': 21,
    'coa': 21,
    'acetyl-coa': 2,
    'succinyl-coa': 4,
}
