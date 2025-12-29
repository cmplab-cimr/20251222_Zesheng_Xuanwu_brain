from scripts.src.core.common.config import CoreConstants
from scripts.src.core.data.data_class import InputMetaboliteData

glucose_unlabeled_dict = {
    "ratio_list":
        [CoreConstants.natural_c13_ratio] * 6,
    "abundance": 1}


def input_mid_data_processor(input_raw_metabolite_dict):
    return {
        input_metabolite_name: InputMetaboliteData(input_metabolite_name, abundance_data_list)
        for input_metabolite_name, abundance_data_list in input_raw_metabolite_dict.items()}
