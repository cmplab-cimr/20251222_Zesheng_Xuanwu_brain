from .basic_shape_elements import Elements as BasicElements, common_legend_generator
from .metabolic_network.elements import Elements as MetabolicNetworkElements
from .data_figure.elements import Elements as DataFigureElements
from .diagrams import Elements as DiagramElements
from .common.common_functions import convert_theta_to_coordinate


class Elements(BasicElements, MetabolicNetworkElements, DataFigureElements, DiagramElements):
    BasicElements = BasicElements
    MetabolicNetworkElements = MetabolicNetworkElements
    DataFigureElements = DataFigureElements
    DiagramElements = DiagramElements

    convert_theta_to_coordinate = convert_theta_to_coordinate


def merge_dict_with_conflict(*dict_list):
    final_dict = {}
    for current_dict in dict_list:
        for key, value in current_dict.items():
            if key in final_dict:
                raise ValueError('Key {} exists in previous dict with value {}\nNew value: {}'.format(
                    key, final_dict[key], value))
            final_dict[key] = value
    return final_dict

