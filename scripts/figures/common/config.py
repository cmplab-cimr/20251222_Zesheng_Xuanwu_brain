import enum

from ..figure_plotting_package.elements import Elements as GeneralElements
from ..figure_plotting_package.common.config import ParameterName as GeneralParameterName
from ..figure_plotting_package.common.third_party_packages import plt, mcolors
from ..figure_plotting_package.common.common_figure_materials import CommonElementConfig
from ..figure_plotting_package.common.figure_data_format import BasicFigureData as RawBasicFigureData
from ..figure_plotting_package.common.classes import Vector, VerticalAlignment, HorizontalAlignment, FontWeight, \
    LineStyle, FontStyle
from ..figure_plotting_package.common.common_functions import default_parameter_extract
from ..figure_plotting_package.common.color import ColorConfig

from ..figure_plotting_package.data_figure.config import (
    merge_axis_format_dict, ParameterName as DataParameterName, DataFigureConfig, Keywords,
    merge_complete_config_dict, generate_violin_config_dict)

from scripts.config import FigureDataKeywords, np
CompositeFigure = GeneralElements.CompositeFigure
TextBox = GeneralElements.TextBox

class Direct(object):
    figure_output_direct = 'figures/output_figure'


class Figure(GeneralElements.Figure):
    figure_output_direct = Direct.figure_output_direct
    top_margin_ratio = 0.01
    side_margin_ratio = 0.02


class ParameterName(DataParameterName):
    branch = 'branch'
    blunt = 'blunt'
    cycle = 'cycle'

    positive = 'positive'
    negative = 'negative'
    medium_data = 'medium_data'


class BasicFigureData(RawBasicFigureData):
    data_direct = Direct.figure_output_direct


class DataName(object):
    fangchao_fly_data = 'fangchao_data_fruit_fly'
