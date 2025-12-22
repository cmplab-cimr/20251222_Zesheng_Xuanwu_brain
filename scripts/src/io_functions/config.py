from .packages import pathlib, np
from scripts.config import Direct as GeneralDirect, FigureDataKeywords
from ..core.common.functions import mid_name_process, tissue_specific_name_constructor
from ..core.common.classes import TransformDict
from ..core.common.config import CoreConstants
from scripts.figures.figure_plotting_package.common.figure_data_format import BasicFigureData as RawBasicFigureData, \
    FigureData as RawFigureData
from scripts.figures.figure_plotting_package.common.classes import Vector
from scripts.figures.common.elements import Elements, ParameterName

class Keywords(object):
    average = 'average'
    experimental = 'experimental'


class FigureData(RawFigureData):
    def __init__(self, data_prefix, data_name):
        super().__init__(GeneralDirect.figure_raw_data_direct, data_prefix, data_name)


def check_and_mkdir_of_direct(direct_str, file_path=False):
    direct_obj = pathlib.Path(direct_str)
    if file_path:
        direct_obj = direct_obj.parent
    dir_stack = []
    while not direct_obj.exists():
        dir_stack.append(direct_obj)
        direct_obj = direct_obj.parent
    while len(dir_stack) != 0:
        missed_direct = dir_stack.pop()
        missed_direct.mkdir()

