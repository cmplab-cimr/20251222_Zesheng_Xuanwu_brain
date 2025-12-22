from .config import GeneralElements, ParameterName
from ..figure_elements.metabolism_volcano_figure import MetabolismVolcanoFigure
from ..figure_elements.mid_comparison_figure import MIDComparisonGridBarWithLegendDataFigure
from ..figure_elements.flux_comparison_figure import FluxComparisonGridBarWithLegendDataFigure


class Elements(GeneralElements):
    MetabolismVolcanoFigure = MetabolismVolcanoFigure
    MIDComparisonGridBarWithLegendDataFigure = MIDComparisonGridBarWithLegendDataFigure
    FluxComparisonGridBarWithLegendDataFigure = FluxComparisonGridBarWithLegendDataFigure
