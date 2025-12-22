from .data_figure import DataFigure
from .heatmap_data_figure import HeatmapConfig, ColorBarDataFigure, BasicHeatmapDataFigure
# from .violin_box_scatter_mix_data_figure import LossDistanceGridBoxScatterMixDataFigure
from .scatter_data_figure import BasicScatterDataFigure
from .bar_data_figure import BasicBarDataFigure, BasicMIDComparisonGridBarDataFigure, \
    BasicFluxErrorBarDataFigure, BasicSingleBarDataFigure
from .violin_box_data_figure import BasicViolinBoxDataFigure
from .histogram_data_figure import BasicHistogramDataFigure


class Elements(object):
    DataFigure = DataFigure
    HeatmapConfig = HeatmapConfig
    ColorBarDataFigure = ColorBarDataFigure
    HeatmapDataFigure = BasicHeatmapDataFigure
    ScatterDataFigure = BasicScatterDataFigure
    BarDataFigure = BasicBarDataFigure
    FluxErrorBarDataFigure = BasicFluxErrorBarDataFigure
    MIDComparisonGridBarDataFigure = BasicMIDComparisonGridBarDataFigure
    BasicSingleBarDataFigure = BasicSingleBarDataFigure
    ViolinBoxDataFigure = BasicViolinBoxDataFigure
    HistogramDataFigure = BasicHistogramDataFigure

