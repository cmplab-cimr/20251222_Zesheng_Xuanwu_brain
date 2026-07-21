from .config import GeneralElements, ParameterName
from ..figure_elements.metabolism_volcano_figure import MetabolismVolcanoFigure
from ..figure_elements.mid_comparison_figure import MIDComparisonGridBarWithLegendDataFigure
from ..figure_elements.flux_comparison_figure import FluxComparisonGridBarWithLegendDataFigure

from ..figure_content.elements import (
    MultiTissueOplsdaGridFigure,
    ClusteredHeatmapFigure,
    MetabolomicsVolcanoFigure,
    PathwayBubbleFigure,
    CrossTissueBarFigure,
    SplotFigure,
    PhenotypeCorrelationFigure,
    DiffCountBarFigure,
    PermutationGridFigure,
    MultiVolcanoGridFigure,
    QCFigure,
    UpSetFigure,
    EnrichmentLollipopFigure,
    OrganClassEnrichmentHeatmap,
    LipidMetaboliteDirectionFigure,
    LipidCrossTissuePhenotypeFigure,
    PathwayDiagramFigure,
    CorrelationScatterGridFigure,
    CognitionCorrelationFigure,
)


class Elements(GeneralElements):
    MetabolismVolcanoFigure = MetabolismVolcanoFigure
    MIDComparisonGridBarWithLegendDataFigure = MIDComparisonGridBarWithLegendDataFigure
    FluxComparisonGridBarWithLegendDataFigure = FluxComparisonGridBarWithLegendDataFigure

    # 代谢组学 Figure 4 Elements
    MultiTissueOplsdaGridFigure = MultiTissueOplsdaGridFigure
    ClusteredHeatmapFigure = ClusteredHeatmapFigure
    MetabolomicsVolcanoFigure = MetabolomicsVolcanoFigure
    PathwayBubbleFigure = PathwayBubbleFigure
    CrossTissueBarFigure = CrossTissueBarFigure
    SplotFigure = SplotFigure
    PhenotypeCorrelationFigure = PhenotypeCorrelationFigure
    DiffCountBarFigure = DiffCountBarFigure
    PermutationGridFigure = PermutationGridFigure
    MultiVolcanoGridFigure = MultiVolcanoGridFigure
    QCFigure = QCFigure
    UpSetFigure = UpSetFigure
    EnrichmentLollipopFigure = EnrichmentLollipopFigure
    OrganClassEnrichmentHeatmap = OrganClassEnrichmentHeatmap
    LipidMetaboliteDirectionFigure = LipidMetaboliteDirectionFigure
    LipidCrossTissuePhenotypeFigure = LipidCrossTissuePhenotypeFigure
    PathwayDiagramFigure = PathwayDiagramFigure
    CorrelationScatterGridFigure = CorrelationScatterGridFigure
    CognitionCorrelationFigure = CognitionCorrelationFigure
