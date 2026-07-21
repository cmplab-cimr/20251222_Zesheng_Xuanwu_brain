# -*- coding: utf-8 -*-
"""代谢组学 Figure 4 的自定义 Element 类"""

from .oplsda_grid_figure import MultiTissueOplsdaGridFigure
from .clustered_heatmap_figure import ClusteredHeatmapFigure
from .volcano_figure import MetabolomicsVolcanoFigure
from .pathway_bubble_figure import PathwayBubbleFigure
from .cross_tissue_bar_figure import CrossTissueBarFigure
from .splot_figure import SplotFigure
from .phenotype_corr_figure import PhenotypeCorrelationFigure
from .diff_count_bar_figure import DiffCountBarFigure
from .permutation_figure import PermutationGridFigure
from .multi_volcano_grid_figure import MultiVolcanoGridFigure
from .qc_figure import QCFigure
from .upset_figure import UpSetFigure
from .enrichment_lollipop_figure import EnrichmentLollipopFigure
from .organ_class_enrichment_heatmap import OrganClassEnrichmentHeatmap
from .lipid_direction_heatmap import LipidMetaboliteDirectionFigure
from .lipid_cross_tissue_phenotype import LipidCrossTissuePhenotypeFigure
from .pathway_diagram_figure import PathwayDiagramFigure
from .correlation_scatter_grid_figure import CorrelationScatterGridFigure
from .cognition_corr_figure import CognitionCorrelationFigure

__all__ = [
    'MultiTissueOplsdaGridFigure',
    'ClusteredHeatmapFigure',
    'MetabolomicsVolcanoFigure',
    'PathwayBubbleFigure',
    'CrossTissueBarFigure',
    'SplotFigure',
    'PhenotypeCorrelationFigure',
    'DiffCountBarFigure',
    'PermutationGridFigure',
    'MultiVolcanoGridFigure',
    'QCFigure',
    'UpSetFigure',
    'EnrichmentLollipopFigure',
    'OrganClassEnrichmentHeatmap',
    'LipidMetaboliteDirectionFigure',
    'LipidCrossTissuePhenotypeFigure',
    'PathwayDiagramFigure',
    'CorrelationScatterGridFigure',
    'CognitionCorrelationFigure',
]
