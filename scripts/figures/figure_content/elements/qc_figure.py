# -*- coding: utf-8 -*-
"""数据质控图 (Figure S4A)

左: QC 样本各代谢物 RSD% 直方图 (30% 虚线阈值)
右: PCA 散点 (QC vs 样本, QC 应聚集 = 仪器稳定)
QC 数据由 SubfigureA 从 ds 算好后传入。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...figure_plotting_package.basic_shape_elements.composite_figure_and_axes import DataFigureAxes
from ...common.common_figure_materials import CommonColorDict as CCD


class RSDDistrAxes(DataFigureAxes):
    """QC 代谢物 RSD% 直方图。"""

    def __init__(self, qc_rsd, bottom_left, size, **kwargs):
        self.qc_rsd = np.asarray(qc_rsd, dtype=float)
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        valid = self.qc_rsd[np.isfinite(self.qc_rsd) & (self.qc_rsd > 0)]
        if len(valid):
            ax.hist(valid, bins=30, color=CCD.teal, alpha=0.75,
                    edgecolor='white', linewidth=0.3)
        ax.axvline(30, color=CCD.gray, ls='--', lw=0.6)
        ax.set_xlabel('RSD of QC (%)', fontsize=6)
        ax.set_ylabel('# metabolites', fontsize=6)
        ax.tick_params(labelsize=5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        return ax, transform


class PCAScatterAxes(DataFigureAxes):
    """PCA 散点: QC (青, 应聚集) vs 样本 (灰)。"""

    def __init__(self, pc1, pc2, is_qc, bottom_left, size, **kwargs):
        self.pc1 = np.asarray(pc1, dtype=float)
        self.pc2 = np.asarray(pc2, dtype=float)
        self.is_qc = np.asarray(is_qc, dtype=bool)
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        ax.scatter(self.pc1[~self.is_qc], self.pc2[~self.is_qc],
                   s=6, c=CCD.gray, alpha=0.4, label='sample', zorder=3)
        ax.scatter(self.pc1[self.is_qc], self.pc2[self.is_qc],
                   s=14, c=CCD.teal, edgecolors='white', linewidths=0.3,
                   label='QC', zorder=4)
        ax.set_xlabel('PC1', fontsize=6)
        ax.set_ylabel('PC2', fontsize=6)
        ax.legend(fontsize=5, frameon=False)
        ax.tick_params(labelsize=5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        return ax, transform


class QCFigure(CompositeFigure):
    """数据质控: RSD 直方图 + PCA 散点

    figure_data_parameter_dict:
        qc_rsd: array  — QC 各代谢物 RSD%
        pca: dict      — {pc1, pc2, is_qc}
    """

    margin_left = 0.12
    margin_bottom = 0.14
    margin_right = 0.04
    margin_top = 0.04
    gap = 0.08
    panel_w = 0.42
    panel_h = 0.78

    def __init__(self, figure_data_parameter_dict, **kwargs):
        qc_rsd = default_parameter_extract(
            figure_data_parameter_dict, 'qc_rsd', np.array([]), pop=True)
        pca = default_parameter_extract(
            figure_data_parameter_dict, 'pca', {}, pop=True)

        total_width = self.margin_left + 2 * self.panel_w + self.gap + self.margin_right
        total_height = self.margin_bottom + self.panel_h + self.margin_top

        rsd_axes = RSDDistrAxes(
            qc_rsd,
            bottom_left=Vector(self.margin_left, self.margin_bottom),
            size=Vector(self.panel_w, self.panel_h))
        pca_axes = PCAScatterAxes(
            pca.get('pc1', np.array([])), pca.get('pc2', np.array([])),
            pca.get('is_qc', np.array([], dtype=bool)),
            bottom_left=Vector(self.margin_left + self.panel_w + self.gap, self.margin_bottom),
            size=Vector(self.panel_w, self.panel_h))

        subfigure_element_dict = {
            'rsd': {'rsd': rsd_axes},
            'pca': {'pca': pca_axes},
        }
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
