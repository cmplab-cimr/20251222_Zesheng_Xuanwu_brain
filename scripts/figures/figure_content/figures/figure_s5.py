# -*- coding: utf-8 -*-
"""Figure S5: ¹³C flux 补充 (验证 + 完整数据 + 透明化)

S5A 全 MID 拟合散点 / S5B flux 单点+MW / S5C 置换占位 / S5D 全 flux 热图
S5E 谷氨酰胺+趋势 TCA / S5G 肝脏酮体 / S5H 神经元脂质合成个体值(占位) / S5I flux-代谢物矩阵(占位)
去 S5F (依赖正向同位素模拟)。复用 figure_5 的 FluxBar/MIDGoodness + 统计 helper。
"""

import numpy as np
import pandas as pd

from ...common.config import (
    ParameterName, Figure as GeneralFigure, Vector, CompositeFigure,
    GeneralElements, default_parameter_extract, calculate_center_bottom_offset)
from ...common.elements import Elements
from ...common.common_figure_materials import (
    CommonColorDict as CCD, rdbu_diverging_cmap, make_colorbar, fdr_star)
from .figure_5 import (
    _flux_stats, _load_mid, move_center, _place, GROUP_COLORS,
    FluxBarFigure, MIDGoodnessFigure)

Subfigure = Elements.Subfigure
HeatmapDataFigure = GeneralElements.HeatmapDataFigure

# ---- 各 panel 的 flux 列表 ----
FLUX_S5B = ['astrocyte--ACITL_c', 'astrocyte--PC_m', 'neuron--PYK_c', 'neuron--GPS_c',
            'liver--LDH_c', 'astrocyte--GLUD_c', 'neuron--GLND_m']
FLUX_S5E = ['neuron--GLND_m', 'astrocyte--GLUD_m', 'astrocyte--ICDH_m',
            'astrocyte--AKGD_m', 'astrocyte--PC_m']
FLUX_S5G = ['liver--HEX_c', 'liver--CS_m', 'liver--ACSS2_m',
            'astrocyte--BHB_usage', 'neuron--BHB_usage']
HEXMAP_FLUXES = [
    'liver--HEX_c', 'liver--PFK_c', 'liver--PYK_c', 'liver--CS_m', 'liver--ICDH_m', 'liver--AKGD_m',
    'astrocyte--HEX_c', 'astrocyte--PYK_c', 'astrocyte--CS_m', 'astrocyte--ICDH_m',
    'astrocyte--AKGD_m', 'astrocyte--PC_m', 'astrocyte--ACITL_c', 'astrocyte--LIPID_c',
    'astrocyte--GLNS_c', 'astrocyte--LDH_c',
    'neuron--PYK_c', 'neuron--CS_m', 'neuron--ICDH_m', 'neuron--LDH_c', 'neuron--GLND_m',
]


# ============================================================
# FluxHeatmapFigure: 全 flux log2FC 热图 + p 星号 (单列)
# ============================================================

class FluxHeatmapFigure(CompositeFigure):
    """反应 log2FC 单列热图 (Sema vs AD) + p 星号。"""
    def __init__(self, figure_data_parameter_dict, **kwargs):
        reactions = default_parameter_extract(figure_data_parameter_dict, 'reactions', [], pop=True)
        stats = _flux_stats(reactions)
        if stats.empty:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return
        log2fc = np.nan_to_num(stats['log2FC'].values.reshape(-1, 1))
        p_arr = stats['p_mw'].values.reshape(-1, 1)
        labels = [str(r).split('--')[-1] + '\n(' + str(r).split('--')[0][:3] + ')'
                  for r in stats['reaction']]

        def _star(x, y, v):
            return fdr_star(p_arr[y, x]) or None

        vmax = float(np.percentile(np.abs(log2fc), 95)) if len(log2fc) else 1.0
        vmax = max(vmax, 1.0)
        hm_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: rdbu_diverging_cmap(),
            ParameterName.data_nested_list: log2fc,
            ParameterName.data_lim_pair: (-vmax, vmax),
            ParameterName.figure_config_dict: {
                ParameterName.data_value_text_format_dict: {
                    ParameterName.basic_number_format_str: _star, 'font_size': 6}},
            ParameterName.hidden_axis_boundary: False,
            ParameterName.x_tick_labels_list: ['Sema/AD'],
            ParameterName.y_tick_labels_list: labels,
        }
        total_width = 0.50
        total_height = max(0.45, min(0.95, 0.042 * len(reactions) + 0.10))
        hm = HeatmapDataFigure(hm_param, bottom_left=Vector(0.28, 0.08),
                               size=Vector(0.10, total_height - 0.16))
        cbar = make_colorbar(rdbu_diverging_cmap(), -vmax, vmax,
                             Vector(0.40, 0.08), Vector(0.02, total_height - 0.16),
                             label='log2FC', orientation='vertical')
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__({'hm': {'hm': hm}, 'cbar': {'cbar': cbar}},
                         Vector(0, 0), Vector(total_width, total_height),
                         background=False, **kwargs)


def _place_fluxbar(fluxes, subfigure_bottom_left, subfigure_size):
    stats = _flux_stats(fluxes)
    if stats.empty:
        return None
    return _place(FluxBarFigure, {'stats_df': stats, 'group_colors': GROUP_COLORS},
                  subfigure_bottom_left, subfigure_size)


# ============================================================
# Panels
# ============================================================

class SubfigureA(Subfigure):
    """S5A: 全代谢物 MID 拟合散点。"""
    subfigure_label = 'a'
    subfigure_title = 'all_mid_goodness'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        mid = _load_mid()
        exp_pred = [(e, p) for s in mid for m in mid[s].values() for (e, p) in m.values()]
        print(f'Panel S5A: {len(exp_pred)} MID 点')
        element = _place(MIDGoodnessFigure, {'exp_pred': exp_pred},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureB(Subfigure):
    subfigure_label = 'b'
    subfigure_title = 'flux_individual_mw'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        element = _place_fluxbar(FLUX_S5B, subfigure_bottom_left, subfigure_size)
        if element is None:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)
            return
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureC(Subfigure):
    """S5C 置换/敏感性 — 占位 (需重跑 flux 拟合, 工作量过大)。"""
    subfigure_label = 'c'
    subfigure_title = 'permutation_placeholder'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        print('Panel S5C: 置换/敏感性本轮占位')
        super().__init__({}, subfigure_bottom_left, subfigure_size,
                         subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureD(Subfigure):
    """S5D: 全 flux log2FC 热图 + 星号。"""
    subfigure_label = 'd'
    subfigure_title = 'all_flux_heatmap'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        element = _place(FluxHeatmapFigure, {'reactions': HEXMAP_FLUXES},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureE(Subfigure):
    """S5E: 谷氨酰胺循环 + 趋势 TCA/回补。"""
    subfigure_label = 'e'
    subfigure_title = 'gln_trend_tca'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        element = _place_fluxbar(FLUX_S5E, subfigure_bottom_left, subfigure_size)
        if element is None:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)
            return
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureG(Subfigure):
    """S5G: 肝脏/酮体代谢。"""
    subfigure_label = 'g'
    subfigure_title = 'liver_ketone'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        element = _place_fluxbar(FLUX_S5G, subfigure_bottom_left, subfigure_size)
        if element is None:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)
            return
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureH(Subfigure):
    """S5H 神经元脂质合成个体值 — 占位 (待 ScatterDataFigure 实现)。"""
    subfigure_label = 'h'
    subfigure_title = 'neuron_lipid_transparency_placeholder'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        print('Panel S5H: 神经元脂质合成个体值 (占位, 待 ScatterDataFigure)')
        super().__init__({}, subfigure_bottom_left, subfigure_size,
                         subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


class SubfigureI(Subfigure):
    """S5I flux–代谢物相关矩阵 — 占位 (虚拟映射, 待 HeatmapDataFigure 矩阵)。"""
    subfigure_label = 'i'
    subfigure_title = 'flux_metabolite_matrix_placeholder'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        print('Panel S5I: flux-代谢物相关矩阵 (虚拟映射占位)')
        super().__init__({}, subfigure_bottom_left, subfigure_size,
                         subfigure_label=self.subfigure_label, subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Figure S5
# ============================================================

class Figure(GeneralFigure):
    figure_label = 'figure_s5'
    figure_title = 'Figure S5'

    def __init__(self):
        subfigure_class_list = [
            SubfigureA, SubfigureB, SubfigureC, SubfigureD,
            SubfigureE, SubfigureG, SubfigureH, SubfigureI,  # 去 S5F
        ]
        figure_layout_list = [
            (0.18, [(0.5, 'a'), (0.5, 'b')]),
            (0.10, [(1.0, 'c')]),
            (0.22, [(1.0, 'd')]),
            (0.18, [(0.5, 'e'), (0.5, 'g')]),
            (0.16, [(0.5, 'h'), (0.5, 'i')]),
        ]
        super().__init__(
            self.figure_label, subfigure_class_list, figure_layout_list, {},
            figure_title=self.figure_title)
