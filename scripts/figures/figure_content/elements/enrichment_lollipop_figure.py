# -*- coding: utf-8 -*-
"""化学类富集棒棒糖图 (Figure 4 Panel B / 设计方案 v2)

横轴 = 富集比, 1.0 竖虚线为无富集参考。
FDR 星号 + k/K 通过 BarDataFigure supplementary_text 机制传递，
1.0 参考线通过 cutoff 机制传递。

继承 CompleteLollipopFigure。
"""

import numpy as np

from ...common.config import ParameterName
from .lollipop_figure import CompleteLollipopFigure
from ...common.common_figure_materials import CommonColorDict as CCD, is_lipid_class, fdr_star


class EnrichmentLollipopFigure(CompleteLollipopFigure):
    """水平富集棒棒糖图，继承 CompleteLollipopFigure。

    figure_data_parameter_dict:
        recs: list[dict]   — 每类 {class, k, K, ratio, p, fdr}, 按 FDR 排序
        all_pct: dict      — {class: 该类全部% (K/N*100)}, 决定点大小
    """

    total_width = 1.35

    def __init__(self, figure_data_parameter_dict, **kwargs):
        recs = figure_data_parameter_dict.get('recs', [])
        recs.reverse()
        all_pct = figure_data_parameter_dict.get('all_pct', {})
        n = len(recs)
        y_lim = 2.1

        # FDR 星号 + k/K: x 接近 1, y 在 (y_edge, 1-y_edge) 等间距
        y_edge = 0.08
        x_right = 0.96
        star_x_offset = 0.05
        sup_texts, sup_locs = [], []
        for i, r in enumerate(recs):
            s = fdr_star(r['fdr'])
            rr = float(r['ratio']) if np.isfinite(r['ratio']) else 1.0
            y = y_edge + (1.0 - 2 * y_edge) * i / max(n - 1, 1)
            if s:
                sup_texts.append(s)
                sup_locs.append((rr / y_lim + star_x_offset, y))
            sup_texts.append(f"{r['k']}/{r['K']}")
            sup_locs.append((x_right, y))

        lipid_color_dict = {
            r['class']: (CCD.lipid_highlight if is_lipid_class(r['class']) else CCD.light_gray)
            for r in recs
        }
        new_figure_data_parameter_dict = {
            **figure_data_parameter_dict,
            ParameterName.figure_height: max(0.5, n * 0.075 + 0.18) if n else 0.5,
            ParameterName.data_nested_list: {
                r['class']: (r['ratio'], all_pct.get(r['class'], 1.0)) for r in recs
            },
            ParameterName.color_dict: figure_data_parameter_dict.get(
                ParameterName.color_dict, lipid_color_dict),
            ParameterName.horizontal: True,
            ParameterName.y_lim_list: [(0, y_lim)],
            ParameterName.x_label_list: '',
            ParameterName.y_label_list: 'Enrichment ratio',
            ParameterName.x_tick_labels_list: [r['class'] for r in recs],
            ParameterName.supplementary_text_list: [sup_texts],
            ParameterName.supplementary_text_loc_list: [sup_locs],
            ParameterName.cutoff: [1.0],
        }

        super().__init__(new_figure_data_parameter_dict, **kwargs)
