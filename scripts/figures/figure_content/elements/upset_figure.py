# -*- coding: utf-8 -*-
"""UpSet 交集图 (Figure S4G 数据表 + 跨器官交集)

用框架 BarDataFigure 画 exclusive 交集大小条 (风格统一, 取 top N 交集)。
x 轴标注交集涉及的组织 (缩写), 重点展示跨器官共享 vs 器官特异。
"""

from itertools import combinations

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...common.common_figure_materials import CommonColorDict as CCD

BarDataFigure = GeneralElements.BarDataFigure

TISSUE_SHORT = {
    'cortex': 'Ctx', 'hippocampus': 'Hipp', 'hypothalamus': 'Hyp',
    'liver': 'Liv', 'kidney': 'Kid', 'fat': 'Fat',
    'gastrocnemius': 'Gas', 'quadriceps': 'Quad', 'ileum': 'Ile', 'serum': 'Ser',
}


class UpSetFigure(CompositeFigure):
    """跨器官差异代谢物交集条 (复用 BarDataFigure)

    figure_data_parameter_dict:
        tissues: list[str]
        sets: dict  — {tissue: set(metabolite_id)} 各组织显著差异代谢物
        top_n: int  — 显示前 N 个最大交集 (默认 15)
    """

    def __init__(self, figure_data_parameter_dict, **kwargs):
        tissues = default_parameter_extract(
            figure_data_parameter_dict, 'tissues', [], pop=True)
        sets = default_parameter_extract(
            figure_data_parameter_dict, 'sets', {}, pop=True)
        top_n = default_parameter_extract(
            figure_data_parameter_dict, 'top_n', 15, pop=True)

        # 算 exclusive 交集 (该 combo 交集 - 其余组织并集)
        valid = [t for t in tissues if t in sets and sets[t]]
        intersections = []
        max_r = min(len(valid), 5)
        for r in range(1, max_r + 1):
            for combo in combinations(valid, r):
                inter = set.intersection(*[sets[t] for t in combo])
                others = [sets[t] for t in valid if t not in combo]
                excl = inter - set.union(*others) if others else inter
                if len(excl) > 0:
                    intersections.append((combo, len(excl)))
        intersections.sort(key=lambda kv: kv[1], reverse=True)
        top = intersections[:top_n]

        if not top:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        sizes = np.array([s for _, s in top], dtype=float)
        labels = [' + '.join(TISSUE_SHORT.get(t, t[:4]) for t in combo) for combo, _ in top]

        bar_param = {
            ParameterName.ax_bottom_left_list: [Vector(0, 0)],
            ParameterName.ax_size_list: [Vector(1, 1)],
            ParameterName.color_dict: {'intersection': CCD.teal},
            ParameterName.data_nested_list: [({'intersection': sizes}, None)],
            ParameterName.array_len_list: [len(top)],
            ParameterName.figure_config_dict: {
                ParameterName.column_width: 0.6,
                ParameterName.edge: 0.1,
            },
            ParameterName.x_tick_labels_list: [labels],
            ParameterName.y_label_list: ['exclusive intersection size'],
            ParameterName.legend: False,
        }
        bar_obj = BarDataFigure(
            bar_param,
            bottom_left=Vector(0.10, 0.20),
            size=Vector(0.86, 0.75))

        subfigure_element_dict = {'upset': {'upset': bar_obj}}
        plot_w = max(0.5, len(top) * 0.09)
        total_width = plot_w + 0.12
        total_height = 0.95
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
