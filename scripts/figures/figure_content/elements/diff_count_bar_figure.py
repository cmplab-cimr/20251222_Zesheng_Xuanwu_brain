# -*- coding: utf-8 -*-
"""各器官差异代谢物计数堆叠条 (Figure 4 Panel B / 设计方案 B)

按代谢物 Class I 种类 (不分上调/下调) 堆叠, 每器官显示全局 top5 种类 + Others,
不同颜色 + PatchLegend 标注各段种类。用框架 BarDataFigure + stack_different_class。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...figure_plotting_package.basic_shape_elements.legend import PatchLegend
from ...common.common_figure_materials import CommonColorDict as CCD

BarDataFigure = GeneralElements.BarDataFigure

# top5 种类的区分色 (Others 固定浅灰)
DEFAULT_CLASS_COLORS = [CCD.teal, '#d62728', '#9467bd', '#2ca02c', '#ff7f0e', '#1f77b4', '#8c564b']
TISSUE_SHORT = {
    'cortex': 'Ctx', 'hippocampus': 'Hipp', 'hypothalamus': 'Hyp',
    'liver': 'Liv', 'kidney': 'Kid', 'fat': 'Fat',
    'gastrocnemius': 'Gas', 'quadriceps': 'Quad', 'ileum': 'Ile', 'serum': 'Ser',
}


class DiffCountBarFigure(CompositeFigure):
    """各器官差异代谢物 Class I 计数堆叠条 + 图例

    figure_data_parameter_dict:
        tissues: list[str]        — 组织 (横轴)
        counts: list[dict]        — 每组织 {class_label: count}, 各组织 keys 一致 (top5 + Others)
        class_colors: dict        — {class_label: color} (可选, 缺省自动)
        class_names: dict         — {class_label: legend 显示名} (可选, 缺省用 label 本身)
    """

    legend_center_y = 0.93
    legend_w = 0.62
    legend_h = 0.06
    bar_bottom = 0.10
    bar_height = 0.78
    bar_left = 0.08
    bar_width = 0.88

    def __init__(self, figure_data_parameter_dict, **kwargs):
        tissues = default_parameter_extract(
            figure_data_parameter_dict, 'tissues', [], pop=True)
        counts = default_parameter_extract(
            figure_data_parameter_dict, 'counts', [], pop=True)
        class_colors = default_parameter_extract(
            figure_data_parameter_dict, 'class_colors', None, pop=True)
        class_names = default_parameter_extract(
            figure_data_parameter_dict, 'class_names', None, pop=True)
        n = len(tissues)
        if n == 0 or len(counts) != n:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        segments = list(counts[0].keys())
        # 颜色: 'Others' 固定浅灰, 其余按 class_colors 或自动循环
        color_dict = {}
        auto_idx = 0
        for s in segments:
            if s == 'Others':
                color_dict[s] = CCD.light_gray
            elif class_colors and s in class_colors:
                color_dict[s] = class_colors[s]
            else:
                color_dict[s] = DEFAULT_CLASS_COLORS[auto_idx % len(DEFAULT_CLASS_COLORS)]
                auto_idx += 1
        name_dict = class_names or {s: s for s in segments}

        seg_arrays = {s: np.array([c.get(s, 0) for c in counts], dtype=float) for s in segments}
        tick_labels = [TISSUE_SHORT.get(t, t[:4]) for t in tissues]

        bar_param = {
            ParameterName.ax_bottom_left_list: [Vector(0, 0)],
            ParameterName.ax_size_list: [Vector(1, 1)],
            ParameterName.color_dict: color_dict,
            ParameterName.data_nested_list: [(seg_arrays, None)],
            ParameterName.array_len_list: [n],
            ParameterName.figure_config_dict: {
                ParameterName.column_width: 0.6,
                ParameterName.edge: 0.1,
            },
            ParameterName.stack_different_class: True,
            ParameterName.x_tick_labels_list: [tick_labels],
            ParameterName.y_label_list: ['# diff. metabolites'],
            ParameterName.legend: False,
        }
        bar_obj = BarDataFigure(
            bar_param,
            bottom_left=Vector(self.bar_left, self.bar_bottom),
            size=Vector(self.bar_width, self.bar_height))

        legend_obj = PatchLegend(
            Vector(0.5, self.legend_center_y), self.legend_w, self.legend_h,
            color_dict, name_dict, config_dict={})

        subfigure_element_dict = {
            'bar': {'bar': bar_obj},
            'legend': {'legend': legend_obj},
        }
        total_width = 1.0
        total_height = 1.0
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
