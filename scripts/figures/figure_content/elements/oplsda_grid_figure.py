# -*- coding: utf-8 -*-
"""多组织 OPLS-DA 得分图网格

Panel A: 各组织分别做 OPLS-DA，得分图排成小网格(small multiples)
横轴 = 预测成分(组间差异)，纵轴 = 正交成分(组内变异)
"""

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, FontWeight,
    np,
)
from ...common.config import CommonElementConfig
from ...common.common_figure_materials import CommonColorDict as CCD

ScatterDataFigure = GeneralElements.ScatterDataFigure

# 组织英文映射
TISSUE_NAMES = {
    'cortex': 'Cortex', 'hippocampus': 'Hippocampus', 'hypothalamus': 'Hypothalamus',
    'liver': 'Liver', 'kidney': 'Kidney', 'fat': 'Adipose',
    'gastrocnemius': 'Gastrocnemius', 'quadriceps': 'Quadriceps', 'ileum': 'Ileum',
    'serum': 'Serum',
}


class OplsdaGridScatterDataFigure(ScatterDataFigure):
    """单个组织的 OPLS-DA 得分散点图"""

    def __init__(self, figure_data_parameter_dict, bottom_left: Vector, size: Vector, **kwargs):
        tissue = default_parameter_extract(
            figure_data_parameter_dict, 'tissue', 'unknown', pop=True)
        r2y = default_parameter_extract(
            figure_data_parameter_dict, 'r2y', 0.0, pop=True)
        q2 = default_parameter_extract(
            figure_data_parameter_dict, 'q2', 0.0, pop=True)

        # 组织名作为 subplot title
        tissue_display = TISSUE_NAMES.get(tissue, tissue)
        title_text = f'{tissue_display}\nR2Y={r2y:.2f} Q2={q2:.2f}'

        title_format_dict = {
            **DataFigureConfig.common_subplot_text_format_dict_generator(),
            ParameterName.font_size: 5,
        }

        figure_config_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_config_dict, {}, pop=True)
        figure_config_dict.setdefault(
            ParameterName.subplot_name_text_format_dict, title_format_dict)
        # 父类 BasicScatterDataFigure 直接索引该键, 处理后须放回
        figure_data_parameter_dict[ParameterName.figure_config_dict] = figure_config_dict

        figure_data_parameter_dict.setdefault(ParameterName.subplot_name_list, [title_text])
        figure_data_parameter_dict.setdefault(ParameterName.text_axis_loc_pair, Vector(0.5, 1.1))

        super().__init__(figure_data_parameter_dict, bottom_left, size, **kwargs)


class MultiTissueOplsdaGridFigure(CompositeFigure):
    """多组织 OPLS-DA 得分图网格

    figure_data_parameter_dict:
        tissue_results: dict[str, OPLSDAResult]  — {tissue: result}
        tissue_order: list[str]  — 组织排列顺序
        group_colors: dict[str, str]  — {group_name: color_hex}
    """

    subplot_width = 0.22
    subplot_height = 0.20
    subplot_gap_x = 0.03
    subplot_gap_y = 0.04
    title_height = 0.04

    def __init__(self, figure_data_parameter_dict, **kwargs):
        tissue_results = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_results', {}, pop=True)
        tissue_order = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_order', list(tissue_results.keys()), pop=True)
        group_colors = default_parameter_extract(
            figure_data_parameter_dict, 'group_colors', {}, pop=True)

        # 筛选有结果的组织
        valid_tissues = [t for t in tissue_order if t in tissue_results]
        n = len(valid_tissues)
        if n == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False)
            return

        # 5 列布局（10 个组织正好 5+5 两行）
        ncols = min(5, n)
        nrows = (n + ncols - 1) // ncols
        total_width = ncols * self.subplot_width + (ncols - 1) * self.subplot_gap_x
        total_height = nrows * self.subplot_height + (nrows - 1) * self.subplot_gap_y + self.title_height

        subfigure_element_dict = {}

        for idx, t in enumerate(valid_tissues):
            result = tissue_results[t]
            row, col = divmod(idx, ncols)
            x0 = col * (self.subplot_width + self.subplot_gap_x)
            y0 = (nrows - 1 - row) * (self.subplot_height + self.subplot_gap_y)

            scores = result.scores
            g1_name, g2_name = result.group_names
            g1_mask = scores['group'] == g1_name
            g2_mask = scores['group'] == g2_name

            # 构建散点数据
            g1_color = group_colors.get(g1_name, CCD.gray)
            g2_color = group_colors.get(g2_name, CCD.teal)

            data_dict = {
                ParameterName.x_value_array: np.concatenate([
                    scores.loc[g1_mask, 'Predictive'].values,
                    scores.loc[g2_mask, 'Predictive'].values,
                ]),
                ParameterName.y_value_array: (
                    np.concatenate([
                        scores.loc[g1_mask, 'Orthogonal'].values,
                        scores.loc[g2_mask, 'Orthogonal'].values,
                    ]),
                    None,
                ),
                ParameterName.marker_size: 6,
                ParameterName.marker_color: (
                    [g1_color] * g1_mask.sum() + [g2_color] * g2_mask.sum()
                ),
                ParameterName.scatter_param_dict: {
                    ParameterName.alpha: 0.85,
                    ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order,
                },
            }

            pred_vals = scores['Predictive'].values
            orth_vals = scores['Orthogonal'].values
            x_pad = max(0.5, (pred_vals.max() - pred_vals.min()) * 0.15)
            y_pad = max(0.5, (orth_vals.max() - orth_vals.min()) * 0.15)
            x_lim = (pred_vals.min() - x_pad, pred_vals.max() + x_pad)
            y_lim = (orth_vals.min() - y_pad, orth_vals.max() + y_pad)

            subplot_param_dict = {
                'tissue': t,
                'r2y': result.r2y,
                'q2': result.q2,
                ParameterName.ax_bottom_left_list: [Vector(0, 0)],
                ParameterName.ax_size_list: [Vector(1, 1)],
                ParameterName.color_dict: None,
                ParameterName.data_nested_list: [data_dict],
                ParameterName.scatter_line: [None],
                ParameterName.figure_config_dict: {},
                ParameterName.x_lim_list: [x_lim],
                ParameterName.x_label_list: ['t[1]'],
                ParameterName.x_ticks_list: [[]],
                ParameterName.y_lim_list: [y_lim],
                ParameterName.y_label_list: ['to[1]'],
                ParameterName.y_ticks_list: [[]],
                ParameterName.legend: False,
            }

            name = f'oplsda_{t}'
            subfigure_element_dict[name] = {
                name: OplsdaGridScatterDataFigure(
                    subplot_param_dict,
                    bottom_left=Vector(x0, y0),
                    size=Vector(self.subplot_width, self.subplot_height),
                ),
            }

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
