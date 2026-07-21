# -*- coding: utf-8 -*-
"""S-plot（基于 OPLS-DA）

Panel F: 横轴 p[1] (协方差), 纵轴 p(corr)[1] (相关系数)
取 S 形两端高贡献+高相关的代谢物作为候选驱动/标志物
"""

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig, LineStyle,
    default_parameter_extract, GeneralElements, ColorConfig, FontWeight,
    np,
)

ScatterDataFigure = GeneralElements.ScatterDataFigure

from ...common.common_figure_materials import CommonColorDict as CCD


class SplotFigure(CompositeFigure):
    """S-plot 分析

    figure_data_parameter_dict:
        p_loadings: np.ndarray  — p[1] 协方差
        p_corr: np.ndarray  — p(corr)[1] 相关系数
        vip_scores: np.ndarray  — VIP
        candidates_mask: np.ndarray  — bool, 是否为候选代谢物
        metabolite_names: list[str]  — 代谢物名称
        p_corr_threshold: float  — |p(corr)| 阈值线
        tissue: str  — 组织名
        r2y: float
        q2: float
    """

    height_to_width_ratio = 0.85
    title_height = 0.04

    def __init__(self, figure_data_parameter_dict, **kwargs):
        p_loadings = default_parameter_extract(
            figure_data_parameter_dict, 'p_loadings', np.array([]), pop=True)
        p_corr = default_parameter_extract(
            figure_data_parameter_dict, 'p_corr', np.array([]), pop=True)
        candidates_mask = default_parameter_extract(
            figure_data_parameter_dict, 'candidates_mask', None, pop=True)
        metabolite_names = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_names', None, pop=True)
        p_corr_threshold = default_parameter_extract(
            figure_data_parameter_dict, 'p_corr_threshold', 0.5, pop=True)
        tissue = default_parameter_extract(
            figure_data_parameter_dict, 'tissue', '', pop=True)
        r2y = default_parameter_extract(
            figure_data_parameter_dict, 'r2y', 0.0, pop=True)
        q2 = default_parameter_extract(
            figure_data_parameter_dict, 'q2', 0.0, pop=True)

        if len(p_loadings) == 0 or len(p_corr) == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False)
            return

        # 分三类着色：候选上调、候选下调、不显著
        if candidates_mask is None:
            candidates_mask = np.abs(p_corr) > p_corr_threshold

        is_up = candidates_mask & (p_loadings > 0)
        is_down = candidates_mask & (p_loadings < 0)
        is_insig = ~candidates_mask

        x_coords = p_loadings
        y_coords = p_corr

        colors = []
        for i in range(len(x_coords)):
            if is_up[i]:
                colors.append(CCD.teal)
            elif is_down[i]:
                colors.append(CCD.gray)
            else:
                colors.append(CCD.light_gray)

        data_dict = {
            ParameterName.x_value_array: x_coords,
            ParameterName.y_value_array: (y_coords, None),
            ParameterName.marker_size: 5,
            ParameterName.marker_color: colors,
            ParameterName.scatter_param_dict: {
                ParameterName.alpha: 0.7,
                ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order,
            },
        }

        x_pad = max(0.01, (x_coords.max() - x_coords.min()) * 0.1)
        y_pad = max(0.1, (y_coords.max() - y_coords.min()) * 0.1)
        x_lim = (x_coords.min() - x_pad, x_coords.max() + x_pad)
        y_lim = (y_coords.min() - y_pad, y_coords.max() + y_pad)

        # 阈值线
        common_dash_config = {
            **DataFigureConfig.common_line_param_dict_generator(),
            ParameterName.edge_color: ColorConfig.gray,
            ParameterName.z_order: DataFigureConfig.line_z_order,
            ParameterName.edge_style: LineStyle.dash,
        }
        scatter_line_pair = [[
            [x_lim[0], x_lim[1]], [p_corr_threshold, p_corr_threshold],
            dict(common_dash_config),
        ], [
            [x_lim[0], x_lim[1]], [-p_corr_threshold, -p_corr_threshold],
            dict(common_dash_config),
        ]]

        # 标注 top 候选代谢物
        supplementary_text, supplementary_text_loc = [], []
        if metabolite_names is not None and candidates_mask.any():
            cand_indices = np.where(candidates_mask)[0]
            # 按 |p_corr| 排序，取 top 10
            top_idx = cand_indices[np.argsort(np.abs(p_corr[cand_indices]))[-10:]]
            x_range = x_lim[1] - x_lim[0]
            y_range = y_lim[1] - y_lim[0]
            for i in top_idx:
                name = str(metabolite_names[i])
                if len(name) > 18:
                    name = name[:18] + '..'
                supplementary_text.append(name)
                loc = np.array([
                    (x_coords[i] - x_lim[0]) / x_range,
                    (y_coords[i] - y_lim[0]) / y_range + 0.02,
                ])
                supplementary_text_loc.append(loc)

        title_format = {
            **DataFigureConfig.common_subplot_text_format_dict_generator(),
            ParameterName.font_size: 6,
        }
        tissue_display = {
            'cortex': 'Cortex', 'hippocampus': 'Hippocampus', 'serum': 'Serum',
        }.get(tissue, tissue)
        subplot_title = f'{tissue_display} S-plot (R2Y={r2y:.2f}, Q2={q2:.2f})'

        figure_width = 0.85
        figure_height = figure_width * self.height_to_width_ratio

        scatter_param = {
            ParameterName.ax_bottom_left_list: [Vector(0, 0)],
            ParameterName.ax_size_list: [Vector(1, 1)],
            ParameterName.color_dict: None,
            ParameterName.data_nested_list: [data_dict],
            ParameterName.figure_config_dict: {
                ParameterName.subplot_name_text_format_dict: title_format,
                ParameterName.supplementary_text_format_dict: {
                    **DataFigureConfig.common_supplementary_text_config_dict,
                    ParameterName.font_size: 4,
                    ParameterName.width: 0.1,
                    ParameterName.height: 0.05,
                },
            },
            ParameterName.scatter_line: [scatter_line_pair],
            ParameterName.x_lim_list: [x_lim],
            ParameterName.x_label_list: ['p[1]'],
            ParameterName.x_ticks_list: [[]],
            ParameterName.y_lim_list: [y_lim],
            ParameterName.y_label_list: ['p(corr)[1]'],
            ParameterName.y_ticks_list: [[]],
            ParameterName.subplot_name_list: [subplot_title],
            ParameterName.text_axis_loc_pair: Vector(0.5, 1.1),
            ParameterName.supplementary_text_list: [supplementary_text],
            ParameterName.supplementary_text_loc_list: [supplementary_text_loc],
            ParameterName.legend: False,
        }

        scatter_obj = ScatterDataFigure(
            scatter_param,
            bottom_left=Vector(0, 0),
            size=Vector(figure_width, figure_height),
        )

        subfigure_element_dict = {'splot': {'splot': scatter_obj}}

        self.total_width = figure_width
        self.height_to_width_ratio = figure_height / figure_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(figure_width, figure_height),
            background=False, **kwargs)
