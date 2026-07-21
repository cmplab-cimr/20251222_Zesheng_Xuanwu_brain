# -*- coding: utf-8 -*-
"""通用棒棒糖图 (lollipop chart): BarDataFigure + each_group_scatter 圆点。

公共绘图类，不含任何脂质/代谢组学特异逻辑。color_dict 须由调用方提供。
"""

from ...common.config import (
    ParameterName, Vector, FontWeight, CompositeFigure, TextBox,
    default_parameter_extract, GeneralElements, CommonElementConfig, np, ColorConfig, DataFigureConfig,
)

BarDataFigure = GeneralElements.BarDataFigure


class LollipopDataFigure(BarDataFigure):
    """棒棒糖图：BarDataFigure 基座 + each_group_scatter 圆点。

    figure_data_parameter_dict:
        data_nested_list: dict       — {label: (stem_value, head_value)}, stem=bar 高度, head=圆点位置
        color_dict: dict             — {label: hex} 每项颜色
        x_label_list: str            — x 轴标签
        y_label_list: str            — y 轴标签
        x_tick_labels_list: list[str] — x 轴刻度标签
        horizontal: bool             — 水平模式
        marker_size: (min, max)      — 点大小范围，按 head_value 线性映射
    """

    def __init__(self, figure_data_parameter_dict, bottom_left: Vector, size: Vector, **kwargs):
        data_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.data_nested_list, [], pop=True)
        color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.color_dict, None, force=True)
        x_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.x_label_list, '', pop=True)
        y_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.y_label_list, '', pop=True)
        x_tick_labels = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.x_tick_labels_list, None, pop=True)
        horizontal = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.horizontal, False, pop=True)
        self.horizontal = horizontal
        size_range = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.marker_size, (80, 350), pop=True)
        cutoff_param_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.cutoff_param_dict, {}, pop=True)
        supplementary_text_format_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.supplementary_text_format_dict, {}, pop=True)

        bar_column_width = 0.2
        bar_edge = 0.09

        total_group_num = len(data_dict)
        stem_value_dict = {}
        head_value_list = []
        color_list = []

        for label, (stem_value, head_value) in data_dict.items():
            stem_value_dict[label] = np.array(stem_value)
            head_value_list.append(head_value)
            color_list.append(color_dict[label])
        head_value_array = np.array(head_value_list)
        max_head_value = head_value_array.max()
        min_head_value = head_value_array.min()
        min_size, max_size = size_range
        if max_head_value > min_head_value:
            size_array = min_size + (head_value_array - min_head_value) / (max_head_value - min_head_value) * (max_size - min_size)
        else:
            mean_size = (min_size + max_size) / 2
            size_array = np.full(total_group_num, mean_size)

        raw_data_scatter_param_dict = {
            ParameterName.random: False,
            ParameterName.marker_size: size_array,
            ParameterName.z_order: 10,
            ParameterName.edge_width: 0.8,
            ParameterName.edge_color: ColorConfig.white_color,
        }

        new_figure_data_parameter_dict = {
            **figure_data_parameter_dict,
            ParameterName.ax_bottom_left_list: [Vector(0, 0)],
            ParameterName.ax_size_list: [Vector(1, 1)],
            ParameterName.color_dict: color_dict,
            ParameterName.data_nested_list: [(stem_value_dict, None, stem_value_dict)],
            ParameterName.array_len_list: [1],
            ParameterName.figure_config_dict: {
                ParameterName.column_width: bar_column_width,
                ParameterName.edge: bar_edge,
                ParameterName.scatter_param_dict: raw_data_scatter_param_dict,
                ParameterName.supplementary_text_format_dict: {
                    **DataFigureConfig.common_supplementary_text_config_dict,
                    ParameterName.font_size: 15,
                    ParameterName.width: 0.1,
                    ParameterName.height: 0.02,
                    **supplementary_text_format_dict,
                },
                ParameterName.cutoff_param_dict: {
                    ParameterName.edge_width: 1,
                    ParameterName.edge_color: ColorConfig.gray,
                    **cutoff_param_dict,
                },
                ParameterName.x_label_format_dict: {ParameterName.font_size: 12},
                ParameterName.y_label_format_dict: {ParameterName.font_size: 12},
                ParameterName.x_tick_label_format_dict: {ParameterName.font_size: 10},
                ParameterName.y_tick_label_format_dict: {ParameterName.font_size: 10},
            },
            ParameterName.x_tick_labels_list: [x_tick_labels],
            ParameterName.y_label_list: [y_label],
            ParameterName.x_label_list: [x_label],
            ParameterName.legend: False,
            ParameterName.horizontal: horizontal,
            ParameterName.isolate_color_bar_mode: True,
            **figure_data_parameter_dict,
        }

        super().__init__(new_figure_data_parameter_dict, bottom_left, size, **kwargs)

    def draw(self, fig=None, parent_ax=None, parent_transformation=None):
        ax_and_transform_list = super().draw(fig, parent_ax, parent_transformation)
        return ax_and_transform_list


class CompleteLollipopFigure(CompositeFigure):
    """带标题的完整棒棒糖图。"""

    height_to_width_ratio = 0.8
    legend_height = 0.06
    title_height = 0.05

    def __init__(self, figure_data_parameter_dict, **kwargs):
        total_width = self.total_width
        figure_height = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_height, 0.8, pop=True)
        figure_title = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_title, '', pop=True)

        figure_bottom_line = 0.05
        figure_left = 0.1
        figure_width = total_width - figure_left
        title_height = self.title_height
        title_center_x = figure_left + figure_width / 2
        title_center_y = figure_height + title_height / 2
        total_height = figure_height + self.title_height
        self.total_height = total_height
        self.height_to_width_ratio = total_height / total_width

        lollipop_config_dict = {
            ParameterName.bottom_left: Vector(figure_left, figure_bottom_line),
            ParameterName.size: Vector(figure_width, figure_height),
            ParameterName.figure_data_parameter_dict: {
                **figure_data_parameter_dict,
            },
        }

        main_title_text_config_dict = {
            **CommonElementConfig.common_text_config,
            ParameterName.font_weight: FontWeight.bold,
            ParameterName.text_box: False,
            ParameterName.string: figure_title,
            ParameterName.font_size: 17,
            ParameterName.width: total_width,
            ParameterName.height: title_height,
            ParameterName.center: Vector(title_center_x, title_center_y),
        }

        subfigure_element_dict = {
            'lollipop': {
                'lollipop': LollipopDataFigure(**lollipop_config_dict)},
            'title': {
                'title': TextBox(**main_title_text_config_dict)},
        }
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)