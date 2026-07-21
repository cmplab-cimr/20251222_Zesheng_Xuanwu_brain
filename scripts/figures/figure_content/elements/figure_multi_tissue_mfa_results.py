from pathlib import Path
import sys
import io
from scipy.stats import mannwhitneyu

# Set stdout encoding to UTF-8 for Chinese characters
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 添加项目根目录到 sys.path
_project_root = Path(__file__).parent.parent.parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

from ..common.config import (
    np, ParameterName, Figure as GeneralFigure, Vector,
    diet_color_dict, diet_display_dict, pickle_load)
from ..common.elements import Elements
from common_data.experimental_data.config import Keyword


# 数据目录
_data_dir = 'common_data/experimental_data/mfa_result/raw_data'


def _parse_flux_expression(flux_expr):
    """
    解析flux表达式，支持比例和reverse参数

    语法：
    - 'flux_name' - 普通flux，自动计算净通量（正-逆）
    - 'flux_name:reverse' - 净通量（逆-正）
    - 'flux_A / flux_B' - 比例计算
    - 'flux_A:reverse / flux_B' - 组合使用

    Args:
        flux_expr: flux表达式字符串

    Returns:
        tuple: (is_ratio, flux_spec_list)
            is_ratio: bool，是否是比例计算
            flux_spec_list: list of dict，每个dict包含 {'name': str, 'reverse': bool}
    """
    # 检查是否包含比例符号
    if ' / ' in flux_expr:
        parts = flux_expr.split(' / ')
        flux_spec_list = []
        for part in parts:
            part = part.strip()
            if part.endswith(':reverse'):
                flux_spec_list.append({'name': part[:-8], 'reverse': True})
            else:
                flux_spec_list.append({'name': part, 'reverse': False})
        return True, flux_spec_list
    else:
        # 单个flux
        if flux_expr.endswith(':reverse'):
            return False, [{'name': flux_expr[:-8], 'reverse': True}]
        else:
            return False, [{'name': flux_expr, 'reverse': False}]


def _get_single_flux_data(mouse_data, flux_spec):
    """
    获取单个flux的数据（处理净通量计算）

    Args:
        mouse_data: 小鼠数据字典 {flux_name: array}
        flux_spec: dict，包含 {'name': str, 'reverse': bool}

    Returns:
        np.array: flux数据，如果不存在返回None
    """
    flux_name = flux_spec['name']
    reverse = flux_spec['reverse']

    if flux_name not in mouse_data:
        return None

    forward_flux = mouse_data[flux_name]
    reverse_flux_name = f'{flux_name}__R'

    if reverse_flux_name in mouse_data:
        reverse_flux = mouse_data[reverse_flux_name]
        if reverse:
            # 逆反应减正反应
            net_flux = reverse_flux - forward_flux
        else:
            # 正反应减逆反应（默认）
            net_flux = forward_flux - reverse_flux
        return net_flux
    else:
        return forward_flux


def _filter_middle_percent(data, lower_percentile=25, upper_percentile=75):
    """
    过滤数据，只保留中间50%的数据（舍去两端各25%分位）

    Args:
        data: np.array，一维数据数组

    Returns:
        np.array: 过滤后的数据
    """
    if len(data) == 0:
        return data

    lower_bound = np.percentile(data, lower_percentile)
    upper_bound = np.percentile(data, upper_percentile)

    # 保留在25%到75%分位之间的数据
    filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    return filtered_data


def _statistical_test(flux_expr, diet1, data1, diet2, data2):
    # 进行两两比较
    statistic, p_value = mannwhitneyu(data1, data2, alternative='two-sided')
    display_name1 = diet_display_dict.get(diet1, diet1)
    display_name2 = diet_display_dict.get(diet2, diet2)
    print(f'{flux_expr} - Mann-Whitney U检验 ({display_name1} vs {display_name2}):')
    print(
        f'  {display_name1}组: n={len(data1) // 75}只小鼠, mean={np.mean(data1):.4f}, std={np.std(data1):.4f}')
    print(
        f'  {display_name2}组: n={len(data2) // 75}只小鼠, mean={np.mean(data2):.4f}, std={np.std(data2):.4f}')
    print(f'  U统计量: {statistic:.2f}, p值: {p_value:.2e}')
    if p_value < 0.05:
        print(f'  结论: 两组间存在显著差异 (p < 0.05)')
    else:
        print(f'  结论: 两组间无显著差异 (p >= 0.05)')
    print()


def _prepare_flux_comparison_data_by_diet(
        flux_data_dict, diet_list, flux_name_list, original_diet_color_dict,
        combine_bio_replicates=True):
    """
    将原始flux数据转换为按饮食类型分组的格式（用于比较不同饮食的同一flux）

    支持以下flux表达式语法：
    - 'flux_name' - 普通flux，自动计算净通量（正-逆）
    - 'flux_name:reverse' - 净通量（逆-正）
    - 'flux_A / flux_B' - 比例计算
    - 'flux_A:reverse / flux_B:reverse' - 组合使用

    Args:
        flux_data_dict: 原始flux数据字典 {diet_type: {mouse_index: {flux_name: array}}}
        diet_list: 需要包含的饮食类型列表
        flux_name_list: 需要绘制的flux表达式列表
        original_diet_color_dict: 饮食类型到颜色的映射字典
        combine_bio_replicates: bool，是否合并生物学重复
            - True: 同一饮食下的所有小鼠数据合并为一个class
            - False: 每个小鼠作为独立的class（如hfd_1, lfd_1）

    Returns:
        complete_data_dict: {flux_expr: {'': {class_name: np.array}}}
        mouse_color_dict: {class_name: color}
    """
    complete_data_dict = {flux_expr: {} for flux_expr in flux_name_list}
    mouse_color_dict = {}

    for flux_expr in flux_name_list:
        is_ratio, flux_spec_list = _parse_flux_expression(flux_expr)

        # 确保空字符串group存在
        if '' not in complete_data_dict[flux_expr]:
            complete_data_dict[flux_expr][''] = {}

        for diet_type in diet_list:
            if diet_type not in flux_data_dict:
                continue

            if combine_bio_replicates:
                # 合并模式：同一饮食下所有小鼠数据合并
                all_flux_data = []
                for mouse_index, mouse_data in flux_data_dict[diet_type].items():
                    if is_ratio:
                        numerator_data = _get_single_flux_data(mouse_data, flux_spec_list[0])
                        denominator_data = _get_single_flux_data(mouse_data, flux_spec_list[1])
                        if numerator_data is not None and denominator_data is not None:
                            ratio_data = numerator_data / (denominator_data + 1e-10)
                            all_flux_data.append(ratio_data)
                    else:
                        flux_data = _get_single_flux_data(mouse_data, flux_spec_list[0])
                        if flux_data is not None:
                            all_flux_data.append(flux_data)

                if all_flux_data:
                    # 合并所有数据
                    combined_data = np.concatenate(all_flux_data)
                    # 过滤中间50%数据
                    filtered_data = _filter_middle_percent(combined_data)
                    class_name = diet_type
                    complete_data_dict[flux_expr][''][class_name] = filtered_data
                    mouse_color_dict[class_name] = original_diet_color_dict[diet_type]
                else:
                    raise ValueError()

            else:
                # 分离模式：每个小鼠作为独立class
                for mouse_index, mouse_data in flux_data_dict[diet_type].items():
                    class_name = f'{diet_type}_{mouse_index}'

                    if is_ratio:
                        numerator_data = _get_single_flux_data(mouse_data, flux_spec_list[0])
                        denominator_data = _get_single_flux_data(mouse_data, flux_spec_list[1])
                        if numerator_data is not None and denominator_data is not None:
                            ratio_data = numerator_data / (denominator_data + 1e-10)
                            # 过滤中间50%数据
                            filtered_data = _filter_middle_percent(ratio_data)
                            complete_data_dict[flux_expr][''][class_name] = filtered_data
                        else:
                            raise ValueError()
                    else:
                        flux_data = _get_single_flux_data(mouse_data, flux_spec_list[0])
                        if flux_data is not None:
                            # 过滤中间50%数据
                            filtered_data = _filter_middle_percent(flux_data)
                            complete_data_dict[flux_expr][''][class_name] = filtered_data
                        else:
                            raise ValueError()

                    # 根据饮食类型分配颜色
                    if class_name not in mouse_color_dict and class_name in complete_data_dict[flux_expr].get('', {}):
                        mouse_color_dict[class_name] = original_diet_color_dict[diet_type]
        if combine_bio_replicates:
            diet1, diet2 = diet_list
            data1 = complete_data_dict[flux_expr][''][diet1]
            data2 = complete_data_dict[flux_expr][''][diet2]
            _statistical_test(flux_expr, diet1, data1, diet2, data2)
    return complete_data_dict, mouse_color_dict


def flux_comparison_data_generator(
        flux_name_nested_list, original_diet_color_dict, combine_bio_replicates=True, diet_list=None):
    """
    生成flux比较图的数据

    Args:
        flux_name_nested_list: 二维列表，每行是一组需要比较的flux名称
        original_diet_color_dict: 饮食类型到颜色的映射字典
        combine_bio_replicates: bool，是否合并生物学重复
        diet_list: 需要包含的饮食类型列表，默认为 [Keyword.ctrl, Keyword.hfd]

    Returns:
        complete_data_dict: 用于FluxComparisonGridBoxDataFigure的数据字典
        mouse_color_dict: {class_name: color}
    """
    if diet_list is None:
        diet_list = [Keyword.ctrl, Keyword.hfd]

    # 加载数据
    pickle_dict = pickle_load(_data_dir)
    flux_data_dict = pickle_dict['flux_data_dict']

    # 展平flux_name_nested_list获取所有需要处理的flux名称
    all_flux_names = []
    for row in flux_name_nested_list:
        all_flux_names.extend(row)

    # 准备数据
    complete_data_dict, mouse_color_dict = _prepare_flux_comparison_data_by_diet(
        flux_data_dict, diet_list, all_flux_names, original_diet_color_dict,
        combine_bio_replicates=combine_bio_replicates)

    return complete_data_dict, mouse_color_dict


class Figure(GeneralFigure):
    figure_output_direct = f'{GeneralFigure.figure_output_direct}/mfa_results'
    figure_label = 'figure_multi_tissue_mfa_results'
    figure_size = Vector(8, 8)

    figure_label_data_name_dict = {
        'mfa_result_combine': ('figure_mfa_result_combine', True),
        'mfa_result_combine_2': ('figure_mfa_result_combine_2', True),
        'mfa_result_separate': ('figure_mfa_result_separate', False),
    }

    def __init__(self):
        # current_label = 'mfa_result_combine'
        current_label = 'mfa_result_combine_2'
        # current_label = 'mfa_result_separate'

        (figure_name, combine_bio_replicates) = self.figure_label_data_name_dict[current_label]

        # 定义要绘制的flux
        # 根据实验结果设计，验证以下假设：
        # 1. 肝脏糖异生流量显著升高
        # 2. 肝脏脂肪酸β-氧化流量无显著变化，TCA循环总流量上升
        # 3. 骨骼肌葡萄糖氧化降低，脂肪酸氧化升高
        # 4. WAT代谢变化
        if current_label == 'mfa_result_combine_2':
            flux_name_nested_list = [
                # 第一行：心脏
                ['liver--G6Pase_c', 'liver--CS_m', 'liver--FA_usage', 'liver--FA_usage / liver--CS_m', ],
                # 第二行：肝脏糖酵解vs TCA循环
                ['muscle--PDH_m', 'muscle--CS_m', 'muscle--ICDH_m',  'muscle--FA_usage / muscle--CS_m', ],
                # 第三行：肝脏脂肪酸氧化（预期无显著变化）
                ['white_adipose--lipid_mobilization_c', 'white_adipose--CS_m',
                 'white_adipose--FA_usage', 'white_adipose--FA_usage / white_adipose--CS_m', ],
            ]
        else:
            flux_name_nested_list = [
                # 第一行：肝脏糖异生（预期HFD显著升高）
                ['liver--G6Pase_c', 'liver--CS_m', 'liver--FA_usage', 'liver--FA_usage / liver--CS_m', ],
                # 第二行：肝脏糖酵解vs TCA循环
                ['muscle--GLC_exchange', 'muscle--GLYCOGENOLYSIS_c', 'muscle--PDH_m', 'muscle--CS_m', ],
                # 第三行：肝脏脂肪酸氧化（预期无显著变化）
                ['muscle--ICDH_m', 'muscle--G6PDH2R_PGL_GND_c', 'white_adipose--CS_m', 'white_adipose--lipid_mobilization_c'],
            ]

        # flux显示名称映射（用于在图表上显示更友好的名称）
        display_flux_dict = {
            # 肝脏糖异生
            # 'liver--PEPCK_c': 'PEPCK\n(Liver)',
            'liver--G6Pase_c': 'G6Pase\n(Liver)',
            # 'liver--FBPase_c': 'FBPase\n(Liver)',
            # 肝脏糖酵解和TCA
            'liver--PYK_c': 'PYK\n(Liver)',
            'liver--CS_m': 'CS\n(Liver)',
            'liver--ICDH_m': 'ICDH\n(Liver)',
            'liver--AKGD_m': 'αKGDH\n(Liver)',
            # 肝脏脂肪酸氧化
            'liver--FA_usage': 'FA oxidation\n(Liver)',
            'liver--PDH_m': 'PDH\n(Liver)',
            # 肌肉葡萄糖氧化
            'muscle--PDH_m': 'PDH\n(Muscle)',
            'muscle--CS_m': 'CS\n(Muscle)',
            'muscle--PDH_m / muscle--CS_m': 'PDH/CS\n(Muscle)',
            # 肌肉脂肪酸氧化
            'muscle--FA_usage': 'FA oxidation\n(Muscle)',
            'muscle--FA_exchange': 'FA exchange\n(Muscle)',
            'muscle--FA_usage / muscle--CS_m': 'FAox/CS\n(Muscle)',
            # WAT
            'white_adipose--FA_exchange': 'FA exchange\n(WAT)',
            'white_adipose--lipid_mobilization_c': 'Lipid mobil.\n(WAT)',
            'white_adipose--CS_m': 'CS\n(WAT)',
        }
        display_flux_dict = {flux_name: '' for row in flux_name_nested_list for flux_name in row}

        # 加载并准备数据
        complete_data_dict, mouse_color_dict = flux_comparison_data_generator(
            flux_name_nested_list, diet_color_dict,
            combine_bio_replicates=combine_bio_replicates)

        figure_data_parameter_dict = {
            ParameterName.figure_data: complete_data_dict,
            ParameterName.name_dict: diet_display_dict,
            ParameterName.legend_color_dict: mouse_color_dict,
            ParameterName.color_dict: mouse_color_dict,
            ParameterName.flux_name_list: flux_name_nested_list,
            ParameterName.display_flux_name_dict: display_flux_dict,
            ParameterName.common_x_label: None,
            ParameterName.common_y_label: None,
            ParameterName.low_height: True,
            ParameterName.legend: False,
            ParameterName.figure_config_dict: {
                # ParameterName.column_width: 0.85,
                ParameterName.error_bar_param_dict: {
                    ParameterName.cap_size: 0.5,
                    ParameterName.edge_width: 0.1,
                },
                ParameterName.box_violin_config_dict: {
                    # ParameterName.alpha: 1,
                    ParameterName.column_width: 0.6,
                    ParameterName.gap_inside_column: 0.5,
                }
            }
        }
        flux_comparison_obj = Elements.FluxComparisonGridBoxDataFigure(**{
            ParameterName.bottom_left_offset: Vector(0.1, 0.05),
            ParameterName.scale: 0.8,
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
        })
        content_obj_list = [flux_comparison_obj]

        super().__init__(figure_name, other_obj_list=content_obj_list)
