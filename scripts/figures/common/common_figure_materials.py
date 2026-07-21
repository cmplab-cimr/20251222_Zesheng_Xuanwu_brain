from .config import plt, mcolors, ColorConfig, Color
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np


class ColorMapName(object):
    diet_color_map = 'diet_color_map'
    tissue_color_map = 'tissue_color_map'
    state_color_map = 'state_color_map'
    blue_white_orange_cmap = 'blue_white_orange_cmap'


set_2_color = plt.get_cmap('Set2').colors
pastel_1_color = plt.get_cmap('Pastel1').colors


class CommonColorDict(object):
    """Figure 4 代谢组学主色调 + 共用颜色"""

    # 语义色（灰+青双色，低饱和极简）
    # 灰: AD / 对照 / 下调 / 负向
    gray = Color(128, 128, 128)       # '#808080'
    # 青: Semaglutide / 处理 / 上调 / 正向
    teal = Color(6, 125, 150)         # '#067d96'
    # 浅灰: 不显著 / 中性背景
    light_gray = Color(217, 217, 217) # '#D9D9D9'
    white = Color(255, 255, 255)      # '#FFFFFF'

    # 脂质强调色
    lipid_highlight = Color(109, 164, 109)  # '#6DA46D'

    # RdBu 发散色图端点
    rdbu_negative = Color(33, 102, 172)   # '#2166AC' 蓝(负/耗竭)
    rdbu_positive = Color(178, 24, 43)    # '#B2182B' 红(正/富集)

    # 分组语义
    up = teal
    down = gray
    insig = light_gray


class CommonFigureString(object):
    pass


# ============================================================
# 发散色图
# ============================================================

# 灰 → 白 → 青，中点(0)=白
_gray_teal_cmap = None


def gray_teal_diverging_cmap():
    global _gray_teal_cmap
    if _gray_teal_cmap is None:
        _gray_teal_cmap = LinearSegmentedColormap.from_list(
            'gray_teal', ['#808080', '#FFFFFF', '#067d96'])
    return _gray_teal_cmap


# 蓝(负/耗竭) → 白(0) → 红(正/富集)
_rdbu_cmap = None


def rdbu_diverging_cmap():
    global _rdbu_cmap
    if _rdbu_cmap is None:
        _rdbu_cmap = LinearSegmentedColormap.from_list(
            'rdbu', ['#2166AC', '#FFFFFF', '#B2182B'])
    return _rdbu_cmap


# ============================================================
# 脂质类判断
# ============================================================

LIPID_EXACT = {'GP', 'SP', 'GL', 'FA', 'Bile acids',
               'Cholesterol and esters', 'Steroids', 'Sterol lipids'}
LIPID_KEYWORDS = ['glycerophospho', 'sphingo', 'glycerolipid',
                  'cholesterol', 'sterol', 'steroid', 'fatty acid', 'bile acid', 'carnitine', '脂']


def is_lipid_class(class_label):
    """判断 Class I 标签是否为脂质类。"""
    if class_label is None:
        return False
    label = str(class_label).strip()
    if label in LIPID_EXACT:
        return True
    low = label.lower()
    return any(kw.lower() in low for kw in LIPID_KEYWORDS)


def lipid_class_color(class_label):
    """脂质类统一 lipid_highlight，非脂质 light_gray。"""
    return CommonColorDict.lipid_highlight if is_lipid_class(class_label) else CommonColorDict.light_gray


# ============================================================
# 工具函数
# ============================================================

def fdr_star(fdr):
    """FDR/q 值 → 星号字符串: <0.001 '***', <0.01 '**', <0.05 '*', 否则 ''。"""
    if fdr < 0.001:
        return '***'
    if fdr < 0.01:
        return '**'
    if fdr < 0.05:
        return '*'
    return ''


def make_colorbar(cmap, vmin, vmax, bottom_left, size, label='', n_ticks=5,
                  orientation='horizontal'):
    """创建框架标准 ColorBarDataFigure。

    orientation: 'horizontal' (底部) 或 'vertical' (右侧)。
    """
    from .config import ParameterName, Vector
    from ..figure_plotting_package.data_figure.heatmap_data_figure import ColorBarDataFigure

    mapper = ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap)
    tick_val = np.linspace(vmin, vmax, n_ticks)
    tick_label = [f'{v:.1f}' for v in tick_val]
    axis_label = label if label else None
    cbar_param = {
        ParameterName.cmap_mapper: mapper,
        ParameterName.cbar_orientation: orientation,
        ParameterName.hidden_axis_boundary: False,
        ParameterName.figure_config_dict: {},
    }
    if orientation == ParameterName.vertical:
        cbar_param[ParameterName.right_y_ticks_list] = tick_val
        cbar_param[ParameterName.right_y_tick_labels_list] = tick_label
        cbar_param[ParameterName.right_y_label_list] = axis_label
    else:
        cbar_param[ParameterName.x_ticks_list] = tick_val
        cbar_param[ParameterName.x_tick_labels_list] = tick_label
        cbar_param[ParameterName.x_label_list] = axis_label
    return ColorBarDataFigure(cbar_param, bottom_left=bottom_left, size=size)
