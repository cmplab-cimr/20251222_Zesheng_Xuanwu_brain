# Figure Plotting Package - 科研图表绑制包

这是一个面向科研论文的Python图表绑制工具包，专门用于生成代谢通量分析（MFA）相关的高质量学术图表。该包使用matplotlib作为底层绑制引擎，提供了一套面向对象的API。

## 目录结构

```
figure_plotting_package/
├── __init__.py              # 包入口，导出统一的 Elements 类
├── main.py                  # 主入口函数 draw_figure()
├── elements.py              # 统一的元素注册类（合并所有子模块元素）
├── common/                  # 通用工具和基础类
│   ├── classes.py           # 基础类：Color, Vector, LineStyle, TransformDict 等
│   ├── config.py            # 参数名称常量定义 ParameterName
│   ├── color.py             # 颜色相关工具
│   ├── core_plotting_functions.py  # 核心绑图函数（colormap、axis设置等）
│   ├── common_functions.py  # 通用辅助函数
│   ├── light_weight_functions.py   # 轻量级工具函数
│   ├── common_figure_materials.py  # 通用图表素材
│   └── third_party_packages.py     # 第三方包导入封装
├── basic_shape_elements/    # 基础形状元素
│   ├── elements.py          # 元素注册
│   ├── shapes.py            # 基础形状：Circle, Rectangle, Ellipse, Polygon 等
│   ├── arrow.py             # 各类箭头：Arrow, ChevronArrow, ArcArrow 等
│   ├── basic_region.py      # 区域（Region）类
│   ├── composite_figure_and_axes.py  # 复合图表：Figure, Subfigure, DataFigureAxes
│   ├── modified_text.py     # 文本框 TextBox
│   └── legend.py            # 图例生成器
├── data_figure/             # 数据图表
│   ├── elements.py          # 数据图表元素注册
│   ├── data_figure.py       # 基础数据图表类 DataFigure
│   ├── scatter_data_figure.py    # 散点图
│   ├── bar_data_figure.py        # 柱状图、MID对比图、通量误差图
│   ├── heatmap_data_figure.py    # 热图、ColorBar
│   ├── violin_box_data_figure.py # 小提琴图/箱线图
│   └── histogram_data_figure.py  # 直方图
├── metabolic_network/       # 代谢网络图
│   ├── elements.py          # 网络元素注册
│   ├── metabolic_network.py # 网络图主类 MetabolicNetwork
│   ├── complex_metabolic_network_figure.py  # 复杂网络图组合
│   ├── metabolic_network_contents/   # 网络内容定义
│   │   ├── metabolite.py    # 代谢物
│   │   ├── reaction.py      # 反应
│   │   └── subnetwork.py    # 子网络
│   ├── metabolic_network_elements/    # 网络可视化元素
│   │   ├── metabolite_element.py
│   │   ├── reaction_element.py
│   │   └── subnetwork_element.py
│   └── layout_generator_functions/    # 布局生成函数
├── diagrams/                # 示意图
│   ├── elements.py          # 示意图元素注册
│   ├── mid_diagram.py       # MID（质量同位素分布）图
│   ├── network_diagram.py   # 网络示意图
│   ├── axis_diagram.py      # 坐标轴示意图
│   ├── carbon_backbone.py   # 碳骨架图
│   ├── simple_algorithm_diagram.py  # 算法流程图
│   └── object_diagrams/     # 生物对象示意图（小鼠、器官等SVG）
│       ├── mice.py          # 小鼠
│       ├── human.py         # 人类
│       ├── fruit_fly.py     # 果蝇
│       └── organs.py        # 器官
```

## 核心类层次

### 1. Elements 统一入口
所有图表元素通过 `Elements` 类统一访问：

```python
from figure_plotting_package import Elements

# 基础形状
Elements.Circle, Elements.Rectangle, Elements.Ellipse, Elements.Polygon
Elements.Arrow, Elements.ChevronArrow, Elements.ArcArrow

# 文本和区域
Elements.TextBox, Elements.Region

# 图表组合
Elements.Figure, Elements.Subfigure, Elements.DataFigureAxes

# 数据图表
Elements.DataFigure, Elements.ScatterDataFigure, Elements.BarDataFigure
Elements.HeatmapDataFigure, Elements.ViolinBoxDataFigure

# 代谢网络
Elements.Metabolite, Elements.Reaction, Elements.Subnetwork
Elements.MetabolicNetwork, Elements.MetabolicNetworkWithLegend

# 示意图
Elements.MIDDiagram, Elements.NetworkDiagram, Elements.AxisDiagram
Elements.Mice, Elements.Human, Elements.FruitFly  # 生物对象
```

### 2. 基础类（common/classes.py）

- **Color**: 颜色类，支持RGB/RGBA多种初始化方式
- **Vector**: 2D向量类，支持几何运算
- **LineStyle, JoinStyle**: 线型和连接样式枚举
- **HorizontalAlignment, VerticalAlignment**: 文本对齐枚举
- **FontStyle, FontWeight**: 字体样式和粗细枚举
- **LinearMapper, SegmentedLinearMappers**: 线性映射器（用于数值到视觉属性的转换）

### 3. 图表基类

- **BasicShape**: 所有形状的基类
- **Region**: 可包含子元素的图形区域
- **CompositeFigure**: 复合图表容器
- **Figure**: 顶级图表类
- **Subfigure**: 子图
- **DataFigureAxes**: 数据绑图坐标轴

## 使用示例

### 绑制基础图表

```python
from figure_plotting_package import Elements, draw_figure

# 创建图表
fig = Elements.Figure(
    figure_label='example',
    # ... 其他参数
)

# 绑制并保存
draw_figure(fig, output_svg=True)
```

### 创建代谢网络图

```python
# 定义代谢物
metabolites = Elements.MetaboliteList(...)
reactions = Elements.ReactionList(...)

# 创建网络图
network = Elements.MetabolicNetwork(
    metabolite_list=metabolites,
    reaction_list=reactions,
    # ... 布局和样式参数
)
```

### 数据图表

```python
# 散点图
scatter = Elements.ScatterDataFigure(
    figure_data=...,
    x_label='X轴',
    y_label='Y轴',
)

# 柱状图
bar = Elements.BarDataFigure(
    figure_data=...,
    # ...
)
```

## 参数配置

所有图表元素使用字典式配置，常用参数名在 `common/config.py` 的 `ParameterName` 类中定义：

- **形状参数**: `center`, `size`, `width`, `height`, `color`, `alpha`, `edge_width`, `edge_color`
- **文本参数**: `text`, `font`, `font_size`, `font_color`, `font_weight`, `angle`
- **箭头参数**: `tail`, `head`, `radius`, `head_width`, `stem_width`
- **图表参数**: `x_lim`, `y_lim`, `x_label`, `y_label`, `x_ticks`, `y_ticks`

## 设计原则

1. **面向对象**: 所有图形元素都是对象，支持嵌套和组合
2. **配置字典化**: 使用字典参数进行样式配置，便于批量管理
3. **分层继承**: Elements 类继承自 BasicElements, MetabolicNetworkElements, DataFigureElements, DiagramElements
4. **坐标系统**: 使用左下角为原点的坐标系统，单位为点（points）

## 依赖

- numpy
- matplotlib
- （其他第三方库见 third_party_packages.py）

## 输出格式

支持 PDF 和 SVG 矢量格式输出，适合科研论文发表要求。