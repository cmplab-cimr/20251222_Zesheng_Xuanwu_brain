# data_figure 模块

`figure_plotting_package` 的数据图实现，包含 `HeatmapDataFigure`、`ColorBarDataFigure`、`BarDataFigure`、`ScatterDataFigure`、`ViolinBoxDataFigure`、`HistogramDataFigure` 等。

## 关键陷阱

### 1. ColorBarDataFigure 的 `x_ticks_list` 必须用数据值

`ParameterName.x_ticks_list` 传入的是**实际数据坐标**（如 `np.linspace(vmin, vmax, 5)`），而非像素/索引位置（如 `np.linspace(0, 255, 5)`）。

- `cbar_plotting()` 内部使用 `mapped_image`（`ScalarMappable`）的 norm/cmap 将数据值映射到颜色，tick 位置对应数据空间而非像素空间
- 传入像素位置会导致 tick 全部挤到图像左角（因为像素值 0~255 在数据空间中被映射到极窄范围）

```python
# 正确
tick_val = np.linspace(vmin, vmax, 5)
cbar_param = {
    ParameterName.x_ticks_list: tick_val,           # 数据值
    ParameterName.x_tick_labels_list: [f'{v:.1f}' for v in tick_val],
    ...
}

# 错误 — tick 会被压缩
cbar_param = {
    ParameterName.x_ticks_list: np.linspace(0, 255, 5),  # 像素位置！
    ...
}
```

### 2. ColorBarDataFigure 作为独立元素时必须传 `cmap_mapper`

`ColorBarDataFigure.draw()` 需要 `mapped_image` 参数来渲染颜色条。当作为独立元素（非 HeatmapDataFigure 的 cbar）使用时，`mapped_image` 默认为 `None`，此时会回退到 `self.cmap_mapper`。因此必须在 `figure_data_parameter_dict` 中提供：

```python
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

mapper = ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap)
cbar_param = {
    ParameterName.cmap_mapper: mapper,
    ...
}
```

### 3. `hidden_axis_tick` vs `hidden_axis_boundary`

两个参数独立控制：
- `hidden_axis_boundary: True` → `axis_format_dict[edge_width] = 0`（隐藏轴线）
- `hidden_axis_tick: True` → `edge_width = 0` 且 `axis_tick_length = 0`（隐藏刻度标记）
- 需要边框但无刻度：`hidden_axis_boundary: False` + `hidden_axis_tick: True`

### 4. HeatmapDataFigure 的 `data_lim_pair` 必须对称

对于发散色图，`data_lim_pair = (-vmax, vmax)` 必须正负对称，否则白色不对应 0 值。