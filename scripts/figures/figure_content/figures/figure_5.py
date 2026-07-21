# -*- coding: utf-8 -*-
"""Figure 5: 肝脏–星胶–神经元 三区室 ¹³C 代谢流分析

Panel A: 三区室模型示意图
Panel B: 实测 vs 预测 MID (拟合优度)
Panel C: 全身 flux 网络 (本轮留空占位)
Panel D: acetyl-CoA 三层 Sankey (来源→枢纽→去向)
Panel E: 关键 flux 柱+单点+p
Panel F: flux vs 认知相关 (虚拟映射占位, 代谢流与代谢组不同小鼠)

数据: scripts/common_data/raw_data/mfa_analysis/data_20251115_ad_sema_neuron_astrocyte_cortex_model
       flux_raw_data.xlsx (12 样本 × 50 次拟合 × 522 反应) + mid_raw_data.xlsx
统计: 每样本 50 次拟合取中位数 → 6v6 Mann-Whitney + log2FC (不画 CI)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ...common.config import (
    ParameterName, Figure as GeneralFigure, Vector, CompositeFigure,
    default_parameter_extract, calculate_center_bottom_offset)
from ...common.elements import Elements
from ...common.common_figure_materials import CommonColorDict as CCD
from ...figure_plotting_package.basic_shape_elements.composite_figure_and_axes import DataFigureAxes
from scripts.src.common_parallel_solver.config import pickle_load
from scripts.src.metabolomics_analysis.metabolomics_main_functions import _cache_path
from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword

Subfigure = Elements.Subfigure

# 分组颜色: AD(ctrl)=灰, Semaglutide=青
GROUP_COLORS = {Keyword.ad: CCD.gray, Keyword.sema: CCD.teal}

# ---- 数据路径 ----
_FLUX_DIR = (Path(__file__).parents[3] / 'common_data' / 'raw_data' / 'mfa_analysis'
             / 'data_20251115_ad_sema_neuron_astrocyte_cortex_model')
_FLUX_XLSX = _FLUX_DIR / 'flux_raw_data.xlsx'
_MID_XLSX = _FLUX_DIR / 'mid_raw_data.xlsx'

# ---- 样本 sheet (12: ctrl 1-6 + semaglutide 1-6) ----
_CTRL_SHEETS = [f'glucose_ad_ctrl_{i}' for i in range(1, 7)]
_SEMA_SHEETS = [f'glucose_ad_semaglutide_{i}' for i in range(1, 7)]
_ALL_SHEETS = _CTRL_SHEETS + _SEMA_SHEETS
_SAMPLE_GROUP = {s: Keyword.ad for s in _CTRL_SHEETS} | {s: Keyword.sema for s in _SEMA_SHEETS}

# flux 通量 sheet → flux 小鼠 ID 由 behavioral_data_loader._flux_sheet_to_mouse_map 动态构建
# (经器官编号中转: 编号对应.xlsx 小鼠ID→器官编号 + specific_data_parameters A-L→condition_index)

# 区室
COMPARTMENTS = ['liver', 'astrocyte', 'neuron', 'serum']


# ============================================================
# 通用 helper (同 figure_4 模式)
# ============================================================

def move_center(element_obj, scale, subfigure_size, calculate_center_param_dict=None, extra_offset=None):
    if calculate_center_param_dict is None:
        calculate_center_param_dict = {}
    if extra_offset is None:
        extra_offset = Vector(0., 0.)
    center = element_obj.calculate_center(element_obj, scale, **calculate_center_param_dict)
    center_bottom_offset = calculate_center_bottom_offset(center, subfigure_size) + extra_offset
    element_obj.move_and_scale(bottom_left_offset=center_bottom_offset)


def _place(figure_class, figure_data_parameter_dict, subfigure_bottom_left, subfigure_size):
    """probe + fit scale + 构造(**kwargs) + move_center。"""
    probe = figure_class(figure_data_parameter_dict=dict(figure_data_parameter_dict))
    scale = float(np.min(subfigure_size / probe.size))
    element = figure_class(**{
        ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
        ParameterName.scale: scale,
        ParameterName.bottom_left_offset: subfigure_bottom_left,
    })
    move_center(element, scale, subfigure_size)
    return element


# ============================================================
# flux / mid 数据加载 + 统计 (内联, 不碰 mfa_analysis junction)
# ============================================================

_FLUX_CACHE: dict = {}
def _load_flux_dict() -> dict:
    """读 12 样本 flux sheet → {sheet: DataFrame(50 行 × 反应列)}。"""
    if not _FLUX_CACHE:
        for s in _ALL_SHEETS:
            _FLUX_CACHE[s] = pd.read_excel(_FLUX_XLSX, sheet_name=s)
    return _FLUX_CACHE


_MEDIAN_CACHE: dict | None = None
def _flux_per_sample_median() -> dict:
    """每样本 50 次拟合取中位数 → {sheet: {reaction: median}}。"""
    global _MEDIAN_CACHE
    if _MEDIAN_CACHE is None:
        flux = _load_flux_dict()
        _MEDIAN_CACHE = {s: df.drop(columns=['Index', 'loss'], errors='ignore')
                            .apply(pd.to_numeric, errors='coerce').median().to_dict()
                         for s, df in flux.items()}
    return _MEDIAN_CACHE


def _flux_stats(reactions) -> pd.DataFrame:
    """每个反应: 6 ctrl vs 6 sema Mann-Whitney + log2FC + 中位数。

    Returns DataFrame [reaction, ctrl_median, sema_median, log2FC, p_mw]
    """
    med = _flux_per_sample_median()
    rows = []
    for r in reactions:
        ctrl = np.array([med[s].get(r, np.nan) for s in _CTRL_SHEETS], dtype=float)
        sema = np.array([med[s].get(r, np.nan) for s in _SEMA_SHEETS], dtype=float)
        ctrl = ctrl[~np.isnan(ctrl)]
        sema = sema[~np.isnan(sema)]
        if len(ctrl) < 2 or len(sema) < 2:
            continue
        cm, sm = float(np.median(ctrl)), float(np.median(sema))
        _, p = stats.mannwhitneyu(ctrl, sema, alternative='two-sided')
        log2fc = float(np.log2(sm / cm)) if (cm > 0 and sm > 0) else np.nan
        rows.append({'reaction': r, 'ctrl_median': cm, 'sema_median': sm,
                     'log2FC': log2fc, 'p_mw': float(p),
                     'ctrl_values': ctrl.tolist(), 'sema_values': sema.tolist()})
    return pd.DataFrame(rows)


_MID_CACHE: dict | None = None
def _load_mid() -> dict:
    """读 mid → {sheet: {metabolite: {mid_index: (exp, pred_mean)}}}。pred = 50 次均值。"""
    global _MID_CACHE
    if _MID_CACHE is None:
        out = {}
        for s in _ALL_SHEETS:
            df = pd.read_excel(_MID_XLSX, sheet_name=s)
            pred_cols = [c for c in df.columns if str(c).isdigit()]
            metab = df[df['metabolite_name'] != 'loss']
            d = {}
            for _, row in metab.iterrows():
                m = str(row['metabolite_name'])
                idx = row['mid_index']
                if pd.isna(idx):
                    continue
                exp = float(row['experiments']) if pd.notna(row['experiments']) else np.nan
                pred = float(pd.to_numeric(row[pred_cols], errors='coerce').mean())
                d.setdefault(m, {})[int(idx)] = (exp, pred)
            out[s] = d
        _MID_CACHE = out
    return _MID_CACHE


# ============================================================
# Panel A: 三区室模型示意图 (TODO: ThreeCompartmentDiagramFigure, task 53)
# ============================================================

class ThreeCompartmentAxes(DataFigureAxes):
    """三区室示意图: Liver / Serum / ISF / Astrocyte / Neuron + ¹³C-glucose 注入 + 主要通路。"""
    def __init__(self, bottom_left, size, **kwargs):
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        nodes = {
            'Liver': (0.10, 0.65), 'Serum': (0.36, 0.65), 'ISF': (0.58, 0.65),
            'Astrocyte': (0.84, 0.78), 'Neuron': (0.84, 0.42),
        }
        for name, (x, y) in nodes.items():
            ax.scatter([x], [y], s=420, c='#EEEEEE', edgecolors='#555555',
                       linewidths=0.6, zorder=3)
            ax.text(x, y, name, fontsize=5, ha='center', va='center',
                    fontfamily='DejaVu Sans', zorder=4)
        edges = [('Liver', 'Serum'), ('Serum', 'ISF'), ('ISF', 'Astrocyte'),
                 ('ISF', 'Neuron'), ('Astrocyte', 'Neuron')]
        for a, b in edges:
            ax.annotate('', xy=nodes[b], xytext=nodes[a],
                        arrowprops=dict(arrowstyle='->', color='#999999', lw=0.9,
                                        shrinkA=13, shrinkB=13), zorder=1)
        # ¹³C-glucose 注入 Serum
        ax.text(0.36, 0.93, '¹³C-glucose', fontsize=5.5, ha='center',
                color='#d85a30', fontfamily='DejaVu Sans')
        ax.annotate('', xy=(0.36, 0.70), xytext=(0.36, 0.88),
                    arrowprops=dict(arrowstyle='->', color='#d85a30', lw=1.2), zorder=1)
        # 模块标注
        for lab, (x, y) in [('Glycolysis / TCA / Lipid syn', (0.10, 0.45)),
                            ('Lactate / Alanine shuttle', (0.47, 0.45)),
                            ('ANLS (MCT lactate)', (0.74, 0.60)),
                            ('Gln–Glu cycle', (0.84, 0.22))]:
            ax.text(x, y, lab, fontsize=4, ha='center', fontfamily='DejaVu Sans',
                    color='#666666')
        ax.set_axis_off()
        return ax, transform


class ThreeCompartmentFigure(CompositeFigure):
    def __init__(self, figure_data_parameter_dict, **kwargs):
        total_width = 1.0
        total_height = 0.55
        axes_obj = ThreeCompartmentAxes(
            Vector(0.02, 0.02), Vector(total_width - 0.04, total_height - 0.04))
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__({'three_comp': {'three_comp': axes_obj}},
                         Vector(0, 0), Vector(total_width, total_height),
                         background=False, **kwargs)


class SubfigureA(Subfigure):
    subfigure_label = 'a'
    subfigure_title = 'three_compartment_diagram'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        element = _place(ThreeCompartmentFigure, {}, subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel B: 实测 vs 预测 MID (TODO: MIDComparisonGridBarWithLegendDataFigure, task 54)
# ============================================================

class MIDGoodnessAxes(DataFigureAxes):
    """实测 vs 预测 MID 散点 (全代谢物全 M+i) + y=x 线 + R²。"""
    def __init__(self, exp_pred, bottom_left, size, **kwargs):
        self.exp = np.array([v[0] for v in exp_pred], dtype=float)
        self.pred = np.array([v[1] for v in exp_pred], dtype=float)
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        valid = np.isfinite(self.exp) & np.isfinite(self.pred)
        x, y = self.exp[valid], self.pred[valid]
        ax.scatter(x, y, s=6, c=CCD.teal, alpha=0.4, edgecolors='none')
        hi = float(max(x.max(), y.max())) * 1.05 if len(x) else 1.0
        ax.plot([0, hi], [0, hi], '--', color='#888888', lw=0.8, zorder=2)
        ax.set_xlim(0, hi)
        ax.set_ylim(0, hi)
        ss_tot = float(np.sum((y - y.mean()) ** 2)) if len(y) else 0
        r2 = 1 - float(np.sum((y - x) ** 2)) / ss_tot if ss_tot > 0 else 0
        ax.text(0.04, 0.96, f'R²={r2:.3f}', transform=ax.transAxes, fontsize=8, va='top')
        ax.set_xlabel('Experimental MID', fontsize=7)
        ax.set_ylabel('Predicted MID', fontsize=7)
        ax.tick_params(labelsize=6)
        return ax, transform


class MIDGoodnessFigure(CompositeFigure):
    def __init__(self, figure_data_parameter_dict, **kwargs):
        exp_pred = default_parameter_extract(figure_data_parameter_dict, 'exp_pred', [], pop=True)
        if not exp_pred:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return
        total_width = 0.7
        total_height = 0.7
        axes_obj = MIDGoodnessAxes(exp_pred, Vector(0.20, 0.15),
                                   Vector(total_width - 0.25, total_height - 0.25))
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__({'mid_go': {'mid_go': axes_obj}},
                         Vector(0, 0), Vector(total_width, total_height), background=False, **kwargs)


class SubfigureB(Subfigure):
    subfigure_label = 'b'
    subfigure_title = 'mid_goodness_of_fit'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        mid = _load_mid()
        exp_pred = [(e, p) for s in _ALL_SHEETS for m in mid[s].values()
                    for (e, p) in m.values()]
        print(f'Panel 5B: {len(exp_pred)} 个 exp/pred MID 点')
        element = _place(MIDGoodnessFigure, {'exp_pred': exp_pred},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel C: 全身 flux 网络 (本轮留空占位)
# ============================================================

class SubfigureC(Subfigure):
    """5C flux 网络图本轮留空占位 (用户决策)。"""
    subfigure_label = 'c'
    subfigure_title = 'flux_network_placeholder'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        print('Panel 5C: flux 网络图本轮留空占位')
        super().__init__({}, subfigure_bottom_left, subfigure_size,
                         subfigure_label=self.subfigure_label,
                         subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel D: acetyl-CoA Sankey (TODO: SankeyDiagramFigure, task 55)
# ============================================================

# acetyl-CoA 当量换算 spec (plan §4): (reaction, label, coef, color); FA C16→8, 酮体×2
ACETYL_SRC_SPEC = [
    ('astrocyte--PDH_m', 'Glucose→PDH', 1, '#4DBBD5'),
    ('astrocyte--FA_usage', 'FA β-oxidation', 8, '#E64B35'),
    ('astrocyte--ACACCOA_lyase_m', 'Ketone body', 2, '#F39B7F'),
    ('astrocyte--ACSS2_m', 'Acetate (ACSS2)', 1, '#91D1C2'),
]
ACETYL_SINK_SPEC = [
    ('astrocyte--ICDH_m', 'TCA oxidation', 1, '#3A7CA5'),
    ('astrocyte--ACITL_c', '→Lipid (ACITL)', 1, '#D85A30'),
]


def _acetyl_coa_flows() -> dict:
    """从 flux 中位数算 ctrl/sema 的 acetyl-CoA 来源/去向当量 (×化学计量系数)。"""
    med = _flux_per_sample_median()
    out = {}
    for grp, sheets in [(Keyword.ad, _CTRL_SHEETS), (Keyword.sema, _SEMA_SHEETS)]:
        src, sink = {}, {}
        for rxn, lab, coef, _ in ACETYL_SRC_SPEC:
            vals = [med[s].get(rxn) for s in sheets if pd.notna(med[s].get(rxn))]
            src[lab] = float(np.median(vals)) * coef if vals else 0.0
        for rxn, lab, coef, _ in ACETYL_SINK_SPEC:
            vals = [med[s].get(rxn) for s in sheets if pd.notna(med[s].get(rxn))]
            sink[lab] = float(np.median(vals)) * coef if vals else 0.0
        out[grp] = {'src': src, 'sink': sink}
    return out


class AcetylCoaSankeyAxes(DataFigureAxes):
    """acetyl-CoA 三层 Sankey (来源→枢纽→去向); ctrl/sema 并排两套。"""
    def __init__(self, ctrl, sema, src_spec, sink_spec, bottom_left, size, **kwargs):
        self.ctrl = ctrl
        self.sema = sema
        self.src_spec = src_spec
        self.sink_spec = sink_spec
        super().__init__(bottom_left, size, **kwargs)

    def _draw_block(self, ax, x0, flows, title):
        from matplotlib.patches import Rectangle
        src_vals = [flows['src'].get(lab, 0) for _, lab, _, _ in self.src_spec]
        sink_vals = [flows['sink'].get(lab, 0) for _, lab, _, _ in self.sink_spec]
        norm = max(sum(src_vals), sum(sink_vals), 1.0)
        H = 0.62
        y0 = 0.18
        yc = y0 + H
        src_ys = []
        for v, (_, lab, _, col) in zip(src_vals, self.src_spec):
            h = v / norm * H
            yc -= h
            ax.add_patch(Rectangle((x0, yc), 0.04, h, facecolor=col, edgecolor='white', lw=0.3))
            ax.text(x0 - 0.01, yc + h / 2, f'{lab} ({v:.0f})', fontsize=3.5,
                    ha='right', va='center', fontfamily='DejaVu Sans')
            src_ys.append((yc, yc + h, col))
        mid_x = x0 + 0.16
        ax.add_patch(Rectangle((mid_x, y0), 0.04, H, facecolor='#CCCCCC', edgecolor='#555555', lw=0.4))
        ax.text(mid_x + 0.02, y0 + H + 0.03, 'acetyl-CoA', fontsize=4.5,
                ha='center', fontfamily='DejaVu Sans', fontweight='bold')
        right_x = x0 + 0.32
        yc = y0 + H
        sink_ys = []
        for v, (_, lab, _, col) in zip(sink_vals, self.sink_spec):
            h = v / norm * H
            yc -= h
            ax.add_patch(Rectangle((right_x, yc), 0.04, h, facecolor=col, edgecolor='white', lw=0.3))
            ax.text(right_x + 0.05, yc + h / 2, f'{lab} ({v:.0f})', fontsize=3.5,
                    ha='left', va='center', fontfamily='DejaVu Sans')
            sink_ys.append((yc, yc + h, col))
        for (y1, y2, col) in src_ys:
            ax.fill([x0 + 0.04, mid_x, mid_x, x0 + 0.04], [y1, y0, y0 + H, y2],
                    facecolor=col, alpha=0.35, edgecolor='none')
        for (y1, y2, col) in sink_ys:
            ax.fill([mid_x + 0.04, right_x, right_x, mid_x + 0.04], [y0, y1, y2, y0 + H],
                    facecolor=col, alpha=0.35, edgecolor='none')
        ax.text(x0 + 0.18, 0.95, title, fontsize=6, ha='center', fontweight='bold',
                fontfamily='DejaVu Sans')

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1)
        self._draw_block(ax, 0.02, self.ctrl, 'AD (ctrl)')
        self._draw_block(ax, 0.55, self.sema, 'Sema')
        ax.set_axis_off()
        return ax, transform


class AcetylCoaSankeyFigure(CompositeFigure):
    def __init__(self, figure_data_parameter_dict, **kwargs):
        ctrl = default_parameter_extract(figure_data_parameter_dict, 'ctrl', {}, pop=True)
        sema = default_parameter_extract(figure_data_parameter_dict, 'sema', {}, pop=True)
        total_width = 1.1
        total_height = 0.65
        axes_obj = AcetylCoaSankeyAxes(
            ctrl, sema, ACETYL_SRC_SPEC, ACETYL_SINK_SPEC,
            Vector(0.02, 0.02), Vector(total_width - 0.04, total_height - 0.04))
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__({'sankey': {'sankey': axes_obj}},
                         Vector(0, 0), Vector(total_width, total_height),
                         background=False, **kwargs)


class SubfigureD(Subfigure):
    subfigure_label = 'd'
    subfigure_title = 'acetyl_coa_sankey'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        flows = _acetyl_coa_flows()
        print('Panel 5D: acetyl-CoA Sankey (acetyl-CoA 当量, plan §4 系数)')
        element = _place(AcetylCoaSankeyFigure,
                         {'ctrl': flows[Keyword.ad], 'sema': flows[Keyword.sema]},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel E: 关键 flux 柱+单点+p (TODO: CrossTissueBarFigure/BarDataFigure, task 54)
# ============================================================

# 按 plan §1 稳健程度排列 (稳→趋势)
KEY_FLUXES_5E = [
    'astrocyte--LIPID_c',     # de novo 脂质合成 2.9→52.2 p=0.026
    'astrocyte--GLNS_c',      # 谷氨酰胺合成 p=0.026
    'neuron--LDH_c',          # 乳酸/糖酵解 p=0.065 趋势
    'liver--PYK_c',           # 肝脏糖酵解 p=0.065 趋势
    'astrocyte--PDH_m',       # glucose 入口对照 (~10 不变)
]


class FluxBarAxes(DataFigureAxes):
    """关键 flux 柱图: 每反应 ctrl vs sema 柱 + 6 只单点 overlay + p 星号。"""
    def __init__(self, stats_df, bottom_left, size, group_colors=None, **kwargs):
        self.stats = stats_df
        self.group_colors = group_colors or GROUP_COLORS
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        df = self.stats
        x = np.arange(len(df))
        w = 0.38
        ctrl = df['ctrl_median'].values
        sema = df['sema_median'].values
        ax.bar(x - w/2, ctrl, w, color=self.group_colors[Keyword.ad], label='AD (ctrl)')
        ax.bar(x + w/2, sema, w, color=self.group_colors[Keyword.sema], label='Sema')
        for i, row in df.iterrows():
            ax.scatter(np.full(len(row['ctrl_values']), i - w/2), row['ctrl_values'],
                       color='black', s=10, zorder=3, alpha=0.7)
            ax.scatter(np.full(len(row['sema_values']), i + w/2), row['sema_values'],
                       color='black', s=10, zorder=3, alpha=0.7)
        ax.set_ylim(0, float(max(ctrl.max(), sema.max())) * 1.30)
        for i, p in enumerate(df['p_mw']):
            star = '**' if p < 0.01 else '*' if p < 0.05 else '.' if p < 0.1 else ''
            if star:
                ax.text(i, max(ctrl[i], sema[i]) * 1.08, star,
                        ha='center', fontsize=8, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([str(r).split('--')[-1] for r in df['reaction']],
                           rotation=25, ha='right', fontsize=7)
        ax.legend(fontsize=6, loc='upper right')
        ax.set_ylabel('Flux', fontsize=7)
        ax.tick_params(labelsize=6)
        return ax, transform


class FluxBarFigure(CompositeFigure):
    """FluxBarAxes 的 CompositeFigure 包装 (供 _place 定位)。"""
    def __init__(self, figure_data_parameter_dict, **kwargs):
        stats_df = default_parameter_extract(figure_data_parameter_dict, 'stats_df', None, pop=True)
        group_colors = default_parameter_extract(figure_data_parameter_dict, 'group_colors', None, pop=True)
        if stats_df is None or stats_df.empty:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return
        total_width = 1.0
        total_height = 0.55
        axes_obj = FluxBarAxes(
            stats_df, Vector(0.10, 0.12), Vector(total_width - 0.12, total_height - 0.18),
            group_colors=group_colors)
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__({'fluxbar': {'fluxbar': axes_obj}},
                         Vector(0, 0), Vector(total_width, total_height), background=False, **kwargs)


class SubfigureE(Subfigure):
    subfigure_label = 'e'
    subfigure_title = 'key_flux_bar'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        stats = _flux_stats(KEY_FLUXES_5E)
        if stats.empty:
            print('Panel 5E: 无 flux 统计数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        print(f'Panel 5E: {len(stats)} 个关键 flux')
        element = _place(FluxBarFigure, {'stats_df': stats, 'group_colors': GROUP_COLORS},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel F: flux vs 认知相关 (虚拟映射占位; TODO: CorrelationScatterGridFigure, task 54)
# ============================================================

class SubfigureF(Subfigure):
    """flux 通量 vs flux 行为学认知 (同拨 flux 小鼠, 真实配对)。
    LIPID_c (flux 通量) vs NOR DI / MWM learn_rate (flux 行为学)。"""
    subfigure_label = 'f'
    subfigure_title = 'flux_cognition_corr'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        from scripts.data.data_loader import common_data_loader, DataSource
        from scripts.data.metabolomics_flux_behaviroal_data_loader.behavioral_data_loader import _flux_sheet_to_mouse_map
        sheet_to_mouse = _flux_sheet_to_mouse_map()
        fbds, _ = common_data_loader(DataSource.flux_behavioral)
        nor = fbds.nor.set_index('ID') if fbds.nor is not None else None
        mwm = fbds.mwm
        med = _flux_per_sample_median()
        flux_key = 'astrocyte--LIPID_c'
        sources = []
        if nor is not None and 'DI_鼻尖' in nor.columns:
            sources.append((nor, 'DI_鼻尖', 'NOR DI (nose)'))
        if mwm is not None and 'learn_rate' in mwm.columns:
            sources.append((mwm, 'learn_rate', 'MWM learn-rate'))
        panels = []
        for src_df, col, label in sources:
            xs, ys, colors = [], [], []
            for s in _ALL_SHEETS:
                mouse = sheet_to_mouse.get(s)
                fv = med[s].get(flux_key)
                if mouse is None or fv is None or pd.isna(fv) or mouse not in src_df.index:
                    continue
                pv = src_df.loc[mouse, col]
                if isinstance(pv, pd.Series):
                    pv = pv.iloc[0]
                if pd.isna(pv):
                    continue
                xs.append(float(fv))
                ys.append(float(pv))
                colors.append(GROUP_COLORS[_SAMPLE_GROUP[s]])
            if len(xs) < 3:
                continue
            x_arr, y_arr = np.array(xs), np.array(ys)
            z = np.polyfit(x_arr, y_arr, 1)
            fit_x = np.linspace(x_arr.min(), x_arr.max(), 20)
            r = float(np.corrcoef(x_arr, y_arr)[0, 1]) if len(x_arr) > 2 else 0
            panels.append({
                'x': x_arr, 'y': y_arr, 'group_colors': colors,
                'fit_x': fit_x, 'fit_y': np.polyval(z, fit_x),
                'title': f'LIPID_c vs {label}\nr={r:+.2f}',
                'x_label': 'Astrocyte LIPID_c flux', 'y_label': label,
            })
        print(f'Panel 5F: flux 通量 vs flux 行为学 (真实配对, {len(panels)} panel)')
        if not panels:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        element = _place(Elements.CorrelationScatterGridFigure, {'panels': panels},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Figure 5
# ============================================================

class Figure(GeneralFigure):
    figure_label = 'figure_5'
    figure_title = 'Figure 5'

    def __init__(self):
        subfigure_class_list = [
            SubfigureA, SubfigureB, SubfigureC, SubfigureD, SubfigureE, SubfigureF,
        ]
        figure_layout_list = [
            (0.22, [(0.5, 'a'), (0.5, 'b')]),  # A 示意图 + B MID
            (0.26, [(1.0, 'c')]),              # C flux 网络 (占位)
            (0.26, [(0.5, 'd'), (0.5, 'e')]),  # D Sankey + E 柱
            (0.26, [(1.0, 'f')]),              # F 相关
        ]
        super().__init__(
            self.figure_label, subfigure_class_list, figure_layout_list, {},
            figure_title=self.figure_title)
