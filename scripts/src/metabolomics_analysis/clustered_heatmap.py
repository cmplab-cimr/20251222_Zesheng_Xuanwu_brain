# -*- coding: utf-8 -*-
"""聚类热图：差异代谢物筛选 + 双向层次聚类热图

筛选标准：VIP > 1 且 p < 0.05 且 |FC| > 阈值（可选 raw p 或 BH-FDR）
可视化：z-score 标准化后双向层次聚类，顶部分组注释条
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance

# 中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 组织中文名
TISSUE_NAMES_ZH = {
    'cortex': '皮层', 'hippocampus': '海马', 'hypothalamus': '下丘脑',
    'liver': '肝', 'kidney': '肾', 'fat': '脂肪',
    'gastrocnemius': '腓肠肌', 'quadriceps': '股四头肌', 'ileum': '回肠',
    'serum': '血清',
}


# ============================================================
# 数据处理
# ============================================================

def _get_vip_scores(ds, tissue, group1, group2):
    """运行 OPLS-DA 并返回 VIP scores (pd.Series)"""
    from .plsda_analysis import run_oplsda

    g1_data = ds.get_data_by_group(group1, tissue=tissue)
    g2_data = ds.get_data_by_group(group2, tissue=tissue)
    common_idx = g1_data.index.intersection(g2_data.index)
    merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)
    g1_cols = list(g1_data.loc[common_idx].columns)
    g2_cols = list(g2_data.loc[common_idx].columns)

    result = run_oplsda(merged, g1_cols, g2_cols,
                        group1_name=group1, group2_name=group2)
    return result.vip_scores


def _get_diff_results(ds, tissue, group1, group2, method):
    """运行差异分析并返回结果 DataFrame"""
    from .pathway_enrichment import identify_differential_metabolites

    g1_data = ds.get_data_by_group(group1, tissue=tissue)
    g2_data = ds.get_data_by_group(group2, tissue=tissue)
    common_idx = g1_data.index.intersection(g2_data.index)
    g1_data, g2_data = g1_data.loc[common_idx], g2_data.loc[common_idx]

    merged = pd.concat([g1_data, g2_data], axis=1)
    return identify_differential_metabolites(
        merged, list(g1_data.columns), list(g2_data.columns),
        method=method, p_threshold=1.0, fc_threshold=1.0,
    )


def filter_metabolites(
        vip_scores: pd.Series,
        diff_results: pd.DataFrame,
        vip_threshold: float = 1.0,
        p_threshold: float = 0.05,
        fc_threshold: float = 1.5,
        use_fdr: bool = False,
) -> pd.Index:
    """三重筛选：VIP > threshold & p < p_threshold & |FC| > fc_threshold

    Args:
        use_fdr: True 用 adj_p_value (BH-FDR)，False 用 raw p_value
    """
    common = vip_scores.index.intersection(diff_results.index)
    vip = vip_scores.loc[common]
    diff = diff_results.loc[common]

    p_col = 'adj_p_value' if use_fdr else 'p_value'
    mask = (
        (vip > vip_threshold)
        & (diff[p_col] < p_threshold)
        & (diff[p_col].notna())
        & (np.abs(diff['log2FC']) >= np.log2(fc_threshold))
    )
    return common[mask]


def _zscore_rows(data: np.ndarray) -> np.ndarray:
    """沿 axis=1 做 z-score（行标准化），处理 NaN"""
    with np.errstate(divide='ignore', invalid='ignore'):
        mu = np.nanmean(data, axis=1, keepdims=True)
        sigma = np.nanstd(data, axis=1, keepdims=True)
        sigma[sigma == 0] = 1.0
        z = (data - mu) / sigma
    return z


# ============================================================
# 聚类
# ============================================================

def _cluster_linkage(data: np.ndarray, axis: int = 0) -> np.ndarray:
    """Ward 层次聚类返回 linkage matrix"""
    if axis == 1:
        data = data.T
    # NaN → 0 用于距离计算
    d = np.nan_to_num(data, nan=0.0)
    dist = distance.pdist(d, metric='euclidean')
    return hierarchy.linkage(dist, method='ward')


def _leaves_order(Z: np.ndarray) -> list[int]:
    """从 linkage matrix 获取叶节点顺序"""
    return hierarchy.leaves_list(Z)


# ============================================================
# 绘图
# ============================================================

# 分组配色
_GROUP_COLORS = {
    'ad': '#E64B35',
    'semaglutide': '#4DBBD5',
}


def plot_clustered_heatmap(
        z_matrix: np.ndarray,
        row_order: list[int],
        col_order: list[int],
        row_linkage: np.ndarray,
        col_linkage: np.ndarray,
        metabolite_labels: list[str],
        sample_labels: list[str],
        groups: list[str],
        output_path: str | Path,
        tissue: str = '',
        figsize: tuple[float, float] | None = None,
):
    """绘制聚类热图（Cell Metab 风格）

    布局：
      右上: 列注释条
      中上: 列树状图
      左侧: 行树状图
      中心: 热图
      右下: colorbar
    """
    output_path = Path(output_path)

    n_rows, n_cols = z_matrix.shape

    # 按聚类顺序重排
    z_ordered = z_matrix[np.ix_(row_order, col_order)]
    groups_ordered = [groups[i] for i in col_order]
    sample_labels_ordered = [sample_labels[i] for i in col_order]
    metabolite_labels_ordered = [metabolite_labels[i] for i in row_order]

    # 自动 figsize
    if figsize is None:
        w = max(3.5 + n_cols * 0.45, 8)
        h = max(2.0 + n_rows * 0.22, 6)
        figsize = (min(w, 16), min(h, 18))

    # 布局比例
    left_dend_frac = 0.12
    top_dend_frac = 0.10
    anno_frac = 0.03
    cbar_frac = 0.04

    fig = plt.figure(figsize=figsize, facecolor='white')

    # gridspec: [行树状图 | 热图+注释 | colorbar]
    #           [行标签    | 热图     | colorbar]
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    gs_main = GridSpec(1, 3, figure=fig,
                       width_ratios=[left_dend_frac, 1 - left_dend_frac - cbar_frac - 0.02, cbar_frac],
                       left=0.01, right=0.99, bottom=0.02, top=0.95,
                       wspace=0.02)

    # 左侧 dendrogram
    gs_left = GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_main[0])
    ax_row_dend = fig.add_subplot(gs_left[0, 0])

    # 中间区域：注释条 + dendrogram + heatmap
    gs_mid = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_main[1],
                                     height_ratios=[anno_frac, top_dend_frac, 1 - anno_frac - top_dend_frac],
                                     hspace=0.0)
    ax_anno = fig.add_subplot(gs_mid[0, 0])
    ax_col_dend = fig.add_subplot(gs_mid[1, 0])
    ax_hm = fig.add_subplot(gs_mid[2, 0])

    # 右侧 colorbar
    gs_right = GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_main[2])
    ax_cbar = fig.add_subplot(gs_right[0, 0])

    # ---- 绘制行 dendrogram ----
    hierarchy.dendrogram(row_linkage, orientation='left', no_labels=True,
                         ax=ax_row_dend, color_threshold=0,
                         link_color_func=lambda *a: '#333333')
    ax_row_dend.invert_yaxis()
    ax_row_dend.set_xticks([])
    ax_row_dend.set_yticks([])
    ax_row_dend.spines[:].set_visible(False)

    # ---- 绘制列 dendrogram ----
    hierarchy.dendrogram(col_linkage, orientation='top', no_labels=True,
                         ax=ax_col_dend, color_threshold=0,
                         link_color_func=lambda *a: '#333333',
                         above_threshold_color='#333333')
    ax_col_dend.set_xticks([])
    ax_col_dend.set_yticks([])
    ax_col_dend.spines[:].set_visible(False)
    ax_col_dend.set_xlim(-0.5, n_cols - 0.5)

    # ---- 绘制分组注释条 ----
    _color_to_val = {c: i for i, c in enumerate(sorted(set(_GROUP_COLORS.values())))}
    anno_vals = np.array([_color_to_val.get(_GROUP_COLORS.get(g, '#999999'), 0)
                          for g in groups_ordered]).reshape(1, -1)
    ax_anno.imshow(anno_vals, aspect='auto', interpolation='nearest',
                    cmap=matplotlib.colors.ListedColormap(
                        ['#999999'] + [c for c in sorted(set(_GROUP_COLORS.values()))]),
                    vmin=-0.5, vmax=len(_GROUP_COLORS) - 0.5)
    ax_anno.set_xticks([])
    ax_anno.set_yticks([])
    ax_anno.spines[:].set_visible(False)
    ax_anno.set_xlim(-0.5, n_cols - 0.5)

    # ---- 绘制热图 ----
    vmax = np.nanmax(np.abs(z_ordered[np.isfinite(z_ordered)])) if np.any(np.isfinite(z_ordered)) else 3
    vmax = min(vmax, 3)
    ax_hm.imshow(z_ordered, aspect='auto', cmap='RdBu_r',
                  vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax_hm.set_yticks(range(n_rows))
    ax_hm.set_yticklabels(metabolite_labels_ordered, fontsize=6, fontfamily='DejaVu Sans')
    ax_hm.set_xticks(range(n_cols))
    ax_hm.set_xticklabels(sample_labels_ordered, fontsize=7, rotation=45,
                           ha='right', fontfamily='DejaVu Sans')
    ax_hm.spines['top'].set_visible(False)
    ax_hm.spines['right'].set_visible(False)
    ax_hm.set_xlim(-0.5, n_cols - 0.5)
    ax_hm.set_ylim(n_rows - 0.5, -0.5)

    # ---- colorbar ----
    norm = Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, cax=ax_cbar)
    cb.set_label('Z-score', fontsize=9)
    cb.ax.tick_params(labelsize=7)

    # 标题
    tissue_zh = TISSUE_NAMES_ZH.get(tissue, tissue)
    n_met = len(metabolite_labels)
    fig.suptitle(f'{tissue_zh} ({tissue}) — {n_met} differential metabolites',
                 fontsize=12, fontweight='bold', y=0.98)

    # 分组图例
    legend_patches = [Patch(facecolor=_GROUP_COLORS[g], label=g)
                      for g in sorted(set(groups_ordered))]
    fig.legend(handles=legend_patches, loc='lower right',
               fontsize=8, frameon=False, ncol=2,
               bbox_to_anchor=(0.95, 0.01))

    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


def plot_clustered_heatmap_grid(
        results: list[dict],
        output_path: str | Path,
        ncols: int = 2,
):
    """多组织聚类热图网格（海马+皮层并排）"""
    if not results:
        return
    output_path = Path(output_path)
    nrows = (len(results) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(6 * ncols, 5 * nrows),
                              facecolor='white')
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1 or ncols == 1:
        axes = axes.reshape(nrows, ncols)

    for idx, r in enumerate(results):
        ax = axes.flat[idx]
        z_ordered = r['z_ordered']
        groups_ordered = r['groups_ordered']
        sample_labels = r['sample_labels']
        metabolite_labels = r['metabolite_labels']
        tissue = r['tissue']

        vmax = min(np.nanmax(np.abs(z_ordered[np.isfinite(z_ordered)])) if np.any(np.isfinite(z_ordered)) else 3, 3)
        ax.imshow(z_ordered, aspect='auto', cmap='RdBu_r',
                  vmin=-vmax, vmax=vmax, interpolation='nearest')
        ax.set_xticks(range(len(sample_labels)))
        ax.set_xticklabels(sample_labels, fontsize=6, rotation=45, ha='right')
        ax.set_yticks(range(len(metabolite_labels)))
        ax.set_yticklabels(metabolite_labels, fontsize=5, fontfamily='DejaVu Sans')

        # 顶部注释条
        anno_colors = [_GROUP_COLORS.get(g, '#999') for g in groups_ordered]
        for j, c in enumerate(anno_colors):
            ax.add_patch(plt.Rectangle((j - 0.5, -2.5), 1, 2, color=c, clip_on=False))

        tissue_zh = TISSUE_NAMES_ZH.get(tissue, tissue)
        n_met = len(metabolite_labels)
        ax.set_title(f'{tissue_zh} ({n_met} metabolites)', fontsize=10, fontweight='bold')

    # 隐藏多余子图
    for idx in range(len(results), nrows * ncols):
        axes.flat[idx].set_visible(False)

    # 图例
    legend_patches = [Patch(facecolor=_GROUP_COLORS[g], label=g)
                      for g in sorted(_GROUP_COLORS.keys())]
    fig.legend(handles=legend_patches, loc='lower right', fontsize=8,
               frameon=False, ncol=2)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


# ============================================================
# 主入口
# ============================================================

def run_clustered_heatmap(
        ds, tissue: str, group1: str, group2: str,
        method: str = 'ttest',
        p_threshold: float = 0.05,
        fc_threshold: float = 1.5,
        vip_threshold: float = 1.0,
        use_fdr: bool = False,
        diff_result=None, oplsda_result=None,
) -> dict | None:
    """对单个组织运行聚类热图分析

    Args:
        diff_result: 预算的差异分析 DataFrame（来自 _run_differential），提供则跳过内部 _get_diff_results
        oplsda_result: 预算的 OPLSDAResult（来自 _run_oplsda_single），提供则取其 vip_scores 跳过内部 _get_vip_scores

    Returns:
        dict with keys: tissue, filtered_index, z_matrix, row_order, col_order,
                        row_linkage, col_linkage, metabolite_labels, sample_labels,
                        groups, z_ordered, groups_ordered, vip_scores, diff_results
        或 None（若无差异代谢物）
    """
    # VIP（优先用预算 OPLSDAResult，避免重复跑 OPLS-DA）
    if oplsda_result is not None and oplsda_result.vip_scores is not None:
        vip_scores = oplsda_result.vip_scores
    else:
        vip_scores = _get_vip_scores(ds, tissue, group1, group2)
    # 差异分析（优先用预算 diff，避免重复计算）
    diff_results = diff_result if diff_result is not None else _get_diff_results(ds, tissue, group1, group2, method)
    # 三重筛选
    filtered_idx = filter_metabolites(vip_scores, diff_results,
                                      vip_threshold, p_threshold, fc_threshold,
                                      use_fdr=use_fdr)

    p_label = 'adj_p' if use_fdr else 'raw_p'
    n_filt = len(filtered_idx)
    print(f'{tissue}: VIP>{vip_threshold} & {p_label}<{p_threshold} & |FC|>{fc_threshold} → {n_filt} 代谢物')
    if n_filt < 3:
        print(f'  差异代谢物不足 3 个，跳过')
        return None

    # 取全部样本的筛选后数据
    g1_data = ds.get_data_by_group(group1, tissue=tissue)
    g2_data = ds.get_data_by_group(group2, tissue=tissue)
    common_idx = g1_data.index.intersection(g2_data.index).intersection(filtered_idx)
    if len(common_idx) < 3:
        print(f'  共有差异代谢物不足 3 个，跳过')
        return None

    merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)
    # 去掉全 NaN 列
    merged = merged.loc[:, merged.notna().any(axis=0)]

    # Z-score
    z_matrix = _zscore_rows(merged.values)

    # 聚类
    row_linkage = _cluster_linkage(z_matrix, axis=0)
    col_linkage = _cluster_linkage(z_matrix, axis=1)
    row_order = _leaves_order(row_linkage)
    col_order = _leaves_order(col_linkage)

    # 分组标签
    g1_cols = set(g1_data.loc[common_idx].columns)
    groups = [group1 if c in g1_cols else group2 for c in merged.columns]

    # 代谢物名称
    info = ds.metabolite_info
    if 'Compounds' in info.columns:
        met_labels = [str(info.loc[i, 'Compounds']) if i in info.index else str(i)
                      for i in merged.index]
    else:
        met_labels = [str(i) for i in merged.index]

    # 样本名简化（只取动物 ID）
    sample_labels = [c.split('-')[0] if '-' in c else c for c in merged.columns]

    # Class I 分类（dict {metabolite_index: class_label}），供 Panel B 左侧分类色条
    if 'Class I' in info.columns:
        class_lookup = info['Class I'].reindex(merged.index)
        class_i = {i: (str(v) if pd.notna(v) else 'Other') for i, v in class_lookup.items()}
    else:
        class_i = {i: 'Other' for i in merged.index}
    # Compounds 名 → index 反查表（供 Panel B highlight 按名称指定代谢物）
    if 'Compounds' in info.columns:
        comp_lookup = info['Compounds'].reindex(merged.index)
        metabolite_name_to_index = {str(v): i for i, v in comp_lookup.items() if pd.notna(v)}
    else:
        metabolite_name_to_index = {}

    # 按 col_order 重排后的分组
    z_ordered = z_matrix[np.ix_(row_order, col_order)]
    groups_ordered = [groups[i] for i in col_order]
    sample_labels_ordered = [sample_labels[i] for i in col_order]
    met_labels_ordered = [met_labels[i] for i in row_order]

    return {
        'tissue': tissue,
        'filtered_index': common_idx,
        'z_matrix': z_matrix,
        'row_order': row_order,
        'col_order': col_order,
        'row_linkage': row_linkage,
        'col_linkage': col_linkage,
        'metabolite_labels': met_labels,
        'sample_labels': sample_labels,
        'groups': groups,
        'z_ordered': z_ordered,
        'groups_ordered': groups_ordered,
        'sample_labels_ordered': sample_labels_ordered,
        'metabolite_labels_ordered': met_labels_ordered,
        'vip_scores': vip_scores,
        'diff_results': diff_results,
        'class_i': class_i,
        'metabolite_name_to_index': metabolite_name_to_index,
    }


def export_clustered_heatmap(result: dict, output_path: str | Path):
    """导出聚类热图结果到 xlsx"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ds_idx = result['filtered_index']
    vip = result['vip_scores']
    diff = result['diff_results']

    common = ds_idx.intersection(vip.index).intersection(diff.index)
    df = pd.DataFrame({
        'VIP': vip.loc[common],
        'log2FC': diff.loc[common, 'log2FC'],
        'p_value': diff.loc[common, 'p_value'],
        'adj_p_value': diff.loc[common, 'adj_p_value'],
    })
    # 按 VIP 降序
    df = df.sort_values('VIP', ascending=False)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='差异代谢物')
        # 聚类后的 z-score 矩阵
        z_df = pd.DataFrame(
            result['z_matrix'][np.ix_(result['row_order'], result['col_order'])],
            index=[result['metabolite_labels'][i] for i in result['row_order']],
            columns=[result['sample_labels'][i] for i in result['col_order']],
        )
        z_df.to_excel(writer, sheet_name='Z-score矩阵')

    print(f'已导出: {output_path}')
