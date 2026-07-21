# -*- coding: utf-8 -*-
"""S-plot 分析：基于 OPLS-DA 的候选驱动代谢物筛选

横轴 = p[1]（预测载荷，协方差贡献），纵轴 = p(corr)[1]（相关系数载荷，可靠性）
S 形两端「高贡献 + 高相关」的代谢物为候选驱动/标志物
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

TISSUE_NAMES_ZH = {
    'cortex': '皮层', 'hippocampus': '海马', 'hypothalamus': '下丘脑',
    'liver': '肝', 'kidney': '肾', 'fat': '脂肪',
    'gastrocnemius': '腓肠肌', 'quadriceps': '股四头肌', 'ileum': '回肠',
    'serum': '血清',
}


def run_splot_oplsda(ds, tissue: str, group1: str, group2: str, method: str = 'ttest'):
    """运行 OPLS-DA 并返回 S-plot 数据"""
    from .plsda_analysis import run_oplsda

    g1_data = ds.get_data_by_group(group1, tissue=tissue)
    g2_data = ds.get_data_by_group(group2, tissue=tissue)
    if g1_data.empty or g2_data.empty:
        return None

    common_idx = g1_data.index.intersection(g2_data.index)
    merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)
    g1_cols = list(g1_data.loc[common_idx].columns)
    g2_cols = list(g2_data.loc[common_idx].columns)

    result = run_oplsda(merged, g1_cols, g2_cols,
                        group1_name=group1, group2_name=group2)
    result.tissue = tissue
    result.metabolite_info = ds.metabolite_info
    return result


def identify_candidates(
        result, p_corr_threshold: float = 0.5, vip_threshold: float = 1.0,
) -> pd.DataFrame:
    """筛选候选驱动代谢物: |p(corr)[1]| > threshold 且 VIP > threshold"""
    vip = result.vip_scores
    pc = result.p_corr
    pl = result.p_loadings

    if vip is None or pc is None:
        return pd.DataFrame()

    common = vip.index.intersection(pc.index).intersection(pl.index)
    df = pd.DataFrame({
        'VIP': vip.loc[common],
        'p[1]': pl.loc[common],
        'p(corr)[1]': pc.loc[common],
    })

    # 候选标志物
    df['candidate'] = (df['p(corr)[1]'].abs() >= p_corr_threshold) & (df['VIP'] >= vip_threshold)
    df['direction'] = np.where(df['p(corr)[1]'] >= 0, 'up', 'down')

    # 附加代谢物信息
    if result.metabolite_info is not None:
        info = result.metabolite_info
        info_cols = [c for c in ['Compounds', '物质', 'Class I', 'Class II', 'kegg_map']
                     if c in info.columns]
        df = df.join(info[info_cols], how='left')

    return df.sort_values('p(corr)[1]', key=abs, ascending=False)


def plot_splot(
        result, candidates: pd.DataFrame | None = None,
        output_path: str | Path | None = None,
        p_corr_threshold: float = 0.5,
        vip_threshold: float = 1.0,
        max_labels: int = 15,
):
    """绘制 S-plot (p[1] vs p(corr)[1])"""
    pc = result.p_corr
    pl = result.p_loadings
    vip = result.vip_scores

    if pc is None or pl is None:
        print('无 S-plot 数据')
        return

    common = pc.index.intersection(pl.index)
    x = pl.loc[common].values
    y = pc.loc[common].values
    vip_vals = vip.loc[common].values if vip is not None else np.ones(len(common))

    fig, ax = plt.subplots(figsize=(8, 7), facecolor='white')

    # 背景点
    mask_bg = (np.abs(y) < p_corr_threshold) | (vip_vals < vip_threshold)
    mask_candidate = ~mask_bg

    ax.scatter(x[mask_bg], y[mask_bg], s=12, c='#cccccc', alpha=0.5, edgecolors='none', zorder=2)

    # 候选点（按方向着色）
    up_mask = mask_candidate & (y >= 0)
    down_mask = mask_candidate & (y < 0)
    ax.scatter(x[up_mask], y[up_mask], s=30, c='#E64B35', alpha=0.8, edgecolors='none',
               label=f'Up in {result.group_names[1]}', zorder=3)
    ax.scatter(x[down_mask], y[down_mask], s=30, c='#4DBBD5', alpha=0.8, edgecolors='none',
               label=f'Up in {result.group_names[0]}', zorder=3)

    # 阈值线
    ax.axhline(p_corr_threshold, color='#888', ls='--', lw=0.8, alpha=0.6)
    ax.axhline(-p_corr_threshold, color='#888', ls='--', lw=0.8, alpha=0.6)
    ax.axhline(0, color='#aaa', ls='-', lw=0.5, alpha=0.4)
    ax.axvline(0, color='#aaa', ls='-', lw=0.5, alpha=0.4)

    # 标注 top 候选代谢物
    if candidates is not None and not candidates.empty:
        cands = candidates[candidates['candidate']].head(max_labels)
        info = result.metabolite_info
        for idx, row in cands.iterrows():
            label = str(info.loc[idx, 'Compounds']) if (info is not None and idx in info.index and 'Compounds' in info.columns) else str(idx)
            if len(label) > 20:
                label = label[:18] + '..'
            ax.annotate(label, (row['p[1]'], row['p(corr)[1]']),
                        fontsize=5, alpha=0.85,
                        xytext=(4, 4), textcoords='offset points')

    ax.set_xlabel('p[1] (Covariance with predictive component)', fontsize=10)
    ax.set_ylabel('p(corr)[1] (Reliability / Correlation)', fontsize=10)
    tissue_zh = TISSUE_NAMES_ZH.get(result.tissue, result.tissue)
    ax.set_title(f'{tissue_zh} S-Plot — {result.group_names[0]} vs {result.group_names[1]}\n'
                 f'R2Y={result.r2y:.3f}, Q2={result.q2:.3f}',
                 fontsize=11, fontweight='bold')
    ax.legend(loc='lower right', fontsize=8, frameon=True)

    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f'已保存: {output_path}')
    plt.close(fig)


def plot_vip_bar(
        candidates: pd.DataFrame, output_path: str | Path | None = None,
        max_bars: int = 25, tissue: str = '',
):
    """绘制候选代谢物 VIP 评分条形图"""
    if candidates.empty:
        return

    cands = candidates[candidates['candidate']].head(max_bars).copy()
    if cands.empty:
        return

    # 按方向着色
    colors = ['#E64B35' if d == 'up' else '#4DBBD5' for d in cands['direction']]
    labels = cands.get('Compounds', cands.index.astype(str)).values
    # 截断长名称
    labels = [l[:25] + '..' if len(str(l)) > 25 else str(l) for l in labels]

    fig, ax = plt.subplots(figsize=(6, max(3, len(cands) * 0.3)), facecolor='white')
    y_pos = range(len(cands) - 1, -1, -1)

    ax.barh(y_pos, cands['VIP'].values, color=colors, edgecolor='#555', linewidth=0.5, height=0.7)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(1.0, color='#888', ls='--', lw=0.8, alpha=0.6)
    ax.set_xlabel('VIP Score', fontsize=10)

    tissue_zh = TISSUE_NAMES_ZH.get(tissue, tissue)
    ax.set_title(f'{tissue_zh} — Top Candidate Metabolites (VIP)', fontsize=11, fontweight='bold')

    from matplotlib.patches import Patch
    legend_patches = [Patch(facecolor='#E64B35', label='Up in Sema'),
                      Patch(facecolor='#4DBBD5', label='Up in AD')]
    ax.legend(handles=legend_patches, loc='lower right', fontsize=7, frameon=True)

    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f'已保存: {output_path}')
    plt.close(fig)


def export_splot(candidates: pd.DataFrame, output_path: str | Path):
    """导出 S-plot 候选代谢物到 xlsx"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_excel(output_path, engine='openpyxl')
    print(f'已导出: {output_path}')
