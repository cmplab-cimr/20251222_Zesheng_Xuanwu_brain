# -*- coding: utf-8 -*-
"""表型-代谢物相关分析

行为学指标（NOR DI、旷场、Y迷宫）vs 关键脂质（脑/血浆 S-plot 候选代谢物）
Pearson / Spearman 相关 + BH-FDR 校正 + 相关热图 + 散点图矩阵
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

TISSUE_NAMES_ZH = {
    'cortex': '皮层', 'hippocampus': '海马', 'hypothalamus': '下丘脑',
    'liver': '肝', 'kidney': '肾', 'fat': '脂肪',
    'gastrocnemius': '腓肠肌', 'quadriceps': '股四头肌', 'ileum': '回肠',
    'serum': '血清',
}


# ============================================================
# 表型指标提取
# ============================================================

def extract_phenotype_metrics(behavioral_ds) -> pd.DataFrame:
    """提取核心行为学指标 (index = 动物 ID 大写), 含合成认知评分 CCS。

    认知成分 (越高 = 认知越好): MWM_learn_rate, NOR_DI_鼻尖
    CCS = zscore(两成分, nan_policy='omit').mean(axis=1);
    NOR_DI_中心 用药后方向相反 (剔出 CCS, 保留列); OF_区域时间 作协变量。
    """
    nor = behavioral_ds.nor.set_index('ID')
    of = behavioral_ds.open_field.set_index('ID')

    mwm_lr = (behavioral_ds.mwm['learn_rate']
              if behavioral_ds.mwm is not None else pd.Series(np.nan, dtype=float))
    metrics = pd.DataFrame({
        'MWM_learn_rate': mwm_lr,
        'NOR_DI_鼻尖': nor['DI_鼻尖'],
        'NOR_DI_中心': nor['DI_中心'],
        'OF_区域时间': of['区域时间(s)'],
    })
    metrics.index = metrics.index.str.upper()

    # 合成认知综合评分 (cognitive composite score); 越高 = 认知越好
    # 仅纳入用药后方向正确的认知成分; NOR_DI_中心 用药后反向 (d=+0.25), 剔出 CCS
    cog_cols = ['MWM_learn_rate', 'NOR_DI_鼻尖']
    metrics['CCS'] = metrics[cog_cols].apply(stats.zscore, nan_policy='omit').mean(axis=1)
    return metrics


# ============================================================
# 关键代谢物提取
# ============================================================

def get_key_metabolites(ds, tissue: str, group1: str, group2: str,
                        p_corr_threshold: float = 0.5,
                        vip_threshold: float = 1.0,
                        max_metabolites: int = 30) -> list[str]:
    """从 S-plot 候选中提取关键代谢物 Index 列表"""
    from .splot import run_splot_oplsda, identify_candidates

    result = run_splot_oplsda(ds, tissue, group1, group2)
    if result is None:
        return []

    cands = identify_candidates(result, p_corr_threshold, vip_threshold)
    sig = cands[cands['candidate']].index.tolist()
    return sig[:max_metabolites]


def get_metabolite_matrix(ds, tissue: str, metabolite_ids: list[str]) -> pd.DataFrame:
    """获取指定组织的代谢物矩阵（代谢物 × 动物ID）"""
    data = ds.get_data_by_tissue(tissue)
    common = [m for m in metabolite_ids if m in data.index]
    if not common:
        return pd.DataFrame()
    sub = data.loc[common]
    # 列名 C11-hippocampus → C11
    sub.columns = [c.split('-')[0].upper() for c in sub.columns]
    return sub


# ============================================================
# 相关分析
# ============================================================

def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]
    adjusted = np.empty(n)
    adjusted[-1] = sorted_p[-1]
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i + 1], sorted_p[i] * n / (i + 1))
    result = np.empty(n)
    result[sorted_idx] = adjusted
    return np.clip(result, 0, 1)


def compute_correlations(
        metabolite_matrix: pd.DataFrame,
        phenotype_matrix: pd.DataFrame,
        method: str = 'spearman',
) -> pd.DataFrame:
    """计算代谢物 × 表型 相关矩阵

    Returns: DataFrame with columns [metabolite, phenotype, r, p_value, adj_p, metabolite_name]
    """
    # 取共有动物
    common_animals = sorted(set(metabolite_matrix.columns) & set(phenotype_matrix.index))
    if len(common_animals) < 4:
        return pd.DataFrame()

    met = metabolite_matrix[common_animals]
    pheno = phenotype_matrix.loc[common_animals]

    records = []
    p_vals = []
    for met_id in met.index:
        x = pd.to_numeric(met.loc[met_id], errors='coerce').values
        for pheno_col in pheno.columns:
            y = pd.to_numeric(pheno[pheno_col], errors='coerce').values
            # mask NaN
            valid = ~(np.isnan(x) | np.isnan(y))
            if valid.sum() < 4:
                r, p = np.nan, np.nan
            else:
                if method == 'pearson':
                    r, p = stats.pearsonr(x[valid], y[valid])
                else:
                    r, p = stats.spearmanr(x[valid], y[valid])
            records.append({
                'metabolite_id': met_id,
                'phenotype': pheno_col,
                'r': r, 'p_value': p,
            })
            p_vals.append(p if not np.isnan(p) else 1.0)

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    adj_p = _benjamini_hochberg(np.array(p_vals))
    df['adj_p'] = adj_p
    df['significant'] = df['adj_p'] < 0.1
    return df


def run_phenotype_correlation(
        ds, behavioral_ds, tissues: list[str],
        group1: str, group2: str,
        method: str = 'spearman',
        p_corr_threshold: float = 0.5,
        vip_threshold: float = 1.0,
        max_metabolites: int = 30,
        extra_metabolite_ids: list[str] | None = None,
) -> dict[str, pd.DataFrame]:
    """多组织表型-代谢物相关分析

    extra_metabolite_ids: 额外的代谢物 Index 列表，会与 S-plot 候选取并集。
                          用于确保特定代谢物（如脂质类）被纳入分析。
                          默认 None 时，自动纳入所有 Class I 脂质类代谢物。

    Returns: {tissue: correlation DataFrame}
    """
    pheno = extract_phenotype_metrics(behavioral_ds)
    results = {}

    # 自动收集所有 Class I 脂质类代谢物 Index（确保脂质不被遗漏）
    if extra_metabolite_ids is None:
        from scripts.figures.common.common_figure_materials import is_lipid_class
        info = ds.metabolite_info
        if 'Class I' in info.columns:
            extra_metabolite_ids = [
                str(i) for i, c in info['Class I'].items()
                if is_lipid_class(c)
            ]

    for tissue in tissues:
        key_mets = get_key_metabolites(ds, tissue, group1, group2,
                                       p_corr_threshold, vip_threshold, max_metabolites)
        # 合并额外代谢物
        if extra_metabolite_ids:
            extra_in_tissue = [m for m in extra_metabolite_ids
                              if m in ds.get_data_by_tissue(tissue).index]
            key_mets = list(dict.fromkeys(key_mets + extra_in_tissue))  # 保序去重

        if not key_mets:
            print(f'{tissue}: 无关键代谢物')
            continue

        met_matrix = get_metabolite_matrix(ds, tissue, key_mets)
        if met_matrix.empty:
            print(f'{tissue}: 代谢物矩阵为空')
            continue

        corr = compute_correlations(met_matrix, pheno, method)
        if corr.empty:
            print(f'{tissue}: 相关分析为空')
            continue

        # 附加代谢物名称
        info = ds.metabolite_info
        if 'Compounds' in info.columns:
            name_map = info['Compounds'].to_dict()
            corr['metabolite_name'] = corr['metabolite_id'].map(name_map)

        results[tissue] = corr
        n_sig = corr['significant'].sum()
        print(f'{tissue}: {len(key_mets)} 关键代谢物 × {pheno.shape[1]} 表型指标, '
              f'{n_sig} 显著相关 (adj_p<0.1)')

    return results


# ============================================================
# 可视化
# ============================================================

def _short_name(name: str, max_len: int = 22) -> str:
    if pd.isna(name) or not name:
        return ''
    name = str(name)
    return name[:max_len] + '..' if len(name) > max_len else name


def plot_correlation_heatmap(
        corr_df: pd.DataFrame,
        output_path: str | Path,
        p_threshold: float = 0.1,
        title: str = '',
):
    """相关热图：代谢物 (y) × 表型 (x)，颜色 = r，标记显著性"""
    if corr_df.empty:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 筛选至少在一个表型上显著的代谢物
    sig_mets = corr_df[corr_df['significant']]['metabolite_id'].unique()
    if len(sig_mets) == 0:
        sig_mets = corr_df.nsmallest(15, 'adj_p')['metabolite_id'].unique()
    plot_df = corr_df[corr_df['metabolite_id'].isin(sig_mets)].copy()

    # 构建矩阵
    pheno_order = ['MWM_learn_rate', 'NOR_DI_鼻尖', 'NOR_DI_中心', 'CCS', 'OF_区域时间']
    pheno_order = [p for p in pheno_order if p in plot_df['phenotype'].unique()]

    # 代谢物按平均 |r| 降序
    met_avg_r = plot_df.groupby('metabolite_id')['r'].apply(lambda x: x.abs().mean())
    met_order = met_avg_r.sort_values(ascending=False).index.tolist()

    r_matrix = plot_df.pivot(index='metabolite_id', columns='phenotype', values='r')
    p_matrix = plot_df.pivot(index='metabolite_id', columns='phenotype', values='adj_p')
    r_matrix = r_matrix.reindex(index=met_order, columns=pheno_order)
    p_matrix = p_matrix.reindex(index=met_order, columns=pheno_order)

    n_met, n_pheno = r_matrix.shape
    fig_h = max(3 + n_met * 0.35, 4)
    fig_w = max(3 + n_pheno * 1.5, 5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='white')

    vmax = 1.0
    im = ax.imshow(r_matrix.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

    # 显著性标记
    for i in range(n_met):
        for j in range(n_pheno):
            p_val = p_matrix.iloc[i, j]
            r_val = r_matrix.iloc[i, j]
            if pd.isna(r_val):
                continue
            marker = ''
            if pd.notna(p_val):
                if p_val < 0.01:
                    marker = '**'
                elif p_val < 0.05:
                    marker = '*'
                elif p_val < 0.1:
                    marker = '.'
            if marker:
                color = 'white' if abs(r_val) > 0.6 else 'black'
                ax.text(j, i, marker, ha='center', va='center', fontsize=7,
                        color=color, fontweight='bold')

    # 标签
    ax.set_xticks(range(n_pheno))
    ax.set_yticks(range(n_met))

    # 代谢物名称
    info_map = {}
    name_col = 'metabolite_name' if 'metabolite_name' in plot_df.columns else None
    if name_col:
        for _, row in plot_df.drop_duplicates('metabolite_id').iterrows():
            info_map[row['metabolite_id']] = str(row.get('metabolite_name', row['metabolite_id']))
    y_labels = [_short_name(info_map.get(m, str(m))) for m in met_order]
    ax.set_yticklabels(y_labels, fontsize=7, fontfamily='DejaVu Sans')

    x_labels = pheno_order
    ax.set_xticklabels(x_labels, fontsize=8, rotation=30, ha='right')
    ax.xaxis.set_ticks_position('top')

    cb = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label('Correlation (r)', fontsize=9)
    cb.ax.tick_params(labelsize=7)

    if title:
        ax.set_title(title, fontsize=11, fontweight='bold', pad=20)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


def plot_scatter_matrix(
        corr_df: pd.DataFrame,
        metabolite_matrix: pd.DataFrame,
        phenotype_matrix: pd.DataFrame,
        output_path: str | Path,
        top_n: int = 8,
        tissue: str = '',
):
    """Top 相关散点图矩阵"""
    if corr_df.empty:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 取 |r| 最大的 top_n 对
    top = corr_df.reindex(corr_df['r'].abs().sort_values(ascending=False).index).head(top_n)

    n_plots = len(top)
    ncols = min(4, n_plots)
    nrows = (n_plots + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.2 * nrows),
                              facecolor='white')
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = np.atleast_2d(axes)

    common_animals = sorted(set(metabolite_matrix.columns) & set(phenotype_matrix.index))

    for idx, (_, row) in enumerate(top.iterrows()):
        ax = axes.flat[idx]
        met_id = row['metabolite_id']
        pheno_col = row['phenotype']

        x = pd.to_numeric(metabolite_matrix.loc[met_id, common_animals], errors='coerce').values
        y = pd.to_numeric(phenotype_matrix.loc[common_animals, pheno_col], errors='coerce').values
        valid = ~(np.isnan(x) | np.isnan(y))

        if valid.sum() < 3:
            ax.set_visible(False)
            continue

        # 分组着色
        from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
        group_map = {
            'C11': Keyword.ad, 'C12': Keyword.ad, 'C13': Keyword.ad,
            'C14': Keyword.ad, 'C15': Keyword.ad, 'C18': Keyword.ad,
            'C42': Keyword.sema, 'C43': Keyword.sema, 'C44': Keyword.sema,
            'C45': Keyword.sema, 'C46': Keyword.sema, 'C47': Keyword.sema,
        }
        colors = ['#E64B35' if group_map.get(a) == Keyword.ad else '#4DBBD5'
                  for a in common_animals]

        ax.scatter(x[valid], y[valid], c=[c for c, v in zip(colors, valid) if v],
                   s=40, edgecolors='#555', linewidths=0.5, alpha=0.85, zorder=3)

        # 回归线
        if valid.sum() >= 3:
            z = np.polyfit(x[valid], y[valid], 1)
            x_line = np.linspace(x[valid].min(), x[valid].max(), 50)
            ax.plot(x_line, np.polyval(z, x_line), 'k-', lw=1, alpha=0.5)

        met_name = _short_name(str(row.get('metabolite_name', met_id)))
        r_val = row['r']
        p_val = row.get('adj_p', row.get('p_value', np.nan))
        stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
        ax.set_title(f'{met_name}\nr={r_val:.2f}{stars}', fontsize=7)
        ax.set_xlabel(met_name, fontsize=6, fontfamily='DejaVu Sans')
        ax.set_ylabel(pheno_col, fontsize=7)
        ax.tick_params(labelsize=6)

    for idx in range(n_plots, nrows * ncols):
        axes.flat[idx].set_visible(False)

    tissue_zh = TISSUE_NAMES_ZH.get(tissue, tissue)
    fig.suptitle(f'{tissue_zh} — Top Metabolite-Phenotype Correlations',
                 fontsize=11, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


# ============================================================
# 导出
# ============================================================

def export_correlation(results: dict[str, pd.DataFrame], output_path: str | Path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for tissue, df in results.items():
            df.sort_values('adj_p').to_excel(writer, sheet_name=tissue[:31], index=False)
    print(f'已导出: {output_path}')
