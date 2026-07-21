# -*- coding: utf-8 -*-
"""多组织 KEGG 通路富集 + MetaboAnalyst 拓扑分析 + 气泡图

分析流程：
1. 各组织分别做 KEGG 通路富集（超几何检验）+ 差异方向
2. 计算通路拓扑影响因子（pathway impact，betweenness centrality）
3. 汇总为 organs × pathways 气泡图：大小=impact, 颜色=-log10(p)
"""

import json
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import networkx as nx

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

TISSUE_NAMES_ZH = {
    'cortex': '皮层', 'hippocampus': '海马', 'hypothalamus': '下丘脑',
    'liver': '肝', 'kidney': '肾', 'fat': '脂肪',
    'gastrocnemius': '腓肠肌', 'quadriceps': '股四头肌', 'ileum': '回肠',
    'serum': '血清',
}

# 拓扑缓存
_TOPO_CACHE = Path(__file__).parent / 'kegg_pathway_topology.json'

# 默认组织顺序（脑区 → 外周 → 血液）
TISSUE_ORDER = [
    'hippocampus', 'cortex', 'hypothalamus',
    'liver', 'kidney', 'fat', 'gastrocnemius', 'quadriceps', 'ileum',
    'serum',
]


# ============================================================
# KEGG API
# ============================================================

def _kegg_get(url: str, timeout: int = 15) -> str:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.read().decode('utf-8')
    except Exception:
        return ''


def _fetch_pathway_reactions(pathway_id: str) -> list[str]:
    """获取通路的所有反应 ID"""
    map_id = pathway_id.replace('ko', 'map')
    text = _kegg_get(f'https://rest.kegg.jp/link/rn/{map_id}')
    reactions = []
    for line in text.strip().split('\n'):
        if '\t' not in line:
            continue
        rn = line.split('\t')[1].strip()
        if rn.startswith('rn:'):
            reactions.append(rn)
    return reactions


def _fetch_reactions_compounds(
        reaction_ids: list[str], batch_size: int = 10,
) -> dict[str, list[str]]:
    """批量获取反应的化合物 {rn_id: [C-number]}"""
    result: dict[str, list[str]] = defaultdict(list)
    for i in range(0, len(reaction_ids), batch_size):
        batch = reaction_ids[i:i + batch_size]
        text = _kegg_get(f'https://rest.kegg.jp/link/cpd/{"+".join(batch)}')
        for line in text.strip().split('\n'):
            if '\t' not in line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rn = parts[0].strip()
            cpd = parts[1].strip().replace('cpd:', '')
            result[rn].append(cpd)
    return dict(result)


# ============================================================
# 拓扑分析
# ============================================================

def _build_graph(rxn_compounds: dict[str, list[str]]) -> nx.Graph:
    """从反应-化合物映射构建通路无向图"""
    G = nx.Graph()
    for compounds in rxn_compounds.values():
        if len(compounds) < 2:
            continue
        G.add_nodes_from(compounds)
        for i in range(len(compounds)):
            for j in range(i + 1, len(compounds)):
                G.add_edge(compounds[i], compounds[j])
    return G


def _compute_relative_betweenness(G: nx.Graph) -> dict[str, float]:
    """计算归一化 betweenness centrality (MetaboAnalyst 用法)"""
    if len(G.nodes) == 0:
        return {}
    bc = nx.betweenness_centrality(G, normalized=True)
    max_bc = max(bc.values()) if bc else 1.0
    if max_bc == 0:
        return {k: 0.0 for k in bc}
    return {k: v / max_bc for k, v in bc.items()}


def _load_topology_cache() -> dict:
    if _TOPO_CACHE.exists():
        return json.loads(_TOPO_CACHE.read_text(encoding='utf-8'))
    return {}


def _save_topology_cache(cache: dict):
    _TOPO_CACHE.write_text(
        json.dumps(cache, ensure_ascii=False, indent=2), encoding='utf-8')


def compute_pathway_topology(pathway_ids: list[str]) -> dict[str, dict[str, float]]:
    """批量计算通路拓扑（带 KEGG 缓存）

    Returns: {pathway_id: {cpd_id: relative_betweenness_centrality}}
    """
    cache = _load_topology_cache()
    results = {}
    new_count = 0

    for pid in pathway_ids:
        if pid in cache:
            results[pid] = cache[pid]
            continue

        map_id = pid.replace('ko', 'map')
        rxns = _fetch_pathway_reactions(pid)
        if not rxns:
            results[pid] = {}
            cache[pid] = {}
            continue

        rc = _fetch_reactions_compounds(rxns)
        G = _build_graph(rc)
        centrality = _compute_relative_betweenness(G)

        results[pid] = centrality
        cache[pid] = centrality
        new_count += 1

    if new_count:
        _save_topology_cache(cache)
        print(f'  新计算 {new_count} 个通路拓扑 (总计缓存 {len(cache)} 个)')

    return results


def pathway_impact(
        pathway_id: str,
        matched_cpd_ids: set[str],
        topology: dict[str, dict[str, float]],
) -> float:
    """pathway impact = Σ(relative BC of matched metabolites)"""
    topo = topology.get(pathway_id, {})
    if not topo:
        return 0.0
    return sum(topo.get(cpd, 0.0) for cpd in matched_cpd_ids)


# ============================================================
# 多组织分析
# ============================================================

def _cpd_id_map(metabolite_info: pd.DataFrame) -> dict[str, str]:
    """Index → KEGG compound ID"""
    if 'cpd_ID' not in metabolite_info.columns:
        return {}
    return {
        str(idx): str(val).strip()
        for idx, val in metabolite_info['cpd_ID'].items()
        if pd.notna(val) and str(val).strip() not in ('', '-')
    }


def run_multi_tissue_pathway(
        ds, tissues: list[str], group1: str, group2: str,
        method: str = 'ttest',
        p_threshold: float = 0.05,
        fc_threshold: float = 1.5,
        diff_results: dict | None = None,
) -> pd.DataFrame:
    """多组织通路富集 + 拓扑 + 方向

    Args:
        diff_results: dict {tissue: diff_df}（来自 _run_differential），提供则跳过内部 identify_differential_metabolites

    Returns: DataFrame[tissue, tissue_zh, pathway_id, pathway_name_zh,
                       p_value, significant_count, pathway_impact, median_log2FC, ...]
    """
    from .pathway_enrichment import (
        identify_differential_metabolites, pathway_enrichment,
        _build_kegg_pathway_map, KEGG_PATHWAY_NAMES,
        KEGG_PATHWAY_NAMES_ZH,
    )

    meta = ds.metabolite_info
    cpd_map = _cpd_id_map(meta)
    pw_map = _build_kegg_pathway_map(meta)

    records = []
    all_pids = set()
    tissue_data = {}

    log2_fc_threshold = np.log2(fc_threshold)
    for tissue in tissues:
        g1 = ds.get_data_by_group(group1, tissue=tissue)
        g2 = ds.get_data_by_group(group2, tissue=tissue)
        if g1.empty or g2.empty:
            print(f'{tissue}: 数据为空，跳过')
            continue

        common = g1.index.intersection(g2.index)
        g1, g2 = g1.loc[common], g2.loc[common]
        merged = pd.concat([g1, g2], axis=1)

        if diff_results is not None and tissue in diff_results:
            diff = diff_results[tissue]
        else:
            diff = identify_differential_metabolites(
                merged, list(g1.columns), list(g2.columns),
                method=method, p_threshold=p_threshold, fc_threshold=fc_threshold)

        # sig_ids 用面板阈值判断（缓存 diff 用宽松阈值，不能直接用其 significant 列）
        sig_mask = (
            (diff['p_value'] < p_threshold)
            & diff['p_value'].notna()
            & (diff['log2FC'].abs() >= log2_fc_threshold)
        )
        sig_ids = set(diff.index[sig_mask])
        all_ids = set(merged.index[merged.notna().any(axis=1)])
        if not sig_ids:
            print(f'{tissue}: 无差异代谢物')
            continue

        enrich = pathway_enrichment(
            sig_ids, all_ids, pw_map, KEGG_PATHWAY_NAMES, metabolite_info=meta)
        if enrich.empty:
            print(f'{tissue}: 无富集通路')
            continue

        tissue_data[tissue] = {'diff': diff, 'enrich': enrich, 'sig_ids': sig_ids}
        all_pids.update(enrich['pathway_id'].tolist())
        n_sig = len(sig_ids)
        print(f'{tissue}: {n_sig} 差异代谢物, {len(enrich)} 富集通路')

    # 批量拓扑
    topology = compute_pathway_topology(list(all_pids)) if all_pids else {}

    # 汇总
    for tissue, data in tissue_data.items():
        diff = data['diff']
        enrich = data['enrich']

        # 差异代谢物的 cpd_ID 集合
        sig_cpd = {cpd_map[i] for i in data['sig_ids'] if i in cpd_map}

        for _, row in enrich.iterrows():
            pid = row['pathway_id']

            # 通路内差异代谢物的 cpd_ID → 计算 impact
            sig_met_str = row.get('significant_metabolites', '')
            met_ids = {s.strip() for s in str(sig_met_str).split(';') if s.strip()} if pd.notna(sig_met_str) else set()
            met_cpd = {cpd_map[m] for m in met_ids if m in cpd_map}
            impact = pathway_impact(pid, met_cpd, topology)

            # 通路方向 = 差异代谢物 median log2FC
            met_in_pw = met_ids & set(diff.index)
            med_fc = diff.loc[list(met_in_pw), 'log2FC'].median() if met_in_pw else 0.0

            records.append({
                'tissue': tissue,
                'tissue_zh': TISSUE_NAMES_ZH.get(tissue, tissue),
                'pathway_id': pid,
                'pathway_name': KEGG_PATHWAY_NAMES.get(pid, pid),
                'pathway_name_zh': KEGG_PATHWAY_NAMES_ZH.get(pid, KEGG_PATHWAY_NAMES.get(pid, pid)),
                'p_value': row['p_value'],
                'significant_count': row['significant_count'],
                'pathway_size': row['pathway_size'],
                'coverage': row['coverage'],
                'pathway_impact': impact,
                'median_log2FC': med_fc,
            })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df['-log10(p)'] = -np.log10(df['p_value'].clip(lower=1e-50))
    return df


# ============================================================
# 气泡图
# ============================================================

def plot_pathway_bubble(
        results: pd.DataFrame,
        output_path: str | Path,
        max_pathways: int = 25,
        p_filter: float = 0.1,
        color_by: str = 'direction',
):
    """绘制 organs × pathways 气泡图

    Args:
        color_by: 'direction'(红蓝=上调下调) 或 'significance'(渐变=-log10p)
    """
    if results.empty:
        print('无数据可绘图')
        return

    output_path = Path(output_path)

    # 筛选显著通路
    sig = results[results['p_value'] < p_filter].copy()
    if sig.empty:
        sig = results.nsmallest(max_pathways, 'p_value').copy()

    # 选出现次数最多的通路（跨组织反复出现 → 系统性证据）
    pw_counts = sig.groupby('pathway_id').size()
    top_pws = pw_counts.nlargest(max_pathways).index.tolist()
    plot_df = sig[sig['pathway_id'].isin(top_pws)].copy()

    # 排序：通路按出现次数降序，组织按固定顺序
    pw_order = (
        plot_df.groupby('pathway_id')['pathway_name_zh'].first()
        .loc[[p for p in top_pws if p in plot_df['pathway_id'].unique()]]
    )
    pw_labels = pw_order.values.tolist()
    pw_ids = pw_order.index.tolist()

    tissue_labels = [TISSUE_NAMES_ZH.get(t, t) for t in TISSUE_ORDER
                     if t in plot_df['tissue'].unique()]
    tissue_ids = [t for t in TISSUE_ORDER if t in plot_df['tissue'].unique()]

    n_pw = len(pw_ids)
    n_tissue = len(tissue_ids)

    # 尺寸
    fig_w = max(3.5 + n_tissue * 1.2, 8)
    fig_h = max(2 + n_pw * 0.42, 5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='white')

    # 映射到数值坐标
    x_map = {t: i for i, t in enumerate(tissue_ids)}
    y_map = {p: i for i, p in enumerate(pw_ids)}

    xs, ys, sizes, colors = [], [], [], []

    for _, row in plot_df.iterrows():
        t = row['tissue']
        p = row['pathway_id']
        if t not in x_map or p not in y_map:
            continue
        xs.append(x_map[t])
        ys.append(y_map[p])

        # size: impact → 点面积
        imp = max(row['pathway_impact'], 0.01)
        sizes.append(imp * 800 + 30)

        if color_by == 'direction':
            colors.append(row['median_log2FC'])
        else:
            colors.append(row['-log10(p)'])

    xs, ys, sizes = np.array(xs), np.array(ys), np.array(sizes)
    colors = np.array(colors)

    if color_by == 'direction':
        vmax = max(abs(colors.min()), abs(colors.max()), 0.5) if len(colors) else 1
        scatter = ax.scatter(xs, ys, s=sizes, c=colors, cmap='RdBu_r',
                             vmin=-vmax, vmax=vmax, edgecolors='#555',
                             linewidths=0.5, alpha=0.85, zorder=3)
        cb = plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.02)
        cb.set_label('Median log2FC (Sema vs AD)', fontsize=8)
    else:
        scatter = ax.scatter(xs, ys, s=sizes, c=colors, cmap='YlOrRd',
                             vmin=0, edgecolors='#555',
                             linewidths=0.5, alpha=0.85, zorder=3)
        cb = plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.02)
        cb.set_label('-log10(p)', fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # 坐标轴
    ax.set_xticks(range(n_tissue))
    ax.set_xticklabels(tissue_labels, fontsize=9, fontweight='bold')
    ax.set_yticks(range(n_pw))
    ax.set_yticklabels(pw_labels, fontsize=8)
    ax.set_xlim(-0.6, n_tissue - 0.4)
    ax.set_ylim(-0.6, n_pw - 0.4)
    ax.invert_yaxis()

    # 网格
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    for i in range(n_pw + 1):
        ax.axhline(i - 0.5, color='#ddd', linewidth=0.5, zorder=1)
    for i in range(n_tissue + 1):
        ax.axvline(i - 0.5, color='#ddd', linewidth=0.5, zorder=1)
    ax.tick_params(axis='both', length=0)

    # impact 图例
    legend_sizes = [0.1, 0.3, 0.5]
    legend_pts = [plt.scatter([], [], s=imp * 800 + 30, c='grey', edgecolors='#555',
                               linewidths=0.5) for imp in legend_sizes]
    legend_labels = [f'{imp:.1f}' for imp in legend_sizes]
    leg = ax.legend(legend_pts, legend_labels, title='Pathway\nImpact',
                    loc='lower right', fontsize=7, title_fontsize=7,
                    frameon=True, framealpha=0.9, borderpad=0.5)

    ax.set_title('KEGG Pathway Enrichment Across Organs', fontsize=11,
                 fontweight='bold', pad=25)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


def plot_metaboanalyst_scatter(
        results: pd.DataFrame,
        output_path: str | Path,
        tissue: str,
):
    """MetaboAnalyst 风格单组织散点图: x=impact, y=-log10(p)"""
    if results.empty:
        return
    output_path = Path(output_path)
    df = results[results['tissue'] == tissue].copy()
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')

    # 阈值线
    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(0.1, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)

    vmax = max(abs(df['median_log2FC'].min()), abs(df['median_log2FC'].max()), 0.5)
    ax.scatter(df['pathway_impact'], df['-log10(p)'],
               s=df['significant_count'] * 30 + 20,
               c=df['median_log2FC'], cmap='RdBu_r',
               vmin=-vmax, vmax=vmax,
               edgecolors='#555', linewidths=0.5, alpha=0.85)

    # 标注 top 通路
    top = df.nsmallest(8, 'p_value')
    for _, row in top.iterrows():
        name = row['pathway_name_zh'] if len(row['pathway_name_zh']) < 25 else row['pathway_name']
        ax.annotate(name, (row['pathway_impact'], row['-log10(p)']),
                    fontsize=6, alpha=0.9,
                    xytext=(5, 5), textcoords='offset points')

    cb = plt.colorbar(cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(-vmax, vmax), cmap='RdBu_r'), ax=ax)
    cb.set_label('Median log2FC', fontsize=8)

    ax.set_xlabel('Pathway Impact', fontsize=10)
    ax.set_ylabel('-log10(p)', fontsize=10)
    tissue_zh = TISSUE_NAMES_ZH.get(tissue, tissue)
    ax.set_title(f'{tissue_zh} — MetaboAnalyst Pathway Analysis', fontsize=11, fontweight='bold')

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'已保存: {output_path}')


# ============================================================
# 导出
# ============================================================

def export_multi_tissue_pathway(results: pd.DataFrame, output_path: str | Path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        results.sort_values(['tissue', 'p_value']).to_excel(writer, sheet_name='多组织通路分析', index=False)
    print(f'已导出: {output_path}')
