# -*- coding: utf-8 -*-
"""代谢物通路富集分析

使用超几何检验进行 KEGG 通路富集分析，整合差异代谢物筛选。

支持的数据库：
- KEGG：基于数据中的 kegg_map 注释（内置）
- Custom：基于用户提供的代谢物-通路映射文件（JSON/CSV）

可扩展数据库：SMPDB, HMDB, Reactome, WikiPathways 等，
只需提供对应的 metabolite -> pathway 映射文件即可。
"""

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# 加载 KEGG 通路中文名映射
_ZH_PATH = Path(__file__).parent / 'kegg_pathway_names_zh.json'
KEGG_PATHWAY_NAMES_ZH = json.loads(_ZH_PATH.read_text(encoding='utf-8')) if _ZH_PATH.exists() else {}


# ---- KEGG 代谢通路 ID -> 名称 ----

KEGG_PATHWAY_NAMES = {
    'ko00010': 'Glycolysis / Gluconeogenesis',
    'ko00020': 'Citrate cycle (TCA cycle)',
    'ko00030': 'Pentose phosphate pathway',
    'ko00040': 'Pentose and glucuronate interconversions',
    'ko00051': 'Fructose and mannose metabolism',
    'ko00052': 'Galactose metabolism',
    'ko00053': 'Ascorbate and aldarate metabolism',
    'ko00061': 'Fatty acid biosynthesis',
    'ko00062': 'Fatty acid elongation',
    'ko00071': 'Fatty acid degradation',
    'ko00072': 'Synthesis and degradation of ketone bodies',
    'ko00100': 'Steroid biosynthesis',
    'ko00120': 'Primary bile acid biosynthesis',
    'ko00140': 'Steroid hormone biosynthesis',
    'ko00190': 'Oxidative phosphorylation',
    'ko00230': 'Purine metabolism',
    'ko00240': 'Pyrimidine metabolism',
    'ko00250': 'Alanine, aspartate and glutamate metabolism',
    'ko00260': 'Glycine, serine and threonine metabolism',
    'ko00270': 'Cysteine and methionine metabolism',
    'ko00280': 'Valine, leucine and isoleucine degradation',
    'ko00290': 'Valine, leucine and isoleucine biosynthesis',
    'ko00300': 'Lysine biosynthesis',
    'ko00310': 'Lysine degradation',
    'ko00330': 'Arginine biosynthesis',
    'ko00340': 'Histidine metabolism',
    'ko00350': 'Tyrosine metabolism',
    'ko00360': 'Phenylalanine metabolism',
    'ko00380': 'Tryptophan metabolism',
    'ko00400': 'Phenylalanine, tyrosine and tryptophan biosynthesis',
    'ko00410': 'beta-Alanine metabolism',
    'ko00430': 'Taurine and hypotaurine metabolism',
    'ko00450': 'Selenocompound metabolism',
    'ko00480': 'Glutathione metabolism',
    'ko00500': 'Starch and sucrose metabolism',
    'ko00520': 'Amino sugar and nucleotide sugar metabolism',
    'ko00561': 'Glycerolipid metabolism',
    'ko00562': 'Inositol phosphate metabolism',
    'ko00564': 'Glycerophospholipid metabolism',
    'ko00565': 'Ether lipid metabolism',
    'ko00590': 'Arachidonic acid metabolism',
    'ko00591': 'Linoleic acid metabolism',
    'ko00592': 'alpha-Linolenic acid metabolism',
    'ko00600': 'Sphingolipid metabolism',
    'ko00620': 'Pyruvate metabolism',
    'ko00630': 'Glyoxylate and dicarboxylate metabolism',
    'ko00640': 'Propanoate metabolism',
    'ko00650': 'Butanoate metabolism',
    'ko00670': 'One carbon pool by folate',
    'ko00730': 'Thiamine metabolism',
    'ko00740': 'Riboflavin metabolism',
    'ko00750': 'Vitamin B6 metabolism',
    'ko00760': 'Nicotinate and nicotinamide metabolism',
    'ko00770': 'Pantothenate and CoA biosynthesis',
    'ko00780': 'Biotin metabolism',
    'ko00785': 'Lipoic acid metabolism',
    'ko00790': 'Folate biosynthesis',
    'ko00830': 'Retinol metabolism',
    'ko00860': 'Porphyrin metabolism',
    'ko00900': 'Terpenoid backbone biosynthesis',
    'ko01040': 'Biosynthesis of unsaturated fatty acids',
    'ko01100': 'Metabolic pathways',
    'ko01200': 'Carbon metabolism',
    'ko01210': '2-Oxocarboxylic acid metabolism',
    'ko01212': 'Fatty acid metabolism',
    'ko01230': 'Biosynthesis of amino acids',
    'ko04080': 'Neuroactive ligand-receptor interaction',
    'ko04973': 'Carbohydrate metabolism (disease)',
    'ko04974': 'Protein digestion and absorption',
    'ko04975': 'Fat digestion and absorption',
    'ko04976': 'Vitamin digestion and absorption',
}


# ---- 内部工具函数 ----

def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正"""
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


def _build_kegg_pathway_map(metabolite_info: pd.DataFrame) -> dict[str, set[str]]:
    """从 kegg_map 列构建 {pathway_id: {metabolite_Index}} 映射"""
    pathway_map: dict[str, set[str]] = {}
    kegg_col = metabolite_info.get('kegg_map')
    if kegg_col is None:
        return pathway_map
    for idx, val in kegg_col.items():
        if pd.isna(val) or val == '-':
            continue
        for pid in str(val).split(','):
            pid = pid.strip()
            if pid and pid != '-':
                pathway_map.setdefault(pid, set()).add(idx)
    return pathway_map


def _fetch_kegg_pathway_names() -> dict[str, str]:
    """从 KEGG API 获取全部通路名称（koXXXXX -> name），带缓存"""
    try:
        with urllib.request.urlopen('https://rest.kegg.jp/list/pathway', timeout=15) as resp:
            text = resp.read().decode('utf-8')
        names = {}
        for line in text.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) == 2:
                # mapXXXXX -> koXXXXX
                ko_id = 'ko' + parts[0].lstrip('map').lstrip('path:')
                # 去掉行首的 'path:map' 前缀
                ko_id = 'ko' + parts[0].split('map')[-1]
                names[ko_id] = parts[1]
        return names
    except Exception:
        return {}


def _build_custom_pathway_map(mapping_file: str | Path) -> dict[str, set[str]]:
    """从外部文件构建 {pathway_name: {metabolite_Index}} 映射

    JSON: {"pathway_name": ["INDEX1", "INDEX2", ...], ...}
    CSV:  pathway_name,metabolite_index（每行一条）
    """
    path = Path(mapping_file)
    if path.suffix == '.json':
        with open(path, encoding='utf-8') as f:
            raw = json.load(f)
        return {name: set(ids) for name, ids in raw.items()}
    df = pd.read_csv(path, header=None, names=['pathway', 'metabolite'])
    pm: dict[str, set[str]] = {}
    for _, row in df.iterrows():
        pm.setdefault(str(row['pathway']), set()).add(str(row['metabolite']))
    return pm


# ---- 差异代谢物分析 ----

def identify_differential_metabolites(
        data: pd.DataFrame,
        group1_cols: list[str], group2_cols: list[str],
        method: str = 'ttest',
        p_threshold: float = 0.1,
        fc_threshold: float = 1.5,
) -> pd.DataFrame:
    """差异代谢物分析

    Args:
        data: 代谢物(Index) × 样本 丰度矩阵
        group1_cols: 组1 列名
        group2_cols: 组2 列名
        method: 'ttest' 或 'mannwhitney'
        p_threshold: 原始 p 值阈值
        fc_threshold: fold change 阈值

    Returns:
        DataFrame(index=代谢物Index)，含 log2FC / p_value / adj_p_value / significant
    """
    g1 = data[group1_cols].to_numpy(dtype=float)
    g2 = data[group2_cols].to_numpy(dtype=float)
    n = data.shape[0]

    # 过滤两组均全 NaN 的代谢物
    valid_mask = ~(np.all(np.isnan(g1), axis=1) & np.all(np.isnan(g2), axis=1))
    g1, g2 = g1[valid_mask], g2[valid_mask]
    valid_index = data.index[valid_mask]
    n = len(valid_index)

    # log2FC
    g1_mean = np.nanmean(g1, axis=1)
    g2_mean = np.nanmean(g2, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        log2fc = np.log2((g2_mean + 1) / (g1_mean + 1))

    # 逐代谢物统计检验
    p_values = np.full(n, np.nan)
    for i in range(n):
        x = g1[i, ~np.isnan(g1[i])]
        y = g2[i, ~np.isnan(g2[i])]
        if len(x) < 2 or len(y) < 2:
            continue
        try:
            if method == 'mannwhitney':
                _, p_values[i] = stats.mannwhitneyu(x, y, alternative='two-sided')
            else:
                _, p_values[i] = stats.ttest_ind(x, y, equal_var=False)
        except ValueError:
            pass

    # FDR 校正
    valid = ~np.isnan(p_values)
    adj_p = np.full(n, np.nan)
    if valid.any():
        adj_p[valid] = _benjamini_hochberg(p_values[valid])

    result = pd.DataFrame({
        'log2FC': log2fc, 'p_value': p_values, 'adj_p_value': adj_p,
        'group1_mean': g1_mean, 'group2_mean': g2_mean,
    }, index=valid_index)
    result['significant'] = (
        (result['p_value'] < p_threshold)
        & (np.abs(result['log2FC']) >= np.log2(fc_threshold))
    )
    return result


# ---- 通路富集 ----

def pathway_enrichment(
        significant_ids: set[str],
        all_measured_ids: set[str],
        pathway_map: dict[str, set[str]],
        pathway_names: dict[str, str] | None = None,
        min_pathway_size: int = 3,
        metabolite_info: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """超几何检验通路富集

    Args:
        significant_ids: 显著差异代谢物 Index 集合
        all_measured_ids: 背景（所有检测到的）代谢物 Index 集合
        pathway_map: {pathway_id: {metabolite_index_set}}
        pathway_names: {pathway_id: 显示名称}
        min_pathway_size: 通路中最少代谢物数（过滤过小的通路）
        metabolite_info: 代谢物元信息 DataFrame（用于查找英文名）

    Returns:
        DataFrame，按 p_value 升序排列
    """
    if pathway_names is None:
        pathway_names = KEGG_PATHWAY_NAMES

    # 检查缺失名称，从 KEGG API 补全
    missing_ids = set(pathway_map.keys()) - set(pathway_names.keys())
    if missing_ids:
        api_names = _fetch_kegg_pathway_names()
        for pid in missing_ids:
            if pid in api_names:
                pathway_names[pid] = api_names[pid]

    # 代谢物 Index -> 英文名 映射
    compound_map = {}
    if metabolite_info is not None and 'Compounds' in metabolite_info.columns:
        compound_map = metabolite_info['Compounds'].to_dict()

    n_total = len(all_measured_ids)
    n_sig = len(significant_ids)
    records = []

    for pid, metabolites in pathway_map.items():
        in_bg = metabolites & all_measured_ids
        size = len(in_bg)
        if size < min_pathway_size:
            continue
        overlap = metabolites & significant_ids
        k = len(overlap)
        if k == 0:
            continue
        p_val = stats.hypergeom.sf(k - 1, n_total, size, n_sig)
        sorted_overlap = sorted(overlap)
        records.append({
            'pathway_id': pid,
            'pathway_name': pathway_names.get(pid, pid),
            'pathway_name_zh': KEGG_PATHWAY_NAMES_ZH.get(pid, pathway_names.get(pid, pid)),
            'pathway_size': size,
            'significant_count': k,
            'significant_ratio': k / n_sig,
            'coverage': k / size,
            'p_value': p_val,
            'significant_metabolites': '; '.join(sorted_overlap),
            'significant_metabolite_names': '; '.join(
                compound_map.get(idx, idx) for idx in sorted_overlap),
        })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df['adj_p_value'] = _benjamini_hochberg(df['p_value'].to_numpy())
    df['significant'] = df['adj_p_value'] < 0.1
    return df.sort_values('p_value').reset_index(drop=True)


# ---- 整合分析类 ----

class MetabolomicsPathwayEnrichment:
    """差异分析 + 通路富集 一站式分析

    Usage:
        analyzer = MetabolomicsPathwayEnrichment(ds)
        analyzer.run('AD', 'AD+Sema', tissue='cortex')
        analyzer.export('result.xlsx')
    """

    def __init__(self, dataset):
        self.ds = dataset
        self.diff_result: pd.DataFrame | None = None
        self.enrichment_result: pd.DataFrame | None = None

    def run(
            self,
            group1: str, group2: str,
            tissue: str,
            method: str = 'ttest',
            p_threshold: float = 0.05,
            fc_threshold: float = 1.5,
            database: str = 'kegg',
            custom_mapping_file: str | None = None,
            min_pathway_size: int = 3,
    ):
        """运行完整流程：差异分析 → 通路富集"""
        g1_data = self.ds.get_data_by_group(group1, tissue=tissue)
        g2_data = self.ds.get_data_by_group(group2, tissue=tissue)
        if g1_data.empty or g2_data.empty:
            raise ValueError(
                f'分组数据为空: {group1}({g1_data.shape[1]}列), '
                f'{group2}({g2_data.shape[1]}列)')

        # 取两组共有代谢物
        common = g1_data.index.intersection(g2_data.index)
        g1_data, g2_data = g1_data.loc[common], g2_data.loc[common]

        # 差异分析
        merged = pd.concat([g1_data, g2_data], axis=1)
        self.diff_result = identify_differential_metabolites(
            merged, list(g1_data.columns), list(g2_data.columns),
            method=method, p_threshold=p_threshold, fc_threshold=fc_threshold)

        # 附加代谢物元信息
        meta = self.ds.metabolite_info
        meta_cols = [c for c in ['Compounds', '物质', 'Class I', 'Class II', 'kegg_map']
                     if c in meta.columns]
        self.diff_result = self.diff_result.join(meta[meta_cols], how='left')

        # 构建通路映射
        if database == 'custom':
            if custom_mapping_file is None:
                raise ValueError('custom 数据库需要提供 custom_mapping_file')
            pm = _build_custom_pathway_map(custom_mapping_file)
            names = {}
        else:
            pm = _build_kegg_pathway_map(meta)
            names = KEGG_PATHWAY_NAMES

        # 富集分析
        sig_ids = set(self.diff_result[self.diff_result['significant']].index)
        valid_mask = merged.notna().any(axis=1)
        all_ids = set(merged.index[valid_mask])

        self.enrichment_result = pathway_enrichment(
            sig_ids, all_ids, pm, names,
            min_pathway_size=min_pathway_size, metabolite_info=meta)

        n_sig = len(sig_ids)
        n_enriched = (len(self.enrichment_result)
                      if self.enrichment_result is not None else 0)
        print(f'{tissue} | {group1} vs {group2}: '
              f'{n_sig} 差异代谢物 / {n_enriched} 富集通路')
        return self

    def export(self, output_path: str | Path):
        """导出到 xlsx（差异代谢物 + 通路富集 两个 sheet）"""
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            if self.diff_result is not None:
                self.diff_result.sort_values('adj_p_value').to_excel(
                    writer, sheet_name='差异代谢物')
            if self.enrichment_result is not None and not self.enrichment_result.empty:
                self.enrichment_result.to_excel(writer, sheet_name='通路富集')
        print(f'已导出: {path}')
