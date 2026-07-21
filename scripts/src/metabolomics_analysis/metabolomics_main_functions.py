# -*- coding: utf-8 -*-
from pathlib import Path

import numpy as np
import pandas as pd

from ..inventory import MetabolomicsRunningMode, MetabolomicsExperimentName
from .config import default_data_model_name

# 数据集目录映射：experiment_name -> 相对于 scripts/data/ 的路径
_DATASET_DIR_MAP = {
    MetabolomicsExperimentName.data_20260412_metabolomics: 'data_20260412_metabolomics',
}

# 缓存目录
CACHE_DIR = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'cache'


def _cache_path(panel: str, tissue: str = '', suffix: str = '') -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    parts = [panel]
    if tissue:
        parts.append(tissue)
    tag = f'_{suffix}' if suffix else ''
    return CACHE_DIR / f'{"_".join(parts)}{tag}.pkl.gz'


def _load_dataset(experiment_name, msi_level='3'):
    """加载数据集，返回 MetabolomicDataset（经 common_data_loader 分流）"""
    from scripts.data.data_loader import common_data_loader, DataSource
    content, _ = common_data_loader(DataSource.data_20260412_metabolomics, min_msi_level=msi_level)
    return content


def metabolomics_analysis_main(
        running_mode, experiment_name=None, suffix=None,
        tissue=None, method='ttest', p_threshold=0.05, fc_threshold=1.5,
        n_permutations=200, vip_threshold=1.0, msi_level='3'):
    if experiment_name is None:
        experiment_name = default_data_model_name

    ds = _load_dataset(experiment_name, msi_level)

    if running_mode == MetabolomicsRunningMode.summary:
        _run_summary(ds, tissue)
    elif running_mode == MetabolomicsRunningMode.export:
        _run_export(ds, tissue, suffix)
    elif running_mode == MetabolomicsRunningMode.enrichment:
        _run_enrichment(ds, tissue, method, p_threshold, fc_threshold, suffix)
    elif running_mode == MetabolomicsRunningMode.plsda:
        _run_plsda(ds, tissue, suffix)
    elif running_mode == MetabolomicsRunningMode.oplsda:
        _run_oplsda(ds, tissue, n_permutations, suffix)
    elif running_mode == MetabolomicsRunningMode.clustered_heatmap:
        _run_clustered_heatmap(ds, tissue, method, p_threshold, fc_threshold,
                               vip_threshold, suffix)
    elif running_mode == MetabolomicsRunningMode.pathway_bubble:
        _run_pathway_bubble(ds, tissue, method, p_threshold, fc_threshold, suffix)
    elif running_mode == MetabolomicsRunningMode.splot:
        _run_splot(ds, tissue, method, suffix)
    elif running_mode == MetabolomicsRunningMode.phenotype_corr:
        _run_phenotype_corr(ds, tissue, method, suffix)
    elif running_mode == MetabolomicsRunningMode.volcano:
        _run_volcano(ds, tissue, method, p_threshold, fc_threshold, suffix)
    else:
        raise ValueError(f'Unknown running mode: {running_mode}')


# ============================================================
# 汇总 & 导出（无需缓存）
# ============================================================

def _run_summary(ds, tissue=None):
    """打印数据集概览"""
    print(ds.summary())
    print()

    if tissue is not None:
        data = ds.get_data_by_tissue(tissue)
        print(f'--- {tissue} 数据 ---')
        print(f'样本数: {data.shape[1]}')
        print(f'代谢物数: {data.shape[0]}')
        print(f'缺失值比例: {data.isna().mean().mean():.2%}')
        print(f'样本列: {list(data.columns)}')
    else:
        for t in ds.list_tissues():
            data = ds.get_data_by_tissue(t)
            na_ratio = data.isna().mean().mean()
            print(f'{t}: {data.shape[1]} 样本, {data.shape[0]} 代谢物, 缺失 {na_ratio:.2%}')


def _run_export(ds, tissue=None, suffix=None):
    """导出整理后的数据为 xlsx"""
    output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics'
    tag = f'_{suffix}' if suffix else ''
    output_dir.mkdir(parents=True, exist_ok=True)

    tissues = [tissue] if tissue else ds.list_tissues()
    for t in tissues:
        data = ds.get_data_by_tissue(t)
        info = ds.metabolite_info

        output_path = output_dir / f'{t}{tag}.xlsx'
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            merged = info.join(data, how='inner')
            merged.to_excel(writer, sheet_name=t)
        print(f'已导出: {output_path}')

    print(f'\n全部导出完成，输出目录: {output_dir}')


# ============================================================
# 以下函数均支持 pickle 缓存
# ============================================================

def _run_enrichment(ds, tissue, method, p_threshold, fc_threshold, suffix):
    """差异代谢物通路富集分析（默认 ad vs semaglutide）"""
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .pathway_enrichment import MetabolomicsPathwayEnrichment

    group1, group2 = Keyword.ad, Keyword.sema

    output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'enrichment'
    output_dir.mkdir(parents=True, exist_ok=True)

    tag = f'_{suffix}' if suffix else ''
    tissues = [tissue] if tissue else ds.list_tissues()
    for t in tissues:
        analyzer = MetabolomicsPathwayEnrichment(ds)
        analyzer.run(group1, group2, tissue=t,
                     method=method, p_threshold=p_threshold, fc_threshold=fc_threshold)
        output_path = output_dir / f'{t}_{group1}_vs_{group2}{tag}.xlsx'
        analyzer.export(output_path)


def _run_plsda(ds, tissue, suffix):
    """PLS-DA 分析（默认 ad vs semaglutide）"""
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .plsda_analysis import run_plsda, plot_plsda, export_plsda

    group1, group2 = Keyword.ad, Keyword.sema

    output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'plsda'
    output_dir.mkdir(parents=True, exist_ok=True)

    tag = f'_{suffix}' if suffix else ''
    tissues = [tissue] if tissue else ds.list_tissues()
    for t in tissues:
        g1_data = ds.get_data_by_group(group1, tissue=t)
        g2_data = ds.get_data_by_group(group2, tissue=t)
        if g1_data.empty or g2_data.empty:
            print(f'{t}: 数据为空，跳过')
            continue
        common_idx = g1_data.index.intersection(g2_data.index)
        merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)

        result = run_plsda(
            merged, list(g1_data.columns), list(g2_data.columns),
            group1_name=group1, group2_name=group2)
        result.metabolite_info = ds.metabolite_info
        print(f'{t}: R²X={result.r2x:.3f}, R²Y={result.r2y:.3f}')

        base = f'{t}_{group1}_vs_{group2}{tag}'
        export_plsda(result, output_dir / f'{base}.xlsx')
        plot_plsda(result, output_dir / f'{base}.png')


def _run_oplsda_single(ds, tissue, suffix, *, force=False):
    """per-tissue 完整 OPLSDAResult 中间缓存（含 scores+p_loadings+p_corr+vip）。
    供 oplsda/clustered_heatmap/splot 共用，避免重复跑 run_oplsda。
    """
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .splot import run_splot_oplsda

    group1, group2 = Keyword.ad, Keyword.sema
    cache = _cache_path('oplsda', tissue=tissue, suffix=suffix)
    if cache.exists() and not force:
        result = pickle_load(cache)
        print(f'{tissue}: 从缓存加载 OPLS-DA')
        return result

    result = run_splot_oplsda(ds, tissue, group1, group2)
    if result is None:
        return None
    pickle_save(result, cache)
    print(f'{tissue}: 已缓存 OPLS-DA 到 {cache.name}')
    return result


def _run_oplsda(ds, tissue, n_permutations, suffix, *, force=False):
    """OPLS-DA 分析 + 置换检验（带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .plsda_analysis import (
        permutation_test_oplsda, plot_oplsda, plot_oplsda_grid,
        plot_permutation_test, plot_permutation_grid, export_oplsda,
        run_combined_pca,
    )

    group1, group2 = Keyword.ad, Keyword.sema
    tag = f'_{suffix}' if suffix else ''

    # 尝试读取缓存
    cache = _cache_path('oplsda_all', suffix=suffix)
    if cache.exists() and not force:
        cached = pickle_load(cache)
        results = cached['results']
        perm_results = cached.get('perm_results', {})
        print(f'OPLS-DA: 从缓存加载 ({len(results)} 个组织)')
    else:
        output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'oplsda'
        output_dir.mkdir(parents=True, exist_ok=True)

        tissues = [tissue] if tissue else ds.list_tissues()
        results, perm_results, tissue_raw = [], {}, {}

        for t in tissues:
            g1_data = ds.get_data_by_group(group1, tissue=t)
            g2_data = ds.get_data_by_group(group2, tissue=t)
            if g1_data.empty or g2_data.empty:
                print(f'{t}: 数据为空，跳过')
                continue
            if g1_data.shape[1] < 3 or g2_data.shape[1] < 3:
                print(f'{t}: 样本数不足，跳过')
                continue

            common_idx = g1_data.index.intersection(g2_data.index)
            merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)
            g1_cols, g2_cols = list(g1_data.loc[common_idx].columns), list(g2_data.loc[common_idx].columns)
            tissue_raw[t] = (merged, g1_cols, g2_cols)

            # OPLSDAResult 从 per-tissue 缓存取（_run_oplsda_single），避免重复跑 run_oplsda
            result = _run_oplsda_single(ds, t, suffix)
            if result is None:
                continue
            results.append(result)
            print(f'{t}: R²Y={result.r2y:.3f}, Q²={result.q2:.3f}')

            base = f'{t}_{group1}_vs_{group2}{tag}'
            plot_oplsda(result, output_dir / f'{base}_oplsda.png')
            export_oplsda(result, output_dir / f'{base}_oplsda.xlsx')

            if n_permutations > 0:
                perm = permutation_test_oplsda(merged, g1_cols, g2_cols,
                                               n_permutations=n_permutations)
                perm_results[t] = perm
                plot_permutation_test(perm, tissue=t,
                                      output_path=output_dir / f'{base}_perm.png')
                print(f'  置换检验: p(R²Y)={perm["p_r2y"]:.4f}, p(Q²)={perm["p_q2"]:.4f}')

        # 保存缓存
        pickle_save({'results': results, 'perm_results': perm_results}, cache)
        print(f'OPLS-DA: 已缓存到 {cache.name}')

        # 生成组合图
        if len(results) > 1:
            plot_oplsda_grid(results, output_dir / f'oplsda_grid{tag}.png')
        if perm_results and len(perm_results) > 1:
            plot_permutation_grid(perm_results, output_dir / f'permutation_grid{tag}.png')
        if len(tissue_raw) > 1:
            run_combined_pca(tissue_raw, group1, group2,
                             output_dir / f'combined_pca{tag}.png')

    return {'results': results, 'perm_results': perm_results}


def _run_clustered_heatmap(ds, tissue, method, p_threshold, fc_threshold,
                           vip_threshold, suffix, *, force=False):
    """聚类热图（带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .clustered_heatmap import (
        run_clustered_heatmap, plot_clustered_heatmap,
        plot_clustered_heatmap_grid, export_clustered_heatmap,
    )

    group1, group2 = Keyword.ad, Keyword.sema
    tag = f'_{suffix}' if suffix else ''
    tissues = [tissue] if tissue else ['hippocampus', 'cortex']

    # 逐组织缓存
    results = []
    for t in tissues:
        cache = _cache_path('heatmap', tissue=t, suffix=suffix)
        if cache.exists() and not force:
            r = pickle_load(cache)
            print(f'{t}: 从缓存加载聚类热图')
        else:
            diff = _run_differential(ds, t, method, suffix)
            oplsda_result = _run_oplsda_single(ds, t, suffix)
            r = run_clustered_heatmap(
                ds, t, group1, group2,
                method=method, p_threshold=p_threshold,
                fc_threshold=fc_threshold, vip_threshold=vip_threshold,
                diff_result=diff, oplsda_result=oplsda_result)
            if r is None:
                continue
            pickle_save(r, cache)
            print(f'{t}: 已缓存聚类热图到 {cache.name}')

            # 同时输出独立图片
            output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'clustered_heatmap'
            output_dir.mkdir(parents=True, exist_ok=True)
            base = f'{t}_{group1}_vs_{group2}{tag}'
            plot_clustered_heatmap(
                r['z_matrix'], r['row_order'], r['col_order'],
                r['row_linkage'], r['col_linkage'],
                r['metabolite_labels'], r['sample_labels'], r['groups'],
                output_path=output_dir / f'{base}_clustered_heatmap.png',
                tissue=t)
            export_clustered_heatmap(r, output_dir / f'{base}_clustered_heatmap.xlsx')

        results.append(r)

    if len(results) > 1:
        output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'clustered_heatmap'
        plot_clustered_heatmap_grid(results, output_dir / f'clustered_heatmap_grid{tag}.png')

    return results


def _run_differential(ds, tissue, method, suffix, *, force=False):
    """per-tissue 差异分析中间缓存（log2FC/p_value/adj_p_value 等）。

    用宽松阈值(p=1.0, fc=1.0)缓存原始统计，各 panel 自行用面板阈值筛选；
    供 volcano/clustered_heatmap/pathway_bubble 共用，避免重复计算 identify_differential_metabolites。
    """
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .pathway_enrichment import identify_differential_metabolites

    group1, group2 = Keyword.ad, Keyword.sema
    cache = _cache_path('diff', tissue=tissue, suffix=suffix)
    if cache.exists() and not force:
        diff = pickle_load(cache)
        print(f'{tissue}: 从缓存加载差异分析')
        return diff

    g1_data = ds.get_data_by_group(group1, tissue=tissue)
    g2_data = ds.get_data_by_group(group2, tissue=tissue)
    if g1_data.empty or g2_data.empty:
        return None
    common_idx = g1_data.index.intersection(g2_data.index)
    merged = pd.concat([g1_data.loc[common_idx], g2_data.loc[common_idx]], axis=1)
    diff = identify_differential_metabolites(
        merged, list(g1_data.loc[common_idx].columns), list(g2_data.loc[common_idx].columns),
        method=method, p_threshold=1.0, fc_threshold=1.0)
    pickle_save(diff, cache)
    print(f'{tissue}: 已缓存差异分析到 {cache.name}')
    return diff


def _run_volcano(ds, tissue, method, p_threshold, fc_threshold, suffix, *, force=False):
    """火山图数据（从 _run_differential 组装，带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load

    cache = _cache_path('volcano', tissue=tissue or 'all', suffix=suffix)
    if cache.exists() and not force:
        result = pickle_load(cache)
        print(f'火山图: 从缓存加载')
        return result

    tissues = [tissue] if tissue else ds.list_tissues()
    name_map = ds.metabolite_info['Compounds'].to_dict() if 'Compounds' in ds.metabolite_info.columns else {}
    result = {}
    for t in tissues:
        diff = _run_differential(ds, t, method, suffix)
        if diff is None or diff.empty:
            continue
        # 过滤 p_value/log2FC 为 NaN/Inf 的代谢物（ttest 在数据不足/全相同时返回 NaN，否则火山图 y_lim 会变 NaN）
        valid = (
            diff['p_value'].notna() & np.isfinite(diff['p_value'])
            & diff['log2FC'].notna() & np.isfinite(diff['log2FC']))
        diff = diff[valid]
        if diff.empty:
            continue
        log2fc = diff['log2FC'].values
        neg_log10p = -np.log10(diff['p_value'].clip(lower=1e-300)).values
        adj_p = diff['adj_p_value'].values if 'adj_p_value' in diff.columns else diff['p_value'].values
        names = [name_map.get(idx, str(idx)) for idx in diff.index]

        result[t] = {
            'log2fc': log2fc,
            'neg_log10p': neg_log10p,
            'p_values': diff['p_value'].values,
            'adj_p': adj_p,
            'thresholds': (fc_threshold, p_threshold),
            'names': names,
            'index': list(diff.index),
        }
        print(f'{t}: {len(log2fc)} 代谢物, '
              f'{(np.abs(log2fc) > np.log2(fc_threshold)).sum()} 超FC阈值, '
              f'{(diff["p_value"] < p_threshold).sum()} 超p阈值')

    pickle_save(result, cache)
    print(f'火山图: 已缓存到 {cache.name}')
    return result


def _run_pathway_bubble(ds, tissue, method, p_threshold, fc_threshold, suffix, *, force=False):
    """多组织通路富集 + 拓扑分析（带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .multi_tissue_pathway import (
        run_multi_tissue_pathway, plot_pathway_bubble,
        plot_metaboanalyst_scatter, export_multi_tissue_pathway,
        TISSUE_ORDER,
    )

    group1, group2 = Keyword.ad, Keyword.sema
    tag = f'_{suffix}' if suffix else ''

    cache = _cache_path('pathway_bubble', suffix=suffix)
    if cache.exists() and not force:
        results = pickle_load(cache)
        print(f'通路气泡图: 从缓存加载 ({len(results)} 行)')
        return results

    tissues = [tissue] if tissue else TISSUE_ORDER

    # 预算各组织 diff（从 _run_differential 缓存），避免 run_multi_tissue_pathway 内部重复计算
    diff_results = {}
    for t in tissues:
        diff = _run_differential(ds, t, method, suffix)
        if diff is not None:
            diff_results[t] = diff

    results = run_multi_tissue_pathway(
        ds, tissues, group1, group2,
        method=method, p_threshold=p_threshold, fc_threshold=fc_threshold,
        diff_results=diff_results)

    if results.empty:
        print('无通路富集结果')
        return results

    pickle_save(results, cache)
    print(f'通路气泡图: 已缓存到 {cache.name}')

    # 同时输出独立图片
    output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'pathway_bubble'
    output_dir.mkdir(parents=True, exist_ok=True)
    base = f'{group1}_vs_{group2}{tag}'

    plot_pathway_bubble(results, output_dir / f'pathway_bubble_direction{tag}.png',
                        color_by='direction')
    plot_pathway_bubble(results, output_dir / f'pathway_bubble_significance{tag}.png',
                        color_by='significance')

    for t in tissues:
        t_data = results[results['tissue'] == t]
        if t_data.empty:
            continue
        plot_metaboanalyst_scatter(
            t_data, output_dir / f'{t}_{base}_metaboanalyst{tag}.png', tissue=t)

    export_multi_tissue_pathway(results, output_dir / f'multi_tissue_pathway{tag}.xlsx')

    return results


def _run_splot(ds, tissue, method, suffix, *, force=False):
    """S-plot 分析（带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from .splot import (
        identify_candidates,
        plot_splot, plot_vip_bar, export_splot,
    )

    group1, group2 = Keyword.ad, Keyword.sema
    tag = f'_{suffix}' if suffix else ''
    tissues = [tissue] if tissue else ['hippocampus', 'cortex']

    results = {}
    for t in tissues:
        cache = _cache_path('splot', tissue=t, suffix=suffix)
        if cache.exists() and not force:
            cached = pickle_load(cache)
            results[t] = cached
            print(f'{t}: 从缓存加载 S-plot')
            continue

        # OPLSDAResult 从 per-tissue 缓存取（_run_oplsda_single），避免重复跑 run_oplsda
        result = _run_oplsda_single(ds, t, suffix)
        if result is None or result.p_corr is None:
            print(f'{t}: 无法计算 S-plot 数据')
            continue

        candidates = identify_candidates(result, p_corr_threshold=0.5, vip_threshold=1.0)
        n_cand = candidates['candidate'].sum() if not candidates.empty else 0
        print(f'{t}: R2Y={result.r2y:.3f}, Q2={result.q2:.3f}, '
              f'{n_cand} 候选驱动代谢物')

        cached = {'oplsda_result': result, 'candidates': candidates}
        pickle_save(cached, cache)
        print(f'{t}: 已缓存 S-plot 到 {cache.name}')
        results[t] = cached

        # 同时输出独立图片
        output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'splot'
        output_dir.mkdir(parents=True, exist_ok=True)
        base = f'{t}_{group1}_vs_{group2}{tag}'
        plot_splot(result, candidates, output_dir / f'{base}_splot.png')
        plot_vip_bar(candidates, output_dir / f'{base}_vip_bar.png', tissue=t)
        export_splot(candidates, output_dir / f'{base}_splot.xlsx')

    return results


def _run_phenotype_corr(ds, tissue, method, suffix, *, force=False):
    """表型-代谢物相关分析（带缓存）"""
    from scripts.src.common_parallel_solver.config import pickle_save, pickle_load
    from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword
    from scripts.data.data_loader import common_data_loader, DataSource
    from .phenotype_correlation import (
        run_phenotype_correlation, extract_phenotype_metrics,
        get_metabolite_matrix, get_key_metabolites,
        plot_correlation_heatmap, plot_scatter_matrix, export_correlation,
    )

    group1, group2 = Keyword.ad, Keyword.sema
    tag = f'_{suffix}' if suffix else ''

    cache = _cache_path('phenotype_corr', suffix=suffix)
    if cache.exists() and not force:
        cached = pickle_load(cache)
        print(f'表型相关: 从缓存加载')
        return cached

    behavioral_ds, _ = common_data_loader(DataSource.metabolomics_behavioral)
    tissues = [tissue] if tissue else ['hippocampus', 'serum']

    correlations = run_phenotype_correlation(
        ds, behavioral_ds, tissues, group1, group2, method=method)

    pheno = extract_phenotype_metrics(behavioral_ds)

    cached = {'correlations': correlations, 'pheno': pheno}
    pickle_save(cached, cache)
    print(f'表型相关: 已缓存到 {cache.name}')

    # 同时输出独立图片
    output_dir = Path(__file__).parents[2] / 'output' / 'metabolomics' / 'phenotype_corr'
    output_dir.mkdir(parents=True, exist_ok=True)

    for tissue_name, corr_df in correlations.items():
        all_mets = corr_df['metabolite_id'].unique().tolist()
        met_matrix = get_metabolite_matrix(ds, tissue_name, all_mets)

        base = f'{tissue_name}_{group1}_vs_{group2}{tag}'
        tissue_zh = {'hippocampus': '海马', 'cortex': '皮层', 'serum': '血清'}.get(tissue_name, tissue_name)

        plot_correlation_heatmap(
            corr_df, output_dir / f'{base}_heatmap.png',
            title=f'{tissue_zh} — Metabolite-Phenotype Correlations')
        plot_scatter_matrix(
            corr_df, met_matrix, pheno,
            output_dir / f'{base}_scatter.png', tissue=tissue_name)

    if correlations:
        export_correlation(correlations, output_dir / f'phenotype_correlation{tag}.xlsx')

    return cached
