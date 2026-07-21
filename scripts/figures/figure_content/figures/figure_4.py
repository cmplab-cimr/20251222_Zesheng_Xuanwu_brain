# -*- coding: utf-8 -*-
"""Figure 4: AD + Semaglutide 多组织系统性代谢重编程

Panel A: 多组织 OPLS-DA 网格
Panel B: 差异代谢物聚类热图
Panel C: 火山图
Panel D: 跨器官通路富集点图
Panel E: 关键代谢物跨组织柱状图
Panel F: S-plot
Panel G: 代谢物–表型相关
"""

import numpy as np
import pandas as pd

from ...common.config import ParameterName, Figure as GeneralFigure, Vector, calculate_center_bottom_offset
from ...common.elements import Elements
from ...common.common_figure_materials import (
    CommonColorDict as CCD, gray_teal_diverging_cmap, rdbu_diverging_cmap, is_lipid_class, lipid_class_color)
from scripts.src.common_parallel_solver.config import pickle_load
from scripts.src.metabolomics_analysis.metabolomics_main_functions import (
    _cache_path, _load_dataset,
)
from scripts.src.inventory import MetabolomicsExperimentName
from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword

Subfigure = Elements.Subfigure

# 分组颜色: AD=灰, Semaglutide=青
GROUP_COLORS = {Keyword.ad: CCD.gray, Keyword.sema: CCD.teal}

# 火山图语义配色: 上调=青, 下调=灰, 不显著=浅灰
VOLCANO_COLOR_DICT = {
    ParameterName.positive: CCD.teal,
    ParameterName.negative: CCD.gray,
    ParameterName.medium_data: CCD.light_gray,
}

# 表型名中文 -> 英文映射 (Panel G 热图 x 轴标签): 中文键用于匹配缓存数据, 英文用于显示
PHENO_NAME_EN = {
    'MWM_learn_rate': 'MWM learn-rate',
    'NOR_DI_鼻尖': 'NOR DI (nose)',
    'NOR_DI_中心': 'NOR DI (center)',
    'OF_区域时间': 'OF center time',
    'CCS': 'Cognitive score',
}


# Panel E 跨组织柱状图展示的代谢物 ID（按故事线手动调整；为空则 Panel E 留空）
# 候选为跨组织系统性变化（多组织 VIP 高）的代谢物, 覆盖 9-10/10 组织
CROSS_TISSUE_BAR_METABOLITES = [
    'MEDN4139',       # 5-Hydroxy-2-methyl-indole-3-carboxylic acid — 色氨酸/吲哚 (肠道-脑轴), 10 组织
    'MEDP2752',       # 3-Indole carboxylic acid glucuronide — 色氨酸/吲哚葡萄糖醛酸, 10 组织
    'MEDP4197',       # 5-Hydroxy-6-methoxyindole glucuronide — 吲哚代谢, 10 组织
    'MEDN0471',       # MHPG sulfate — 去甲肾上腺素代谢 (神经递质), 9 组织
    'MEDN3905',       # N2-(1-Carboxyethyl)-L-arginine — 氨基酸衍生, 10 组织
    'MEDP1493N',      # N6-Methyladenosine (m6A) — 表观转录修饰, 9 组织
]

# Panel B 聚类热图主图最多显示的代谢物行数（按 VIP 排名截取，缓存仍存全量）
HEATMAP_TOP_N = 60

# Panel B 右侧强调的代谢物（按故事线手动填，可为 metabolite index 或 Compounds 名称）
# 仅命中主图 top-N 的代谢物才画连线；候选为 hippocampus top-60 内 VIP 最高的脂质
HEATMAP_HIGHLIGHT_METABOLITES = [
    'MEDN1430*207',   # Prostaglandin A2 (FA) — 炎症介质
    'MEDP1171P',      # Butenoyl-PAF (GP) — 血小板活化因子, 神经炎症
    'MEDN1414N',      # 12S-HHT (FA) — 血栓素通路
    'MEDN1597',       # PE lyso alkenyl 16:0 (GP) — 缩醛磷脂, AD 相关
    'MEDP2746',       # LPC(O-18:3) (GP) — 醚键溶血磷脂
    'MADN3058',       # 5,6-DiHETrE (FA) — EET 衍生物, 血管/神经
    'MEDN0752',       # 15-HEPE (FA) — EPA 抗炎衍生物
    'ME0057544HN',    # PC(22:4/15:0) (GP) — 多不饱和磷脂
]

# Panel D/E/F 共用代谢物：驱动 GP/FA 富集的已注释脂质代谢物
# 初版由 select_lipid_metabolites_for_DEF.py 筛选后填入，用户后续手动调整
# 覆盖胆碱通路（胆碱、磷酸胆碱、GPC）、lysoPC/lysoPE、3-磷酸甘油、关键 FFA、酰基肉碱等
LIPID_METABOLITES_FOR_DEF = [
    'MEDN0391*206',   # DHA (FA) — 5 tissues, 脑脂质/抗炎
    'MEDP1414',       # Carnitine C12:1 (FA) — 4 tissues, 酰基肉碱/FA氧化
    'MEDN1430*207',   # Prostaglandin A2 (FA) — 3 tissues, Level 1a, 炎症介质
    'MEDP1328*282',   # LPC(20:2/0:0) (GP) — 3 tissues, 溶血磷脂
    'MEDN1414N',      # 12S-HHT (FA) — 3 tissues, Level 1a, 血栓素通路
    'ME0012956HN',    # PC(16:0/20:4) (GP) — 3 tissues, 磷脂周转
    'MEDP1693*265',   # LPC(12:0/0:0) (GP) — 3 tissues, 溶血磷脂
    'MEDN0356HN',     # LPS(22:6) (GP) — 3 tissues, Level 1a, 溶血磷脂
    'MEDN0370*212',   # LPA(0:0/18:0) (GP) — 3 tissues, Level 1a, 磷脂代谢
    'MEDN0398*204',   # AA (FA) — 3 tissues, Level 1a, 花生四烯酸
    'MEDN0106',       # Taurochenodesoxycholic Acid (Bile acids) — 2 tissues, Level 1a
    'MEDP1171P',      # Butenoyl-PAF (GP) — 2 tissues, 血小板活化因子
    'MADN3058',       # 5,6-DiHETrE (FA) — 2 tissues, Level 1b, EET衍生物
    'MEDP1399',       # Carnitine C16:2 (FA) — 2 tissues, 酰基肉碱
    'ME0054539N',     # 1-Eicosanoyl-sn-glycero-3-phosphocholine (GP) — 3 tissues, Level 1b
    'MEDP1902*270',   # LPC(16:2/0:0) (GP) — 2 tissues, 溶血磷脂
    'MEDN1277',       # LPE(0:0/20:5) (GP) — 2 tissues, 溶血磷脂
    'ME0057016HN',    # PC(18:1/16:0) (GP) — 2 tissues, 磷脂周转
]


def move_center(element_obj, scale, subfigure_size, calculate_center_param_dict=None, extra_offset=None):
    if calculate_center_param_dict is None:
        calculate_center_param_dict = {}
    if extra_offset is None:
        extra_offset = Vector(0., 0.)
    center = element_obj.calculate_center(element_obj, scale, **calculate_center_param_dict)
    center_bottom_offset = calculate_center_bottom_offset(center, subfigure_size) + extra_offset
    element_obj.move_and_scale(bottom_left_offset=center_bottom_offset)


def _place(figure_class, figure_data_parameter_dict, subfigure_bottom_left, subfigure_size):
    """probe + fit scale + 构造(**kwargs) + move_center, 各 Subfigure 复用。"""
    probe = figure_class(figure_data_parameter_dict=dict(figure_data_parameter_dict))
    scale = float(np.min(subfigure_size / probe.size))
    element = figure_class(**{
        ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
        ParameterName.scale: scale,
        ParameterName.bottom_left_offset: subfigure_bottom_left,
    })
    move_center(element, scale, subfigure_size)
    return element


_DATASET_CACHE = {}
def _get_ds():
    """共享数据集实例 (SubfigureC/E 共用, 避免重复 _load_dataset 打印两次黑名单过滤)。"""
    if 'ds' not in _DATASET_CACHE:
        _DATASET_CACHE['ds'] = _load_dataset(MetabolomicsExperimentName.data_20260412_metabolomics)
    return _DATASET_CACHE['ds']


def _class_enrichment(ds, diff_ids, background_ids=None, min_K=15):
    """化学类富集分析 (超几何 over-representation + BH-FDR)。

    diff_ids: 显著差异代谢物 ID; background_ids: 背景代谢物 (默认 metabolite_info 全集)。
    min_K: 背景代谢物数 K <= min_K 的类不参与 BH 校正 (小类随机性大、降低校正惩罚项 m),
           其 FDR 记为 1.0。
    返回 list[dict] {class, k, K, ratio, p, fdr}, 按 FDR 升序。
    ratio = (k/n) / (K/N) = 差异占比 / 背景占比, 去除注释密度偏差。
    """
    from scipy.stats import hypergeom
    from collections import Counter
    cls = ds.metabolite_info['Class I']
    diff_cls = cls.reindex(list(diff_ids)).dropna().astype(str)
    if background_ids is None:
        bg_cls = cls.dropna().astype(str)
    else:
        bg_cls = cls.reindex(list(background_ids)).dropna().astype(str)
    diff_count = Counter(diff_cls.tolist())
    bg_count = Counter(bg_cls.tolist())
    n = len(diff_cls)
    N = len(bg_cls)
    classes = sorted(set(diff_count) | set(bg_count))
    recs = []
    for c in classes:
        k = diff_count.get(c, 0)
        K = bg_count.get(c, 0)
        ratio = (k / n) / (K / N) if (K > 0 and n > 0) else float('nan')
        p = float(hypergeom.sf(k - 1, N, K, n)) if (K > 0 and k > 0) else 1.0
        recs.append({'class': c, 'k': k, 'K': K, 'ratio': ratio, 'p': p, 'fdr': 1.0})
    # BH-FDR 只对 K > min_K 的类 (小类不参与, 降低惩罚项 m)
    big = [r for r in recs if r['K'] > min_K]
    if big:
        ps = np.array([r['p'] for r in big])
        order = np.argsort(ps)
        m = len(ps)
        fdr_r = ps[order] * m / np.arange(1, m + 1)
        fdr_r = np.minimum.accumulate(fdr_r[::-1])[::-1]
        fdr = np.empty(m)
        fdr[order] = np.clip(fdr_r, 0, 1)
        for r, f in zip(big, fdr):
            r['fdr'] = float(f)
    recs.sort(key=lambda r: r['fdr'])
    return recs



# ============================================================
# Panel A: 多组织 OPLS-DA 网格
# ============================================================

class SubfigureA(Subfigure):
    subfigure_label = 'a'
    subfigure_title = 'multi_tissue_oplsda_grid'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('oplsda_all')
        if not cache.exists():
            print('Panel A: 未找到 OPLS-DA 缓存，请先运行 oplsda 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        cached = pickle_load(cache)
        results = cached['results']
        tissue_results = {r.tissue: r for r in results if hasattr(r, 'tissue')}

        tissue_order = [
            'cortex', 'hippocampus', 'hypothalamus',
            'liver', 'kidney', 'fat',
            'gastrocnemius', 'quadriceps', 'ileum',
            'serum',
        ]

        scale = 0.5
        element = Elements.MultiTissueOplsdaGridFigure(**{
            ParameterName.figure_data_parameter_dict: {
                'tissue_results': tissue_results,
                'tissue_order': tissue_order,
                'group_colors': GROUP_COLORS,
            },
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel B: 差异代谢物聚类热图
# ============================================================

def _build_clustered_heatmap_params(r):
    """从全量热图缓存构造 Panel B 主图参数。

    在 figure 层做行数限制（按 VIP 取 top-N），子集重聚类行，列聚类沿用缓存；
    三重映射 filtered_index↔z_matrix 行↔metabolite_labels 统一用行下标索引。
    返回 ClusteredHeatmapFigure 的 figure_data_parameter_dict，无数据返回 None。
    """
    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    filtered_index = list(r['filtered_index'])
    z_full = r['z_matrix']            # (n_full, n_cols)
    n_full = len(filtered_index)
    if n_full == 0 or z_full.size == 0:
        return None

    # —— 按 VIP 取 top-N ——
    vip_sub = r['vip_scores'].reindex(filtered_index).fillna(0)
    top_mids = (vip_sub.nlargest(HEATMAP_TOP_N).index.tolist()
                if n_full > HEATMAP_TOP_N else filtered_index)
    mid_to_local = {mid: i for i, mid in enumerate(filtered_index)}
    top_local_rows = [mid_to_local[mid] for mid in top_mids]

    z_sub = z_full[top_local_rows, :]                       # (N, n_cols)
    sub_index = [filtered_index[i] for i in top_local_rows]

    # —— 子集重聚类行；列聚类沿用缓存 ——
    d = np.nan_to_num(z_sub, nan=0.0)
    if len(top_local_rows) > 1:
        row_linkage_sub = hierarchy.linkage(distance.pdist(d, metric='euclidean'), method='ward')
        row_order_sub = [int(x) for x in hierarchy.leaves_list(row_linkage_sub)]
    else:
        row_order_sub = [0]
    col_order = list(r['col_order'])
    col_linkage = r['col_linkage']

    # —— 排序到「行×列聚类」顺序 ——
    z_ordered = z_sub[np.ix_(row_order_sub, col_order)]

    groups = r['groups']
    groups_ordered = [groups[c] for c in col_order]

    met_labels = r['metabolite_labels']
    met_labels_sub = [met_labels[i] for i in top_local_rows]
    met_labels_ordered = [met_labels_sub[r] for r in row_order_sub]

    class_i_map = r.get('class_i', {})
    class_i_sub = [class_i_map.get(mid, 'Other') for mid in sub_index]
    class_i_ordered = [class_i_sub[r] for r in row_order_sub]

    sample_labels = r['sample_labels']
    sample_labels_ordered = [sample_labels[c] for c in col_order]

    present_classes = sorted(set(class_i_ordered))
    class_color_map = {c: lipid_class_color(c) for c in present_classes}

    # —— 右侧强调代谢物 → z_ordered 行下标 ——
    name_to_index = r.get('metabolite_name_to_index', {})
    sub_index_pos = {mid: i for i, mid in enumerate(sub_index)}   # mid → z_sub 行位置
    highlight_indices, highlight_labels = [], []
    for h in HEATMAP_HIGHLIGHT_METABOLITES:
        mid = h if h in sub_index_pos else name_to_index.get(h)
        if mid is None or mid not in sub_index_pos:
            continue
        pos = row_order_sub.index(sub_index_pos[mid])            # z_ordered 行下标
        highlight_indices.append(pos)
        highlight_labels.append(met_labels_ordered[pos])

    # —— z-score 色标范围（自适应，截断到 [1.5, 3]）——
    finite = z_ordered[np.isfinite(z_ordered)]
    z_abs_max = float(np.max(np.abs(finite))) if finite.size else 2.0
    z_bound = min(max(z_abs_max, 1.5), 3.0)
    z_lim = (-z_bound, z_bound)

    return {
        'z_matrix_ordered': z_ordered,
        'col_linkage': col_linkage,
        'groups_ordered': groups_ordered,
        'group_colors': GROUP_COLORS,
        'metabolite_labels': met_labels_ordered,
        'sample_labels': sample_labels_ordered,
        'class_i_ordered': class_i_ordered,
        'class_color_map': class_color_map,
        'highlight_indices': highlight_indices,
        'highlight_labels': highlight_labels,
        'z_lim': z_lim,
        'cmap': gray_teal_diverging_cmap(),
    }


class SubfigureBOld(Subfigure):
    subfigure_label = 'b'
    subfigure_title = 'clustered_heatmap'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('heatmap', tissue='hippocampus')
        if not cache.exists():
            print('Panel B: 未找到聚类热图缓存，请先运行 clustered_heatmap 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        params = _build_clustered_heatmap_params(pickle_load(cache))
        if params is None:
            print('Panel B: 热图缓存无有效数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        probe = Elements.ClusteredHeatmapFigure(figure_data_parameter_dict=dict(params))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.ClusteredHeatmapFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel BNew: 各器官差异代谢物计数 + 脂质占比 (设计方案 B)
# ============================================================

class SubfigureB(Subfigure):
    subfigure_label = 'b'
    subfigure_title = 'class_enrichment_lollipop'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        ds = _get_ds()
        tissues = ['cortex', 'hippocampus', 'hypothalamus', 'liver', 'kidney',
                   'fat', 'gastrocnemius', 'quadriceps', 'ileum', 'serum']
        fc_thr = np.log2(1.5)
        p_thr = 0.05
        scale = 0.37
        diff_set = set()
        for t in tissues:
            cache = _cache_path('diff', tissue=t)
            if not cache.exists():
                continue
            diff = pickle_load(cache)
            if diff is None or diff.empty:
                continue
            sig = diff[(np.abs(diff['log2FC']) >= fc_thr) & (diff['p_value'] < p_thr)
                       & diff['p_value'].notna()]
            diff_set |= set(sig.index)
        if not diff_set:
            print('Panel B: 未找到差异分析缓存，请先运行 volcano 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 全局化学类富集分析 (超几何 + BH-FDR), 按 FDR 排序
        recs_all = _class_enrichment(ds, diff_set)
        # 图中只保留 K>=15 且 FDR<0.99 的类 (去掉小类 + 完全不显著的, 即 SP 及之前)
        recs = [r for r in recs_all if r['K'] >= 15 and r['fdr'] < 0.9]
        # all_pct (K/N*100) 决定棒棒糖点大小 = 类规模
        from collections import Counter
        cls_all = ds.metabolite_info['Class I'].dropna().astype(str)
        N_all = len(cls_all)
        bg_count = Counter(cls_all.tolist())
        all_pct = {c: bg_count.get(c, 0) / N_all * 100 for c in {r['class'] for r in recs}}
        # 打印图中各类的 p / FDR 完整列表
        from ...common.common_figure_materials import fdr_star
        print('=== Panel B 图中各类 (K>=15) p / FDR ===')
        for r in recs:
            print(f"  {r['class']:42s} k/K={r['k']:>3}/{r['K']:<4} ratio={r['ratio']:.2f}  "
                  f"p={r['p']:.4f}  FDR={r['fdr']:.4f}  {fdr_star(r['fdr'])}")

        params = {'recs': recs, 'all_pct': all_pct}
        element = Elements.EnrichmentLollipopFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        # 导出补充材料: Panel B 全局化学类富集表 (全组织合并 diff_set 的超几何 + BH-FDR)
        from pathlib import Path as _P
        _src = _P(__file__).parents[3] / 'output' / 'metabolomics' / 'figure_source'
        _src.mkdir(parents=True, exist_ok=True)
        _b = pd.DataFrame(recs_all)[['class', 'k', 'K', 'ratio', 'p', 'fdr']]
        _b['in_figure'] = ((_b['K'] >= 15) & (_b['fdr'] < 0.9))  # 图中画出的类
        _b.sort_values('fdr').to_excel(_src / 'panel_b_global_class_enrichment.xlsx', index=False)
        print(f'已导出 Panel B 富集表: {_src / "panel_b_global_class_enrichment.xlsx"}')
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel C: 火山图
# ============================================================

class SubfigureCOld(Subfigure):
    subfigure_label = 'c'
    subfigure_title = 'volcano_plot'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('volcano', tissue='hippocampus')
        if not cache.exists():
            print('Panel C: 未找到火山图缓存，请先运行 volcano 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        vol_data = pickle_load(cache)
        t = 'hippocampus'
        if t not in vol_data:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        d = vol_data[t]

        # 取 Class I（决定点颜色/大小）+ 重要脂质下标（与 Panel B 同列表）
        ds = _get_ds()
        if 'Class I' in ds.metabolite_info.columns:
            class_series = ds.metabolite_info['Class I'].reindex(d['index'])
            classes = [class_series.iloc[i] for i in range(len(d['index']))]
        else:
            classes = ['Other'] * len(d['index'])

        index_to_pos = {idx: i for i, idx in enumerate(d['index'])}
        highlight_indices = [index_to_pos[h] for h in HEATMAP_HIGHLIGHT_METABOLITES
                             if h in index_to_pos]

        figure_data_parameter_dict = {
            'log2fc': d['log2fc'],
            'neg_log10p': d['neg_log10p'],
            'classes': classes,
            'names': d['names'],
            'highlight_indices': highlight_indices,
            'thresholds': d['thresholds'],
        }
        probe = Elements.MetabolomicsVolcanoFigure(figure_data_parameter_dict=dict(figure_data_parameter_dict))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.MetabolomicsVolcanoFigure(**{
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel D: 跨器官通路富集点图
# ============================================================

# ============================================================
# Panel CNew: 器官 × 化学类 富集热图 (设计方案 v2, HERO)
# ============================================================

class SubfigureC(Subfigure):
    subfigure_label = 'c'
    subfigure_title = 'organ_class_enrichment_heatmap'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        ds = _get_ds()
        organs = ['cortex', 'hippocampus', 'hypothalamus', 'liver', 'kidney',
                  'fat', 'gastrocnemius', 'quadriceps', 'ileum', 'serum']
        fc_thr = np.log2(1.5)
        p_thr = 0.05
        lipid_classes = ['GP', 'GL', 'SP', 'FA']
        nonlipid_classes = [
            'Amino acid and Its metabolites', 'Organic acid and Its derivatives',
            'Carbohydrates and Its metabolites', 'Nucleotide and Its metabolites']
        classes = lipid_classes + nonlipid_classes
        scale = 0.3

        ratio_matrix, fdr_matrix, valid_organs, organ_recs = [], [], [], []
        for t in organs:
            cache = _cache_path('diff', tissue=t)
            if not cache.exists():
                continue
            diff = pickle_load(cache)
            if diff is None or diff.empty:
                continue
            sig = diff[(np.abs(diff['log2FC']) >= fc_thr) & (diff['p_value'] < p_thr)
                       & diff['p_value'].notna()]
            # 背景 = 该组织测到的代谢物 (非空行)
            tissue_data = ds.get_data_by_tissue(t)
            bg_ids = tissue_data.index[tissue_data.notna().any(axis=1)]
            recs = _class_enrichment(ds, set(sig.index), background_ids=bg_ids)
            rec_by_class = {r['class']: r for r in recs}
            ratio_matrix.append([rec_by_class.get(c, {}).get('ratio', 1.0) for c in classes])
            fdr_matrix.append([rec_by_class.get(c, {}).get('fdr', 1.0) for c in classes])
            valid_organs.append(t)
            organ_recs.append(rec_by_class)

        if not valid_organs:
            print('Panel C: 未找到差异分析缓存，请先运行 volcano 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        params = {
            'organs': valid_organs, 'classes': classes,
            'ratio_matrix': np.array(ratio_matrix, dtype=float),
            'fdr_matrix': np.array(fdr_matrix, dtype=float),
            'lipid_n': len(lipid_classes),
        }
        element = Elements.OrganClassEnrichmentHeatmap(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        # 导出补充材料: Panel C 器官×类富集长表
        from pathlib import Path as _P
        _src = _P(__file__).parents[3] / 'output' / 'metabolomics' / 'figure_source'
        _src.mkdir(parents=True, exist_ok=True)
        _rows = []
        for _t, _rbc in zip(valid_organs, organ_recs):
            for _c in classes:
                _r = _rbc.get(_c, {})
                _rows.append({'organ': _t, 'class': _c, 'k': _r.get('k', 0), 'K': _r.get('K', 0),
                              'ratio': _r.get('ratio', 1.0), 'fdr': _r.get('fdr', 1.0)})
        pd.DataFrame(_rows).sort_values(['organ', 'fdr']).to_excel(
            _src / 'panel_c_organ_class_enrichment.xlsx', index=False)
        print(f'已导出 Panel C 富集表: {_src / "panel_c_organ_class_enrichment.xlsx"}')
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


class SubfigureDOld(Subfigure):
    subfigure_label = 'd'
    subfigure_title = 'pathway_bubble'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('pathway_bubble')
        if not cache.exists():
            print('Panel D: 未找到通路气泡图缓存，请先运行 pathway_bubble 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        pathway_df = pickle_load(cache)

        tissue_order = [
            'cortex', 'hippocampus', 'hypothalamus',
            'liver', 'kidney', 'fat',
            'gastrocnemius', 'quadriceps', 'ileum',
            'serum',
        ]

        figure_data_parameter_dict = {
            'pathway_df': pathway_df,
            'max_pathways': 20,
            'color_by': 'direction',
            'tissue_order': tissue_order,
        }
        probe = Elements.PathwayBubbleFigure(figure_data_parameter_dict=dict(figure_data_parameter_dict))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.PathwayBubbleFigure(**{
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel DNew: 已注释脂质代谢物跨器官方向热图 (设计方案 v3 D1)
# ============================================================

class SubfigureD(Subfigure):
    subfigure_label = 'd'
    subfigure_title = 'lipid_direction_heatmap'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        if not LIPID_METABOLITES_FOR_DEF:
            print('Panel D: LIPID_METABOLITES_FOR_DEF 为空，请运行 select_lipid_metabolites_for_DEF.py 后填入')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        scale = 0.37
        ds = _get_ds()
        tissues = [
            'cortex', 'hippocampus', 'hypothalamus',
            'liver', 'kidney', 'fat',
            'gastrocnemius', 'quadriceps', 'ileum',
            'serum',
        ]

        # 从各组织 diff cache 提取 LIPID_METABOLITES_FOR_DEF 的 log2FC
        met_ids = LIPID_METABOLITES_FOR_DEF
        direction_rows = []
        valid_mets = []
        for mid in met_ids:
            row = []
            for t in tissues:
                cache = _cache_path('diff', tissue=t)
                if not cache.exists():
                    row.append(np.nan)
                    continue
                diff = pickle_load(cache)
                if diff is None or diff.empty or mid not in diff.index:
                    row.append(np.nan)
                else:
                    val = diff.loc[mid, 'log2FC']
                    row.append(float(val) if pd.notna(val) else np.nan)
            if any(np.isfinite(v) for v in row):
                direction_rows.append(row)
                valid_mets.append(mid)

        if not direction_rows:
            print('Panel D: 指定代谢物在各组织均无有效 log2FC')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        direction_matrix = np.array(direction_rows, dtype=float)

        name_map = ds.metabolite_info['Compounds'].to_dict() if 'Compounds' in ds.metabolite_info.columns else {}
        met_labels = [str(name_map.get(m, m)) for m in valid_mets]
        tissue_short = {
            'cortex': 'Ctx', 'hippocampus': 'Hipp', 'hypothalamus': 'Hyp',
            'liver': 'Liv', 'kidney': 'Kid', 'fat': 'Fat',
            'gastrocnemius': 'Gas', 'quadriceps': 'Quad', 'ileum': 'Ile',
            'serum': 'Ser',
        }
        tissue_labels = [tissue_short.get(t, t[:4]) for t in tissues]

        params = {
            'direction_matrix': direction_matrix,
            'metabolite_labels': met_labels,
            'tissue_labels': tissue_labels,
        }

        element = Elements.LipidMetaboliteDirectionFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel E: 关键代谢物跨组织柱状图 (旧版)
# ============================================================

class SubfigureEOld(Subfigure):
    subfigure_label = 'e'
    subfigure_title = 'cross_tissue_bar'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        metabolite_ids = CROSS_TISSUE_BAR_METABOLITES
        if not metabolite_ids:
            print('Panel E: CROSS_TISSUE_BAR_METABOLITES 为空，请按故事线填入代谢物 ID')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 渲染时直接读取数据集内联生成（特定代谢物比较图，无需预计算缓存）
        ds = _get_ds()
        group1, group2 = Keyword.ad, Keyword.sema
        tissues = ds.list_tissues()
        name_map = ds.metabolite_info['Compounds'].to_dict() if 'Compounds' in ds.metabolite_info.columns else {}

        values, metabolite_names = {}, {}
        for met_id in metabolite_ids:
            met_values = {}
            for t in tissues:
                g1_data = ds.get_data_by_group(group1, tissue=t)
                g2_data = ds.get_data_by_group(group2, tissue=t)
                if met_id not in g1_data.index or met_id not in g2_data.index:
                    continue
                g1_vals = pd.to_numeric(g1_data.loc[met_id], errors='coerce').dropna().values
                g2_vals = pd.to_numeric(g2_data.loc[met_id], errors='coerce').dropna().values
                if len(g1_vals) == 0 or len(g2_vals) == 0:
                    continue
                met_values[t] = {group1: g1_vals, group2: g2_vals}
            if met_values:
                values[met_id] = met_values
                metabolite_names[met_id] = name_map.get(met_id, str(met_id))

        if not values:
            print('Panel E: 指定代谢物在数据集中均无有效数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        figure_data_parameter_dict = {
            'metabolites': list(values.keys()),
            'values': values,
            'metabolite_names': metabolite_names,
            'tissue_order': tissues,
            'group_colors': GROUP_COLORS,
            'group_names': (group1, group2),
        }
        probe = Elements.CrossTissueBarFigure(figure_data_parameter_dict=dict(figure_data_parameter_dict))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.CrossTissueBarFigure(**{
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel ENew: 海马差异代谢物聚类热图 (设计方案 v3, 高亮脂质代谢物)
# ============================================================

def _build_panel_e_heatmap_params(r):
    """从海马热图缓存构造 Panel E 的 ClusteredHeatmapFigure 参数。

    行限制 top-N (50)，子集重聚类行，列聚类沿用缓存。
    高亮 LIPID_METABOLITES_FOR_DEF 中的代谢物。
    """
    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    filtered_index = list(r['filtered_index'])
    z_full = r['z_matrix']
    n_full = len(filtered_index)
    if n_full == 0 or z_full.size == 0:
        return None

    PANEL_E_TOP_N = 50
    vip_sub = r['vip_scores'].reindex(filtered_index).fillna(0)
    top_mids = (vip_sub.nlargest(PANEL_E_TOP_N).index.tolist()
                if n_full > PANEL_E_TOP_N else filtered_index)
    mid_to_local = {mid: i for i, mid in enumerate(filtered_index)}
    top_local_rows = [mid_to_local[mid] for mid in top_mids]

    z_sub = z_full[top_local_rows, :]
    sub_index = [filtered_index[i] for i in top_local_rows]

    d = np.nan_to_num(z_sub, nan=0.0)
    if len(top_local_rows) > 1:
        row_linkage_sub = hierarchy.linkage(distance.pdist(d, metric='euclidean'), method='ward')
        row_order_sub = [int(x) for x in hierarchy.leaves_list(row_linkage_sub)]
    else:
        row_order_sub = [0]
    col_order = list(r['col_order'])
    col_linkage = r['col_linkage']

    z_ordered = z_sub[np.ix_(row_order_sub, col_order)]
    groups_ordered = [r['groups'][c] for c in col_order]

    met_labels = r['metabolite_labels']
    met_labels_sub = [met_labels[i] for i in top_local_rows]
    met_labels_ordered = [met_labels_sub[r] for r in row_order_sub]

    class_i_map = r.get('class_i', {})
    class_i_sub = [class_i_map.get(mid, 'Other') for mid in sub_index]
    class_i_ordered = [class_i_sub[r] for r in row_order_sub]

    sample_labels = r['sample_labels']
    sample_labels_ordered = [sample_labels[c] for c in col_order]

    present_classes = sorted(set(class_i_ordered))
    class_color_map = {c: lipid_class_color(c) for c in present_classes}

    # 高亮 LIPID_METABOLITES_FOR_DEF
    sub_index_pos = {mid: i for i, mid in enumerate(sub_index)}
    highlight_indices, highlight_labels = [], []
    for h in LIPID_METABOLITES_FOR_DEF:
        if h not in sub_index_pos:
            continue
        pos = row_order_sub.index(sub_index_pos[h])
        highlight_indices.append(pos)
        highlight_labels.append(met_labels_ordered[pos])

    finite = z_ordered[np.isfinite(z_ordered)]
    z_abs_max = float(np.max(np.abs(finite))) if finite.size else 2.0
    z_bound = min(max(z_abs_max, 1.5), 3.0)
    z_lim = (-z_bound, z_bound)

    return {
        'z_matrix_ordered': z_ordered,
        'col_linkage': col_linkage,
        'groups_ordered': groups_ordered,
        'group_colors': GROUP_COLORS,
        'metabolite_labels': met_labels_ordered,
        'sample_labels': sample_labels_ordered,
        'class_i_ordered': class_i_ordered,
        'class_color_map': class_color_map,
        'highlight_indices': highlight_indices,
        'highlight_labels': highlight_labels,
        'z_lim': z_lim,
        'cmap': rdbu_diverging_cmap(),
    }


class SubfigureE(Subfigure):
    subfigure_label = 'e'
    subfigure_title = 'hippocampus_clustered_heatmap'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('heatmap', tissue='hippocampus')
        if not cache.exists():
            print('Panel E: 未找到聚类热图缓存，请先运行 clustered_heatmap 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        params = _build_panel_e_heatmap_params(pickle_load(cache))
        if params is None:
            print('Panel E: 海马热图缓存无有效数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        element = Elements.ClusteredHeatmapFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: 0.26,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, 0.26, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel F: S-plot (旧版)

class SubfigureFOld(Subfigure):
    subfigure_label = 'f'
    subfigure_title = 'splot'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('splot', tissue='hippocampus')
        if not cache.exists():
            print('Panel F: 未找到 S-plot 缓存，请先运行 splot 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        cached = pickle_load(cache)
        result = cached['oplsda_result']
        candidates = cached['candidates']

        figure_data_parameter_dict = {
            'p_loadings': result.p_loadings.values,
            'p_corr': result.p_corr.values,
            'vip_scores': result.vip_scores.values,
            'candidates_mask': candidates['candidate'].values if 'candidate' in candidates.columns else None,
            'metabolite_names': (
                candidates['Compounds'].tolist()
                if 'Compounds' in candidates.columns
                else [str(i) for i in candidates.index]
            ),
            'p_corr_threshold': 0.5,
            'tissue': 'hippocampus',
            'r2y': result.r2y,
            'q2': result.q2,
        }
        probe = Elements.SplotFigure(figure_data_parameter_dict=dict(figure_data_parameter_dict))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.SplotFigure(**{
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel F: 关键脂质代谢物跨器官柱图 (设计方案 v3)
# ============================================================

class SubfigureFOld2(Subfigure):
    subfigure_label = 'f'
    subfigure_title = 'lipid_cross_tissue_bar'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        met_ids = LIPID_METABOLITES_FOR_DEF
        if not met_ids:
            print('Panel F: LIPID_METABOLITES_FOR_DEF 为空')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        scale = 0.22
        ds = _get_ds()
        group1, group2 = Keyword.ad, Keyword.sema
        tissues = ds.list_tissues()
        name_map = ds.metabolite_info['Compounds'].to_dict() if 'Compounds' in ds.metabolite_info.columns else {}

        values, metabolite_names = {}, {}
        for met_id in met_ids:
            met_values = {}
            for t in tissues:
                g1_data = ds.get_data_by_group(group1, tissue=t)
                g2_data = ds.get_data_by_group(group2, tissue=t)
                if met_id not in g1_data.index or met_id not in g2_data.index:
                    continue
                g1_vals = pd.to_numeric(g1_data.loc[met_id], errors='coerce').dropna().values
                g2_vals = pd.to_numeric(g2_data.loc[met_id], errors='coerce').dropna().values
                if len(g1_vals) == 0 or len(g2_vals) == 0:
                    continue
                met_values[t] = {group1: g1_vals, group2: g2_vals}
            if met_values:
                values[met_id] = met_values
                metabolite_names[met_id] = name_map.get(met_id, str(met_id))

        if not values:
            print('Panel F: 指定代谢物在数据集中均无有效数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        params = {
            'metabolites': list(values.keys()),
            'values': values,
            'metabolite_names': metabolite_names,
            'tissue_order': tissues,
            'group_colors': GROUP_COLORS,
            'group_names': (group1, group2),
        }

        element = Elements.CrossTissueBarFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel F2: 收窄柱图 — 6 焊点分子 (= G 主角集), 限定脑(海马/皮层)+血浆
# 设计方案"改法二": F 分子集 = G 分子集, 形成 E→F→G 机制单元
# ============================================================

F2_FOCUS_METABOLITES = [
    'ME0054539N',    # 1-Eicosanoyl-GPC — 醚磷脂信号
    'MEDN0106',      # Taurochenodesoxycholic Acid — 胆汁酸
    'MEDN0370*212',  # LPA(0:0/18:0) — 溶血磷脂酸
    'MEDN1277',      # LPE(0:0/20:5) — 溶血磷脂乙醇胺
    'MEDP1171P',     # Butenoyl-PAF — 血小板活化因子
    'MEDP1693*265',  # LPC(12:0/0:0) — 溶血磷脂胆碱
]
F2_TISSUES = ['hippocampus', 'cortex', 'serum']


class SubfigureF2(Subfigure):
    """Panel F 改法二: 收窄到 6 个焊点分子, 限定脑+血浆, 标方向/星号。
    分子集 = G 的主角, 与 G 形成机制单元。复用 CrossTissueBarFigure。"""
    subfigure_label = 'f'
    subfigure_title = 'lipid_focus_cross_tissue_bar'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        met_ids = F2_FOCUS_METABOLITES
        if not met_ids:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        scale = 0.22
        ds = _get_ds()
        group1, group2 = Keyword.ad, Keyword.sema
        tissues = [t for t in F2_TISSUES if t in ds.list_tissues()]
        name_map = (ds.metabolite_info['Compounds'].to_dict()
                    if 'Compounds' in ds.metabolite_info.columns else {})
        values, metabolite_names = {}, {}
        for met_id in met_ids:
            met_values = {}
            for t in tissues:
                g1_data = ds.get_data_by_group(group1, tissue=t)
                g2_data = ds.get_data_by_group(group2, tissue=t)
                if met_id not in g1_data.index or met_id not in g2_data.index:
                    continue
                g1_vals = pd.to_numeric(g1_data.loc[met_id], errors='coerce').dropna().values
                g2_vals = pd.to_numeric(g2_data.loc[met_id], errors='coerce').dropna().values
                if len(g1_vals) == 0 or len(g2_vals) == 0:
                    continue
                met_values[t] = {group1: g1_vals, group2: g2_vals}
            if met_values:
                values[met_id] = met_values
                metabolite_names[met_id] = name_map.get(met_id, str(met_id))
        if not values:
            print('Panel F2: 指定代谢物无有效数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        params = {
            'metabolites': list(values.keys()),
            'values': values,
            'metabolite_names': metabolite_names,
            'tissue_order': tissues,
            'group_colors': GROUP_COLORS,
            'group_names': (group1, group2),
        }
        element = Elements.CrossTissueBarFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel F1: 磷脂重塑通路示意图 (设计方案"改法一")
# PC → PLA2 → 溶血磷脂(LPC/LPE/LPS/LPA) + PUFA(AA/DHA) → COX/LOX → 氧化脂质
# 节点按海马 log2FC 方向上色 (红=上调, 蓝=下调, 灰=酶/无数据)
# ============================================================

# (metabolite_id, label, x, y, is_enzyme); x/y 为归一化坐标
PATHWAY_NODE_DEFS = [
    ('ME0012956HN', 'PC 16:0/20:4', 0.08, 0.86, False),
    ('ME0057016HN', 'PC 18:1/16:0', 0.08, 0.72, False),
    (None, 'PLA2', 0.26, 0.79, True),
    ('MEDP1328*282', 'LPC 20:2', 0.48, 0.92, False),
    ('MEDP1693*265', 'LPC 12:0', 0.48, 0.80, False),
    ('MEDP1902*270', 'LPC 16:2', 0.48, 0.68, False),
    ('MEDN1277', 'LPE 20:5', 0.48, 0.56, False),
    ('MEDN0356HN', 'LPS 22:6', 0.48, 0.44, False),
    ('MEDN0370*212', 'LPA 18:0', 0.48, 0.32, False),
    ('MEDN0398*204', 'AA', 0.48, 0.20, False),
    ('MEDN0391*206', 'DHA', 0.48, 0.10, False),
    (None, 'COX/LOX', 0.68, 0.15, True),
    ('MEDN1430*207', 'PGA2', 0.90, 0.24, False),
    ('MEDN1414N', '12S-HHT', 0.90, 0.12, False),
    ('MADN3058', '5,6-DiHETrE', 0.90, 0.06, False),
    ('ME0054539N', '1-Eic-GPC', 0.68, 0.50, False),
    ('MEDP1171P', 'PAF', 0.68, 0.62, False),
]
PATHWAY_EDGES = [
    (0, 2), (1, 2),                              # PC → PLA2
    (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),  # PLA2 → lysoPL
    (2, 9), (2, 10),                             # PLA2 → PUFA
    (9, 11),                                     # AA → COX/LOX
    (11, 12), (11, 13), (11, 14),               # COX/LOX → oxylipin
]


class SubfigureF1(Subfigure):
    """Panel F 改法一: 磷脂重塑通路示意图, 节点按海马 log2FC 方向上色。
    描述性着色 (非富集检验, 无偏倚), 把 PC→PLA2→溶血磷脂+PUFA→氧化脂质 画成一句话。"""
    subfigure_label = 'f'
    subfigure_title = 'lipid_remodeling_pathway'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('diff', tissue='hippocampus')
        if not cache.exists():
            print('Panel F1: 未找到 hippocampus diff 缓存')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        diff = pickle_load(cache)
        log2fc_map = diff['log2FC'].to_dict()
        nodes = []
        for met_id, label, x, y, is_enzyme in PATHWAY_NODE_DEFS:
            fc = log2fc_map.get(met_id) if met_id else None
            nodes.append({'label': label, 'x': x, 'y': y, 'log2fc': fc, 'is_enzyme': is_enzyme})
        element = _place(
            Elements.PathwayDiagramFigure,
            {'nodes': nodes, 'edges': PATHWAY_EDGES},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel G: 代谢物 × 表型 相关散点网格 (散点+回归线, 按组上色)
# ============================================================

def _build_corr_panels(ds, corr_df, pheno_df, tissue, selections, group_colors):
    """从 phenotype_corr 缓存 + ds 构造散点 panels (CorrelationScatterGridFigure 格式)。
    selections: list[(metabolite_id, phenotype)]。"""
    tissue_data = ds.get_data_by_tissue(tissue)
    sample_animals = {col: ds.sample_info.get(col, {}).get('animal_id') for col in tissue_data.columns}
    name_col = 'metabolite_name' if 'metabolite_name' in corr_df.columns else 'metabolite_id'
    panels = []
    for met_id, phenotype in selections:
        x_list, y_list, colors = [], [], []
        for col in tissue_data.columns:
            animal = sample_animals.get(col)
            if animal is None or animal not in pheno_df.index or met_id not in tissue_data.index:
                continue
            abundance = pd.to_numeric(tissue_data.loc[met_id, col], errors='coerce')
            if pd.isna(abundance):
                continue
            pheno_val = pheno_df.loc[animal, phenotype]
            if pd.isna(pheno_val):
                continue
            x_list.append(float(abundance))
            y_list.append(float(pheno_val))
            group = ds._group_map.get(animal)
            colors.append(group_colors.get(group, CCD.gray))
        if len(x_list) < 3:
            continue
        x_arr, y_arr = np.array(x_list), np.array(y_list)
        z = np.polyfit(x_arr, y_arr, 1)
        fit_x = np.linspace(x_arr.min(), x_arr.max(), 20)
        fit_y = np.polyval(z, fit_x)
        corr_row = corr_df[(corr_df['metabolite_id'] == met_id) & (corr_df['phenotype'] == phenotype)]
        r_val = float(corr_row['r'].iloc[0]) if len(corr_row) else 0.0
        p_val = float(corr_row['p_value'].iloc[0]) if len(corr_row) else 1.0
        met_name = (str(corr_row[name_col].iloc[0]) if len(corr_row)
                    else str(ds.metabolite_info['Compounds'].get(met_id, met_id)))
        pheno_en = PHENO_NAME_EN.get(phenotype, phenotype)
        panels.append({
            'x': x_arr, 'y': y_arr, 'group_colors': colors,
            'fit_x': fit_x, 'fit_y': fit_y,
            'title': f'{met_name[:25]}\nr={r_val:+.2f} p={p_val:.3f}',
            'x_label': 'Abundance', 'y_label': pheno_en,
        })
    return panels


def _top_corr_selections(corr_df, rank=1, top_n_phenos=4):
    """每表型取第 rank 的代谢物 (by |r|), 按 |r| 排序取 top_n_phenos 个表型。
    返回 list[(metabolite_id, phenotype)]。"""
    top_per_pheno = []
    for ph in sorted(corr_df['phenotype'].unique()):
        sub = corr_df[corr_df['phenotype'] == ph].copy()
        sub['abs_r'] = sub['r'].abs()
        sorted_sub = sub.sort_values('abs_r', ascending=False)
        if len(sorted_sub) >= rank:
            row = sorted_sub.iloc[rank - 1]
            top_per_pheno.append((row['metabolite_id'], ph, row['abs_r']))
    top_per_pheno.sort(key=lambda x: -x[2])
    return [(mid, ph) for mid, ph, _ in top_per_pheno[:top_n_phenos]]


# Panel G 认知列 (热图主体): [水迷宫 learn-rate | NOR DI 鼻尖 | 合成认知评分 CCS]
PANEL_G_PHENOTYPES = ['MWM_learn_rate', 'NOR_DI_鼻尖', 'CCS']
PANEL_G_TOP_ROWS = 8


class SubfigureG(Subfigure):
    """Panel G: 海马脂质 × 认知表型相关热图 + 侧边头条×CCS 散点。
    列 [MWM learn-rate | NOR DI | CCS], Spearman r (rdbu 发散色) + raw p 星号; 合并 treated+veh。
    行 = 脂质代谢物中至少一认知列 raw p<0.05 者, 按 mean|r| 取 top N。
    n=12 pilot: BH-FDR 后多不显著, 星号用 nominal p, 图注声明。"""
    subfigure_label = 'g'
    subfigure_title = 'cognition_lipid_correlation'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('phenotype_corr')
        if not cache.exists():
            print('Panel G: 未找到表型相关缓存，请先运行 phenotype_corr 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        cached = pickle_load(cache)
        corr = cached.get('correlations', {}).get('hippocampus')
        pheno = cached.get('pheno')
        if corr is None or corr.empty or pheno is None:
            print('Panel G: 表型相关数据为空')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        ds = _get_ds()

        _sig_p = corr[corr['p_value'] < 0.05]
        _sig_adj = corr[corr['adj_p'] < 0.1]
        print(f'=== Panel G: hippocampus 认知相关 (raw p<0.05: {len(_sig_p)} 对/'
              f'{_sig_p["metabolite_id"].nunique()} 代谢物; adj_p<0.1: {len(_sig_adj)} 对 — '
              f'pilot 下 BH 多不显著, 星号用 nominal p) ===')

        pheno_cols = [p for p in PANEL_G_PHENOTYPES if p in corr['phenotype'].unique()]
        plot_df = corr[corr['phenotype'].isin(pheno_cols)].copy()
        if plot_df.empty:
            print('Panel G: 无认知列表型数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 仅脂质代谢物
        info = ds.metabolite_info
        if 'Class I' in info.columns:
            lipid_ids = {str(i) for i, c in info['Class I'].items() if is_lipid_class(c)}
            plot_df = plot_df[plot_df['metabolite_id'].isin(lipid_ids)]

        # 至少一认知列 raw p<0.05; 按 mean|r| 取 top N
        sig_ids = plot_df[plot_df['p_value'] < 0.05]['metabolite_id'].unique()
        if len(sig_ids) > 0:
            plot_df = plot_df[plot_df['metabolite_id'].isin(sig_ids)]
        mean_abs_r = plot_df.groupby('metabolite_id')['r'].apply(lambda s: s.abs().mean())
        top_ids = mean_abs_r.sort_values(ascending=False).head(PANEL_G_TOP_ROWS).index.tolist()
        plot_df = plot_df[plot_df['metabolite_id'].isin(top_ids)]
        if plot_df.empty:
            print('Panel G: 脂质代谢物无认知相关数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # r / p 矩阵 (行=top_ids, 列=pheno_cols)
        name_col = 'metabolite_name' if 'metabolite_name' in plot_df.columns else 'metabolite_id'
        name_map = (plot_df.drop_duplicates('metabolite_id')
                    .set_index('metabolite_id')[name_col].to_dict())
        r_pivot = plot_df.pivot(index='metabolite_id', columns='phenotype', values='r')
        p_pivot = plot_df.pivot(index='metabolite_id', columns='phenotype', values='p_value')
        r_pivot = r_pivot.reindex(index=top_ids, columns=pheno_cols).fillna(0)
        p_pivot = p_pivot.reindex(index=top_ids, columns=pheno_cols).fillna(1)
        met_labels = [str(name_map.get(m, m)) for m in top_ids]
        pheno_labels = [PHENO_NAME_EN.get(p, p) for p in pheno_cols]

        for mid in top_ids:
            row_r = r_pivot.loc[mid]
            best = row_r.abs().idxmax()
            print(f'  {str(mid):16s} {str(name_map.get(mid, mid))[:34]:34s} '
                  f'mean|r|={row_r.abs().mean():.2f} best={best} r={row_r[best]:+.2f}')

        # 侧边散点: CCS |r| 最大 (优先 raw p<0.05) 的代谢物
        inset_panel = None
        if 'CCS' in pheno_cols:
            ccs_df = plot_df[plot_df['phenotype'] == 'CCS']
            ccs_sig = ccs_df[ccs_df['p_value'] < 0.05]
            pick_df = ccs_sig if not ccs_sig.empty else ccs_df
            if not pick_df.empty:
                ccs_top = pick_df.reindex(pick_df['r'].abs().sort_values(ascending=False).index).iloc[0]
                panels = _build_corr_panels(ds, corr, pheno, 'hippocampus',
                                            [(ccs_top['metabolite_id'], 'CCS')], GROUP_COLORS)
                if panels:
                    inset_panel = panels[0]

        figure_data_parameter_dict = {
            'r_matrix': r_pivot.values,
            'p_matrix': p_pivot.values,
            'metabolite_labels': met_labels,
            'phenotype_labels': pheno_labels,
            'inset_panel': inset_panel,
            'tissue': 'hippocampus',
        }
        element = _place(Elements.CognitionCorrelationFigure, figure_data_parameter_dict,
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


class SubfigureGOld3(Subfigure):
    """[旧版] Panel G 散点网格: 5 表型各取 top1 代谢物画散点+回归线。
    已被认知热图版 (SubfigureG) 取代, 保留备查。"""
    subfigure_label = 'g_old3'
    subfigure_title = 'correlation_scatter_old3'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('phenotype_corr')
        if not cache.exists():
            print('Panel G: 未找到表型相关缓存，请先运行 phenotype_corr 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        cached = pickle_load(cache)
        corr = cached.get('correlations', {}).get('hippocampus')
        pheno = cached.get('pheno')
        if corr is None or corr.empty or pheno is None:
            print('Panel G: 表型相关数据为空')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return
        ds = _get_ds()
        _sig_p = corr[corr['p_value'] < 0.05]
        print(f'=== Panel G: hippocampus 表型相关 (raw p<0.05: {len(_sig_p)} 对, '
              f'{_sig_p["metabolite_id"].nunique()} 代谢物) ===')
        # 主图: 每表型 top1, 取 4 个最好
        main_selections = _top_corr_selections(corr, rank=1, top_n_phenos=4)
        # 次相关 (top2) → figure s4
        sub_selections = _top_corr_selections(corr, rank=2, top_n_phenos=4)
        print(f'  主图 top1: {main_selections}')
        print(f'  s4 top2:   {sub_selections}')
        panels = _build_corr_panels(ds, corr, pheno, 'hippocampus', main_selections, GROUP_COLORS)
        if not panels:
            print('Panel G: 无有效散点数据')
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
# Panel G (旧版2): 脂质代谢物 × 表型相关热图 (设计方案 v3, 合并 treated+veh)
# ============================================================

class SubfigureGOld2(Subfigure):
    subfigure_label = 'g_old2'
    subfigure_title = 'lipid_phenotype_correlation'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        met_ids = LIPID_METABOLITES_FOR_DEF
        if not met_ids:
            print('Panel G: LIPID_METABOLITES_FOR_DEF 为空')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        cache_corr = _cache_path('phenotype_corr')
        if not cache_corr.exists():
            print('Panel G: 未找到表型相关缓存，请先运行 phenotype_corr 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        cached = pickle_load(cache_corr)
        hippocampus_corr = cached.get('correlations', {}).get('hippocampus')
        if hippocampus_corr is None or hippocampus_corr.empty:
            print('Panel G: 海马表型相关数据为空')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 打印全量表型相关统计 (画图时输出显著分子数 + top 10)
        _full = hippocampus_corr
        _sig_p = _full[_full['p_value'] < 0.05]
        _adj_col = ('adj_p' if 'adj_p' in _full.columns
                    else ('adj_p_value' if 'adj_p_value' in _full.columns else None))
        _sig_adj = _full[_full[_adj_col] < 0.05] if _adj_col else _full.iloc[0:0]
        print(f'=== Panel G: hippocampus 表型相关 (全量) ===')
        print(f'  总: {len(_full)} 对 ({_full["metabolite_id"].nunique()} 代谢物 '
              f'× {_full["phenotype"].nunique()} 表型)')
        print(f'  raw p<0.05: {len(_sig_p)} 对, {_sig_p["metabolite_id"].nunique()} 代谢物')
        if _adj_col:
            print(f'  {_adj_col}<0.05: {len(_sig_adj)} 对, {_sig_adj["metabolite_id"].nunique()} 代谢物')
        _sort_col = _adj_col or 'p_value'
        _name_col = 'metabolite_name' if 'metabolite_name' in _full.columns else 'metabolite_id'
        _top = _full.nsmallest(10, _sort_col)
        print(f'  top 10 (by {_sort_col}):')
        for _, _r in _top.iterrows():
            print(f'    {_r["metabolite_id"]:16s} {str(_r[_name_col])[:38]:38s} '
                  f'{str(_r["phenotype"])[:18]:18s} r={_r["r"]:+.3f} p={_r["p_value"]:.4f} '
                  f'{_sort_col}={_r[_sort_col]:.4f}')

        # 过滤到 LIPID_METABOLITES_FOR_DEF
        plot_df = hippocampus_corr[hippocampus_corr['metabolite_id'].isin(met_ids)]
        if plot_df.empty:
            print('Panel G: 脂质代谢物无表型相关数据')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 仅保留至少有一个显著相关 (p<0.05) 的代谢物
        sig_met_ids = plot_df[plot_df['p_value'] < 0.05]['metabolite_id'].unique()
        if len(sig_met_ids) > 0:
            plot_df = plot_df[plot_df['metabolite_id'].isin(sig_met_ids)]
        else:
            top_met_ids = (plot_df.groupby('metabolite_id')['p_value'].min()
                          .nsmallest(10).index.tolist())
            plot_df = plot_df[plot_df['metabolite_id'].isin(top_met_ids)]

        pheno_order = ['NOR_DI', 'NOR_探索时间', '旷场_总路程', '旷场_中心时间', 'Y迷宫_交替%']
        pheno_order = [p for p in pheno_order if p in plot_df['phenotype'].unique()]
        if not pheno_order:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        r_pivot = plot_df.pivot(index='metabolite_id', columns='phenotype', values='r')
        r_pivot = r_pivot.reindex(columns=pheno_order).fillna(0)
        r_matrix = r_pivot.values

        met_order = r_pivot.index.tolist()
        name_col = ('metabolite_name' if 'metabolite_name' in plot_df.columns else 'metabolite_id')
        name_map_corr = (plot_df.drop_duplicates('metabolite_id')
                        .set_index('metabolite_id')[name_col].to_dict())
        met_labels_corr = [str(name_map_corr.get(m, m)) for m in met_order]

        pheno_labels_en = [PHENO_NAME_EN.get(p, p) for p in pheno_order]

        params = {
            'r_matrix': r_matrix,
            'metabolite_labels': met_labels_corr,
            'phenotype_labels': pheno_labels_en,
            'tissue': 'hippocampus',
        }

        scale = 0.57
        element = Elements.PhenotypeCorrelationFigure(**{
            ParameterName.figure_data_parameter_dict: params,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Panel G (旧版): 代谢物–表型相关
# ============================================================

class SubfigureGOld(Subfigure):
    subfigure_label = 'g'
    subfigure_title = 'phenotype_correlation'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('phenotype_corr')
        if not cache.exists():
            print('Panel G: 未找到表型相关缓存，请先运行 phenotype_corr 分析')
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        cached = pickle_load(cache)
        correlations = cached['correlations']

        # 取海马的相关结果
        tissue = 'hippocampus' if 'hippocampus' in correlations else next(iter(correlations), None)
        if tissue is None:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        corr_df = correlations[tissue]
        if corr_df.empty:
            super().__init__({}, subfigure_bottom_left, subfigure_size,
                             subfigure_label=self.subfigure_label,
                             subfigure_title=self.subfigure_title, background=False)
            return

        # 构建 r_matrix 和标签
        sig_mets = corr_df[corr_df['significant']]['metabolite_id'].unique()
        if len(sig_mets) == 0:
            sig_mets = corr_df.nsmallest(15, 'adj_p')['metabolite_id'].unique()
        plot_df = corr_df[corr_df['metabolite_id'].isin(sig_mets)]

        pheno_order = ['NOR_DI', 'NOR_探索时间', '旷场_总路程', '旷场_中心时间', 'Y迷宫_交替%']
        pheno_order = [p for p in pheno_order if p in plot_df['phenotype'].unique()]

        met_avg_r = plot_df.groupby('metabolite_id')['r'].apply(lambda x: x.abs().mean())
        met_order = met_avg_r.sort_values(ascending=False).index.tolist()

        r_matrix = plot_df.pivot(index='metabolite_id', columns='phenotype', values='r')
        r_matrix = r_matrix.reindex(index=met_order, columns=pheno_order).fillna(0).values

        name_col = 'metabolite_name' if 'metabolite_name' in plot_df.columns else 'metabolite_id'
        name_map = plot_df.drop_duplicates('metabolite_id').set_index('metabolite_id')[name_col].to_dict()
        met_labels = [str(name_map.get(m, m)) for m in met_order]

        pheno_labels_en = [PHENO_NAME_EN.get(p, p) for p in pheno_order]
        figure_data_parameter_dict = {
            'r_matrix': r_matrix,
            'metabolite_labels': met_labels,
            'phenotype_labels': pheno_labels_en,
            'tissue': tissue,
        }
        probe = Elements.PhenotypeCorrelationFigure(figure_data_parameter_dict=dict(figure_data_parameter_dict))
        scale = float(np.min(subfigure_size / probe.size))
        element = Elements.PhenotypeCorrelationFigure(**{
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
            ParameterName.scale: scale,
            ParameterName.bottom_left_offset: subfigure_bottom_left,
        })
        move_center(element, scale, subfigure_size)

        subfigure_element_dict = {element.name: element}
        super().__init__(
            subfigure_element_dict, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


def _export_used_differential(ds):
    """导出 differential_volcano.xlsx: 每组织被 figure 真正用到的代谢物
    (Panel B/C 显著差异 ∪ Panel D LIPID_METABOLITES_FOR_DEF ∪ Panel F F2_FOCUS/PATHWAY_NODE)。
    非全量; 含手动列表中不显著的 (用到即收)。"""
    from pathlib import Path as _P
    manual = (set(LIPID_METABOLITES_FOR_DEF) | set(F2_FOCUS_METABOLITES)
              | {m for m, *_ in PATHWAY_NODE_DEFS if m})
    fc_log2 = np.log2(1.5)
    p_thr = 0.05
    info = ds.metabolite_info
    info_cols = [c for c in ['Compounds', '物质', 'Class I', 'Class II', 'Level', 'HMDB', 'KEGG (cpd_ID)']
                 if c in info.columns]
    tissues = ['cortex', 'hippocampus', 'hypothalamus', 'liver', 'kidney', 'fat',
               'gastrocnemius', 'quadriceps', 'ileum', 'serum']
    out_dir = _P(__file__).parents[3] / 'output' / 'metabolomics' / 'volcano'
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / 'differential_volcano.xlsx'
    n_total = 0
    with pd.ExcelWriter(out, engine='openpyxl') as writer:
        for t in tissues:
            cache = _cache_path('diff', tissue=t)
            if not cache.exists():
                continue
            diff = pickle_load(cache)
            if diff is None or diff.empty:
                continue
            sig_arr = ((diff['p_value'] < p_thr) & (diff['log2FC'].abs() >= fc_log2)
                       & diff['p_value'].notna()).values
            manual_arr = diff.index.isin(manual)
            used = sig_arr | manual_arr
            if not used.any():
                continue
            sub = diff.loc[used].copy()
            if info_cols:
                sub = sub.join(info[info_cols], how='left')
            sub['significant'] = sig_arr[used]
            sub['in_manual_list'] = manual_arr[used]
            front = ['Compounds', '物质', 'Class I', 'log2FC', 'p_value', 'adj_p_value',
                     'significant', 'in_manual_list']
            ordered = [c for c in front if c in sub.columns] + [c for c in sub.columns if c not in front]
            sub = sub.sort_values(['significant', 'in_manual_list', 'p_value'],
                                  ascending=[False, False, True], na_position='last')
            sub[ordered].to_excel(writer, sheet_name=t[:31])
            n_total += len(sub)
            print(f'  {t}: 用到 {len(sub)} (显著 {int(sig_arr[used].sum())}, '
                  f'手动列表 {int(manual_arr[used].sum())})')
    print(f'已导出 (仅 figure 用到的): {out.name}  共 {n_total} 行')


# ============================================================
# Figure 4
# ============================================================

class Figure(GeneralFigure):
    figure_label = 'figure_4'
    figure_title = 'Figure 4'

    def __init__(self):
        subfigure_class_list = [
            SubfigureA,
            SubfigureB,       # b — 化学类富集棒棒糖
            SubfigureC,       # c — 器官×类富集热图
            SubfigureD,       # d — 脂质方向热图
            SubfigureE,       # e — 海马聚类热图
            SubfigureF1,      # f — 收窄柱图 (6 焊点分子, 脑+血浆); F1 通路图备用
            SubfigureG,       # g — 脂质×表型相关
        ]

        figure_layout_list = [
            (0.30, [(1.0, 'a')]),              # A: OPLS-DA 网格
            (0.25, [(0.5, 'b'), (0.5, 'c')]),  # B: 富集棒棒糖 + C: 器官×类热图
            (0.25, [(0.5, 'd'), (0.5, 'e')]),  # D: 脂质方向热图 + E: 海马聚类热图
            (0.20, [(0.5, 'f'), (0.5, 'g')]),  # F2: 收窄柱图 + G: 脂质×表型相关
        ]
        super().__init__(
            self.figure_label, subfigure_class_list, figure_layout_list, {},
            figure_title=self.figure_title)
        _export_used_differential(_get_ds())
