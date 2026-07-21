# -*- coding: utf-8 -*-
"""Figure S4: 代谢组学补充图

S4A 数据质控 / S4B 模型验证(置换检验) / S4C 皮层深度 / S4D 其余器官广度 / S4G UpSet 交集
跳过 S4E/S4F/S4H (无数据)。模式与 figure_4.py 完全一致 (Subfigure + Figure + probe/move_center)。
"""

import numpy as np
import pandas as pd

from ...common.config import (
    ParameterName, Figure as GeneralFigure, Vector, calculate_center_bottom_offset)
from ...common.elements import Elements
from ...common.common_figure_materials import (
    CommonColorDict as CCD, gray_teal_diverging_cmap, is_lipid_class, lipid_class_color)
from scripts.src.common_parallel_solver.config import pickle_load
from scripts.src.metabolomics_analysis.metabolomics_main_functions import (
    _cache_path, _load_dataset)
from scripts.src.inventory import MetabolomicsExperimentName
from .figure_4 import _build_corr_panels, _top_corr_selections, GROUP_COLORS, PHENO_NAME_EN

Subfigure = Elements.Subfigure

TISSUES = ['cortex', 'hippocampus', 'hypothalamus', 'liver', 'kidney',
           'fat', 'gastrocnemius', 'quadriceps', 'ileum', 'serum']


def move_center(element_obj, scale, subfigure_size, calculate_center_param_dict=None, extra_offset=None):
    if calculate_center_param_dict is None:
        calculate_center_param_dict = {}
    if extra_offset is None:
        extra_offset = Vector(0., 0.)
    center = element_obj.calculate_center(element_obj, scale, **calculate_center_param_dict)
    center_bottom_offset = calculate_center_bottom_offset(center, subfigure_size) + extra_offset
    element_obj.move_and_scale(bottom_left_offset=center_bottom_offset)


_DATASET_CACHE = {}
def _get_ds():
    if 'ds' not in _DATASET_CACHE:
        _DATASET_CACHE['ds'] = _load_dataset(MetabolomicsExperimentName.data_20260412_metabolomics)
    return _DATASET_CACHE['ds']


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


def _empty_subfigure(self, subfigure_bottom_left, subfigure_size):
    super(self.__class__, self).__init__(
        {}, subfigure_bottom_left, subfigure_size,
        subfigure_label=self.subfigure_label,
        subfigure_title=self.subfigure_title, background=False)


# ============================================================
# S4A: 数据质控 — QC RSD + PCA
# ============================================================

class SubfigureA(Subfigure):
    subfigure_label = 'a'
    subfigure_title = 'qc'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        ds = _get_ds()
        qc_cols = [c for c in ds.sample_data.columns
                   if ds.sample_info.get(c, {}).get('is_qc')]
        if not qc_cols:
            print('S4A: 无 QC 样本')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        qc_data = ds.sample_data[qc_cols].apply(pd.to_numeric, errors='coerce')
        with np.errstate(divide='ignore', invalid='ignore'):
            mean = qc_data.mean(axis=1, skipna=True)
            std = qc_data.std(axis=1, skipna=True)
            rsd = std / mean * 100
        qc_rsd = rsd.replace([np.inf, -np.inf], np.nan).dropna().values
        # PCA (numpy SVD, QC 应聚集)
        all_cols = ds.sample_data.columns.tolist()
        M = ds.sample_data.apply(pd.to_numeric, errors='coerce').fillna(0).T
        is_qc = np.array([ds.sample_info.get(c, {}).get('is_qc', False) for c in all_cols])
        X = M - M.mean(axis=0)
        try:
            U, S, Vt = np.linalg.svd(X, full_matrices=False)
            pc1 = U[:, 0] * S[0]
            pc2 = U[:, 1] * S[1] if len(S) > 1 else np.zeros_like(pc1)
        except Exception:
            pc1, pc2 = np.zeros(len(all_cols)), np.zeros(len(all_cols))
        element = _place(
            Elements.QCFigure,
            {'qc_rsd': qc_rsd, 'pca': {'pc1': pc1, 'pc2': pc2, 'is_qc': is_qc}},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# S4B: 模型验证 — OPLS-DA 置换检验网格
# ============================================================

class SubfigureB(Subfigure):
    subfigure_label = 'b'
    subfigure_title = 'permutation_test'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('oplsda_all')
        if not cache.exists():
            print('S4B: 未找到 oplsda_all 缓存，请先运行 oplsda 分析')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        perm_results = pickle_load(cache).get('perm_results', {})
        if not perm_results:
            print('S4B: oplsda_all 无 perm_results')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        element = _place(
            Elements.PermutationGridFigure,
            {'perm_results': perm_results, 'tissue_order': TISSUES},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# S4C: 皮层深度 — cortex 火山图 (复用 MetabolomicsVolcanoFigure)
# ============================================================

class SubfigureC(Subfigure):
    subfigure_label = 'c'
    subfigure_title = 'cortex_volcano'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('volcano', tissue='cortex')
        if not cache.exists():
            print('S4C: 未找到 cortex 火山图缓存，请先运行 volcano -t cortex')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        vol_data = pickle_load(cache)
        t = 'cortex'
        if t not in vol_data:
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        d = vol_data[t]
        ds = _get_ds()
        if 'Class I' in ds.metabolite_info.columns:
            class_series = ds.metabolite_info['Class I'].reindex(d['index'])
            classes = [class_series.iloc[i] for i in range(len(d['index']))]
        else:
            classes = ['Other'] * len(d['index'])
        element = _place(
            Elements.MetabolomicsVolcanoFigure,
            {'log2fc': d['log2fc'], 'neg_log10p': d['neg_log10p'],
             'classes': classes, 'names': d['names'],
             'highlight_indices': [], 'thresholds': d['thresholds']},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Figure S4
# ============================================================

# ============================================================
# S4D: 其余器官广度 — 下丘脑/回肠/肾/股四头肌 火山图网格
# ============================================================

class SubfigureD(Subfigure):
    subfigure_label = 'd'
    subfigure_title = 'other_tissue_volcano'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        s4d_tissues = ['hypothalamus', 'ileum', 'kidney', 'quadriceps']
        volcano_data, thresholds = {}, None
        for t in s4d_tissues:
            cache = _cache_path('volcano', tissue=t)
            if not cache.exists():
                continue
            vd = pickle_load(cache)
            if t not in vd:
                continue
            d = vd[t]
            volcano_data[t] = {'log2fc': d['log2fc'], 'neg_log10p': d['neg_log10p']}
            thresholds = d['thresholds']
        if not volcano_data:
            print('S4D: 未找到其余组织火山图缓存')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        element = _place(
            Elements.MultiVolcanoGridFigure,
            {'volcano_data': volcano_data, 'tissues': s4d_tissues, 'thresholds': thresholds},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# S4G: 数据表 + 跨器官 UpSet 交集
# ============================================================

class SubfigureG(Subfigure):
    subfigure_label = 'g'
    subfigure_title = 'upset_intersection'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        fc_thr = np.log2(1.5)
        p_thr = 0.05
        sets = {}
        for t in TISSUES:
            cache = _cache_path('diff', tissue=t)
            if not cache.exists():
                continue
            diff = pickle_load(cache)
            if diff is None or diff.empty:
                continue
            sig = diff[(diff['p_value'] < p_thr) & (np.abs(diff['log2FC']) >= fc_thr)
                       & diff['p_value'].notna()]
            if len(sig):
                sets[t] = set(sig.index)
        if not sets:
            print('S4G: 未找到各组织差异分析缓存')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        element = _place(
            Elements.UpSetFigure,
            {'tissues': TISSUES, 'sets': sets, 'top_n': 15},
            subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


# ============================================================
# Figure S4
# ============================================================

# ============================================================
# S4H: 次相关散点 (4 个最好表型的 top2 代谢物)
# ============================================================

class SubfigureH(Subfigure):
    """S4H: 主图 G 的 4 个表型各自的 top2 (次相关) 代谢物散点+回归线。"""
    subfigure_label = 'h'
    subfigure_title = 'correlation_scatter_s4'

    def __init__(self, subfigure_bottom_left=None, subfigure_size=None):
        cache = _cache_path('phenotype_corr')
        if not cache.exists():
            print('S4H: 未找到表型相关缓存')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        cached = pickle_load(cache)
        corr = cached.get('correlations', {}).get('hippocampus')
        pheno = cached.get('pheno')
        if corr is None or corr.empty or pheno is None:
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        ds = _get_ds()
        selections = _top_corr_selections(corr, rank=2, top_n_phenos=4)
        panels = _build_corr_panels(ds, corr, pheno, 'hippocampus', selections, GROUP_COLORS)
        if not panels:
            print('S4H: 无有效散点数据')
            _empty_subfigure(self, subfigure_bottom_left, subfigure_size)
            return
        element = _place(Elements.CorrelationScatterGridFigure, {'panels': panels},
                         subfigure_bottom_left, subfigure_size)
        super().__init__(
            {element.name: element}, subfigure_bottom_left, subfigure_size,
            subfigure_label=self.subfigure_label,
            subfigure_title=self.subfigure_title, background=False)


class Figure(GeneralFigure):
    figure_label = 'figure_s4'
    figure_title = 'Figure S4'

    def __init__(self):
        subfigure_class_list = [
            SubfigureA,
            SubfigureB,
            SubfigureC,
            SubfigureD,
            SubfigureG,
            SubfigureH,
        ]
        figure_layout_list = [
            (0.18, [(1.0, 'a')]),
            (0.22, [(1.0, 'b')]),
            (0.15, [(1.0, 'c')]),
            (0.18, [(1.0, 'd')]),
            (0.12, [(1.0, 'g')]),
            (0.15, [(1.0, 'h')]),
        ]
        super().__init__(
            self.figure_label, subfigure_class_list, figure_layout_list, {},
            figure_title=self.figure_title)
