# -*- coding: utf-8 -*-
"""代谢组学数据加载器

读取 data_20260412_metabolomics 文件夹下的所有 xlsx 文件，
整合代谢物信息和样本数据，方便后续分析。
"""

import re
from pathlib import Path

import pandas as pd

from ..data_20251115_ad_sema.specific_data_parameters import Keyword

# ── 代谢物注释黑名单 (合并版: 中文版 + 补充, 中英对照, 跨 Tier 去重) ──
_BLACKLIST_FILENAME = '代谢物注释黑名单_Metabolite_Annotation_Blacklist.xlsx'
_BLACKLIST_NAME_COL = 'Compounds (英文注释名) / Compound'
_BLACKLIST_CORRECT_NAME_COL = '正确中文名 / Correct Name (CN)'
_BLACKLIST_SHEETS = [
    '01_Tier1_建议移除_Remove',
    '02_Tier2_污染物_Contam',
    '03_Tier3_饮食微生物_Diet',
    '04_Tier4_低可信_LowConf',
]
_EXEMPT_SHEET = '04b_脂质豁免_Lipids_Exempt'
_NAME_FIX_SHEET = '05_名称错配_Name_Mismatch'
# 不作为数据文件加载的 xlsx 文件名关键字
_EXCLUDE_FILE_KEYWORDS = ('行为学', 'Blacklist', '黑名单', '使用说明', 'backup')

# 元数据列名（所有文件共享的前 19 列）
METADATA_COLUMNS = [
    'Index', 'Mode', 'Compounds', '物质',
    'Class I', '物质一级分类', 'Class II', '物质二级分类',
    'Q1 (Da)', 'Molecular weight (Da)', 'Ionization model', 'Formula',
    'Level', 'CAS', 'PubChem CID', 'HMDB', 'ChEBI', 'Metlin', 'cpd_ID',
]

QC_PATTERN = re.compile(r'^QC\d+$')
SHARED_SAMPLE_PATTERN = re.compile(r'^(C\d+)-(\w+)-(C\d+)-\2$')
NORMAL_SAMPLE_PATTERN = re.compile(r'^(C\d+)-(.+)$')


def _parse_sample_name(sample_name: str):
    """解析样本列名

    Returns:
        (animal_id, tissue, is_qc, shared_animals)

    Examples:
        'C11-serum'                              -> ('C11', 'serum', False, [])
        'C11-gastrocnemius-C12-gastrocnemius'     -> ('C11', 'gastrocnemius', False, ['C12'])
        'QC01'                                    -> (None, None, True, [])
        'a3_all_sample_data__QC01'                -> (None, None, True, [])
    """
    # 去掉来源文件前缀（QC 列去重时添加的 {file_stem}__ 前缀）
    clean = sample_name.rsplit('__', 1)[-1] if '__' in sample_name else sample_name
    if QC_PATTERN.match(clean):
        return None, None, True, []

    m = SHARED_SAMPLE_PATTERN.match(sample_name)
    if m:
        return m.group(1), m.group(2), False, [m.group(3)]

    m = NORMAL_SAMPLE_PATTERN.match(sample_name)
    if m:
        return m.group(1), m.group(2), False, []

    return None, sample_name, False, []


class MetabolomicDataset:
    """代谢组学数据集

    整合多个 xlsx 文件的代谢组学数据，提供统一的数据访问接口。

    Usage:
        ds = MetabolomicDataset()
        print(ds.summary())

        # 设置分组后按组筛选
        ds.set_groups({'C11': 'AD', ..., 'C42': 'AD+Sema', ...})
        cortex_ad = ds.get_data_by_group('AD', tissue='cortex')

        # 或直接按组织/动物筛选
        serum = ds.get_data_by_tissue('serum')
        subset = ds.get_data_by_animals(['C11', 'C12'], tissue='cortex')
    """

    # MSI 鉴定等级数值化映射（数值越小置信越高）
    _LEVEL_RANK = {'1a': 1, '1b': 2, '2': 3, '3': 4}

    def __init__(self, data_dir: str | Path | None = None, min_msi_level: str = '3'):
        self._data_dir = Path(data_dir) if data_dir else Path(__file__).parent
        self.metabolite_info: pd.DataFrame = None  # type: ignore
        self.sample_data: pd.DataFrame = None  # type: ignore
        self.sample_info: dict[str, dict] = {}
        self.source_map: dict[str, str] = {}
        self._group_map: dict[str, str] = {
            'C11': Keyword.ad, 'C12': Keyword.ad, 'C13': Keyword.ad,
            'C14': Keyword.ad, 'C15': Keyword.ad, 'C18': Keyword.ad,
            'C42': Keyword.sema, 'C43': Keyword.sema, 'C44': Keyword.sema,
            'C45': Keyword.sema, 'C46': Keyword.sema, 'C47': Keyword.sema,
        }
        self._load_all()
        self.apply_blacklist()
        if min_msi_level != '3':
            self.filter_by_msi_level(min_msi_level)

    def _load_all(self):
        xlsx_files = sorted(
            f for f in self._data_dir.glob('*.xlsx')
            if not f.name.startswith('~')
            and not any(kw in f.name for kw in _EXCLUDE_FILE_KEYWORDS)
        )
        if not xlsx_files:
            raise FileNotFoundError(f'{self._data_dir} 下未找到 xlsx 文件')

        info_frames = []
        data_frames = []

        all_sample_cols: set[str] = set()

        for fpath in xlsx_files:
            df = pd.read_excel(fpath, sheet_name=0)
            meta_cols = [c for c in df.columns if c in METADATA_COLUMNS or c == 'kegg_map']
            sample_cols = [c for c in df.columns if c not in meta_cols]

            # 对已有文件中出现过的列名加上来源前缀（避免 QC 列冲突）
            rename_map = {}
            for col in sample_cols:
                if col in all_sample_cols:
                    new_name = f'{fpath.stem}__{col}'
                    rename_map[col] = new_name
                    col = new_name
                else:
                    all_sample_cols.add(col)

                animal_id, tissue, is_qc, shared = _parse_sample_name(col)
                self.sample_info[col] = {
                    'animal_id': animal_id,
                    'tissue': tissue,
                    'is_qc': is_qc,
                    'shared_animals': shared,
                }
                self.source_map[col] = fpath.name

            info_frames.append(df[meta_cols].copy())

            data = df[sample_cols].rename(columns=rename_map)
            data.index = df['Index']
            data = data[~data.index.duplicated(keep='first')]
            data_frames.append(data)

        # 合并代谢物元信息
        self.metabolite_info = (
            pd.concat(info_frames, ignore_index=True)
            .drop_duplicates(subset='Index', keep='first')
            .set_index('Index')
        )

        # 合并样本丰度（按 Index 外连接，保留所有代谢物）
        self.sample_data = data_frames[0]
        for df in data_frames[1:]:
            self.sample_data = self.sample_data.join(df, how='outer')

    def __repr__(self):
        tissues = ', '.join(self.list_tissues())
        return (
            f'MetabolomicDataset('
            f'metabolites={len(self.metabolite_info)}, '
            f'samples={sum(1 for i in self.sample_info.values() if not i["is_qc"])}, '
            f'tissues=[{tissues}])'
        )

    # ---- 分组管理 ----

    def set_groups(self, group_map: dict[str, str]):
        """设置动物分组映射 {animal_id: group_name}"""
        self._group_map = dict(group_map)

    def get_sample_groups(self, include_qc: bool = False) -> pd.Series:
        """返回每个样本的分组标签"""
        records = {}
        for name, info in self.sample_info.items():
            if info['is_qc']:
                records[name] = 'QC' if include_qc else None
            else:
                records[name] = self._group_map.get(info['animal_id'])
        return pd.Series(records, name='group').dropna()

    # ---- 数据子集 ----

    def get_data_by_tissue(self, tissue: str, exclude_qc: bool = True) -> pd.DataFrame:
        """获取指定组织类型的数据

        Args:
            tissue: 组织名，如 'serum', 'cortex', 'hippocampus'
            exclude_qc: 是否排除 QC 样本
        """
        cols = [
            name for name, info in self.sample_info.items()
            if info['tissue'] == tissue
            and (not info['is_qc'] or not exclude_qc)
            and name in self.sample_data.columns
        ]
        return self.sample_data[cols]

    def get_data_by_group(self, group: str, tissue: str | None = None) -> pd.DataFrame:
        """获取指定分组（和可选组织）的数据"""
        return self.sample_data[self._select_columns(group=group, tissue=tissue)]

    def get_data_by_animals(self, animal_ids: list[str], tissue: str | None = None) -> pd.DataFrame:
        """获取指定动物的数据（含共享样本）"""
        return self.sample_data[self._select_columns(animal_ids=animal_ids, tissue=tissue)]

    def get_data_by_source(self, source_file: str) -> pd.DataFrame:
        """获取指定来源文件的数据"""
        cols = [
            name for name, src in self.source_map.items()
            if src == source_file and name in self.sample_data.columns
        ]
        return self.sample_data[cols]

    def _select_columns(self, *, group: str | None = None,
                        animal_ids: list[str] | None = None,
                        tissue: str | None = None) -> list[str]:
        animal_set = set(animal_ids) if animal_ids else None
        cols = []
        for name, info in self.sample_info.items():
            if info['is_qc'] or name not in self.sample_data.columns:
                continue
            if animal_set is not None:
                primary = info['animal_id'] in animal_set
                shared = any(a in animal_set for a in info['shared_animals'])
                if not (primary or shared):
                    continue
            if group is not None:
                if self._group_map.get(info['animal_id']) != group:
                    continue
            if tissue is not None and info['tissue'] != tissue:
                continue
            cols.append(name)
        return cols

    # ---- 代谢物筛选 ----

    def filter_metabolites(self, compound: str | None = None,
                           class_i: str | None = None,
                           class_ii: str | None = None,
                           cpd_id: str | None = None,
                           kegg_map: str | None = None) -> pd.DataFrame:
        """按代谢物属性筛选，返回匹配的代谢物元信息"""
        mask = pd.Series(True, index=self.metabolite_info.index)
        info = self.metabolite_info
        if compound is not None:
            mask &= info['Compounds'].str.contains(compound, case=False, na=False)
        if class_i is not None:
            mask &= info['Class I'].str.contains(class_i, case=False, na=False)
        if class_ii is not None:
            mask &= info['Class II'].str.contains(class_ii, case=False, na=False)
        if cpd_id is not None:
            mask &= info['cpd_ID'].str.contains(cpd_id, na=False)
        if kegg_map is not None:
            mask &= info['kegg_map'].str.contains(kegg_map, na=False)
        return info[mask]

    def apply_blacklist(self) -> 'MetabolomicDataset':
        """应用代谢物注释黑名单

        从 代谢物注释黑名单.xlsx (合并中文版+补充, 中英对照) 读取黑名单，
        移除不准确注释的代谢物，并修正 05 名称错配 sheet 中记录的名字错配。
        04b 脂质豁免 sheet 中的脂质代谢物不会被移除。
        """
        bl_path = self._data_dir / _BLACKLIST_FILENAME
        if not bl_path.exists():
            return self

        # 收集黑名单代谢物名称（Tier 1-4，已合并补充并跨 Tier 去重）
        blacklisted: set[str] = set()
        for sheet in _BLACKLIST_SHEETS:
            df = pd.read_excel(bl_path, sheet_name=sheet)
            blacklisted.update(df[_BLACKLIST_NAME_COL].dropna().astype(str))

        # 收集豁免代谢物名称（04b 脂质除外）
        exempt_df = pd.read_excel(bl_path, sheet_name=_EXEMPT_SHEET)
        exempt: set[str] = set(exempt_df[_BLACKLIST_NAME_COL].dropna().astype(str))

        # 有效黑名单 = 黑名单 - 豁免
        effective = blacklisted - exempt

        # 移除黑名单代谢物
        keep_mask = ~self.metabolite_info['Compounds'].isin(effective)
        keep_idx = self.metabolite_info.index[keep_mask]
        before = len(self.metabolite_info)
        self.metabolite_info = self.metabolite_info.loc[keep_idx]
        self.sample_data = self.sample_data.loc[self.sample_data.index.intersection(keep_idx)]
        after = len(self.metabolite_info)
        n_removed = before - after
        if n_removed:
            print(f'黑名单过滤: {before} → {after} 代谢物 (移除 {n_removed} 个注释不准确的代谢物)')

        # 修正名字错配
        name_fix_df = pd.read_excel(bl_path, sheet_name=_NAME_FIX_SHEET)
        fix_map = dict(zip(
            name_fix_df[_BLACKLIST_NAME_COL].astype(str),
            name_fix_df[_BLACKLIST_CORRECT_NAME_COL].astype(str),
        ))
        n_fixed = 0
        for compound, correct_name in fix_map.items():
            mask = self.metabolite_info['Compounds'] == compound
            if mask.any():
                self.metabolite_info.loc[mask, '物质'] = correct_name
                n_fixed += 1
        if n_fixed:
            print(f'名字修正: 修正了 {n_fixed} 个代谢物的中文名')

        return self

    def filter_by_msi_level(self, max_level: str = '2') -> 'MetabolomicDataset':
        """按 MSI 鉴定置信等级过滤代谢物

        保留 Level 数值 ≤ max_level 的代谢物。
        等级排序: 1a(1) > 1b(2) > 2(3) > 3(4)

        Args:
            max_level: 最高允许等级 ('1a', '1b', '2', '3')

        Returns:
            self（支持链式调用）
        """
        if 'Level' not in self.metabolite_info.columns:
            return self
        max_rank = self._LEVEL_RANK.get(max_level, 4)
        level_ranks = self.metabolite_info['Level'].map(self._LEVEL_RANK).fillna(4)
        keep = level_ranks <= max_rank
        keep_idx = self.metabolite_info.index[keep]
        before = len(self.metabolite_info)
        self.metabolite_info = self.metabolite_info.loc[keep_idx]
        self.sample_data = self.sample_data.loc[self.sample_data.index.intersection(keep_idx)]
        after = len(self.metabolite_info)
        print(f'MSI Level 过滤 (≤{max_level}): {before} → {after} 代谢物 '
              f'(移除 {before - after} 个 Level > {max_level})')
        return self

    # ---- 汇总信息 ----

    def list_tissues(self) -> list[str]:
        return sorted({i['tissue'] for i in self.sample_info.values() if i['tissue']})

    def list_animals(self) -> list[str]:
        return sorted({i['animal_id'] for i in self.sample_info.values() if i['animal_id']})

    def list_sources(self) -> dict[str, list[str]]:
        sources: dict[str, list[str]] = {}
        for name, src in self.source_map.items():
            sources.setdefault(src, []).append(name)
        return sources

    def summary(self) -> str:
        n_metabolites = len(self.metabolite_info)
        n_samples = sum(1 for i in self.sample_info.values() if not i['is_qc'])
        n_qc = sum(1 for i in self.sample_info.values() if i['is_qc'])
        tissues = ', '.join(self.list_tissues())
        animals = ', '.join(self.list_animals())
        sources = ', '.join(self.list_sources().keys())
        return (
            f'代谢物数量: {n_metabolites}\n'
            f'样本数量: {n_samples} (QC: {n_qc})\n'
            f'组织类型: {tissues}\n'
            f'动物 ID: {animals}\n'
            f'来源文件: {sources}'
        )


def load_metabolomic_data(data_dir: str | Path | None = None) -> MetabolomicDataset:
    """快捷函数：加载代谢组学数据"""
    return MetabolomicDataset(data_dir)
