# -*- coding: utf-8 -*-
"""行为学实验数据加载器（共享: 代谢组小鼠 + 代谢流小鼠 两拨）

通过 load_behavioral_data(data_source) 分流:
- 'metabolomics': 代谢组小鼠 (C11-C18/C42-C47), NOR/旷场/Y迷宫/水迷宫 + derived 全
- 'flux': 代谢流小鼠 (C2-C61), 仅 Y迷宫 + 新物体识别 (多级表头, 复用 _read_multilevel);
         旷场/水迷宫/derived 后续补 (None)

多级表头兼容中英双语 (_zh_part 取中文)。数据目录用 Path(__file__).parent 定位。
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scripts.data.data_20251115_ad_sema.specific_data_parameters import Keyword

_BASE_DIR = Path(__file__).parent
_METABOLOMICS_DIR = _BASE_DIR / 'metabolomics_behavioral_data'
_FLUX_DIR = _BASE_DIR / 'flux_behavioral_data'

# 代谢组小鼠分组 (C11-C18=AD, C42-C47=Sema)
METABOLOMICS_GROUP_MAP: dict[str, str] = {
    'C11': Keyword.ad, 'C12': Keyword.ad, 'C13': Keyword.ad,
    'C14': Keyword.ad, 'C15': Keyword.ad, 'C18': Keyword.ad,
    'C42': Keyword.sema, 'C43': Keyword.sema, 'C44': Keyword.sema,
    'C45': Keyword.sema, 'C46': Keyword.sema, 'C47': Keyword.sema,
}


# ============================================================
# 通用工具
# ============================================================

def _zh_part(val) -> str:
    """表头单元格取中文部分 (双语 '中文\\nEnglish' → '中文'); 无换行则原值。"""
    if pd.isna(val):
        return ''
    return str(val).split('\n')[0].strip()


def _find_col(df: pd.DataFrame, *keywords) -> str | None:
    """按多个关键词模糊匹配列名 (所有关键词都出现); 返回首个匹配, 无则 None。"""
    for c in df.columns:
        cs = str(c)
        if all(k in cs for k in keywords):
            return c
    return None


def _read_multilevel(xlsx_path, sheet_name: str, header_rows: int = 3) -> pd.DataFrame:
    """读多级表头 sheet → 扁平列名 DataFrame (双语表头取中文部分, 合并单元格 ffill)。"""
    raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)
    hdr = raw.iloc[:header_rows]
    data = raw.iloc[header_rows:].reset_index(drop=True)
    n_cols = data.shape[1]
    rows_ff = [[_zh_part(v) for v in hdr.iloc[i].ffill()] if i < len(hdr) else [''] * n_cols
               for i in range(header_rows)]
    parent = rows_ff[header_rows - 2] if header_rows >= 2 else [''] * n_cols
    leaf = rows_ff[header_rows - 1]
    cols = []
    for i in range(n_cols):
        p = parent[i] if i < len(parent) else ''
        lf = leaf[i] if i < len(leaf) else ''
        cols.append(f'{p}_{lf}' if p and p != lf else lf)
    seen: dict = {}
    uniq = []
    for c in cols:
        n = seen.get(c, 0)
        uniq.append(f'{c}#{n}' if n else c)
        seen[c] = n + 1
    data.columns = uniq
    return data


def _normalize_id(raw: pd.DataFrame) -> pd.DataFrame:
    id_col = _find_col(raw, 'ID') or 'ID'
    raw[id_col] = raw[id_col].astype(str).str.strip().str.upper()
    raw = raw.rename(columns={id_col: 'ID'})
    return raw


def _add_group(raw: pd.DataFrame, group_map: dict) -> pd.DataFrame:
    raw['group'] = raw['ID'].apply(lambda x: group_map.get(x.upper(), 'unknown'))
    return raw


def _numeric_cols(raw: pd.DataFrame, exclude: tuple[str, ...] = ('结果', '试验', 'ID', 'group')) -> pd.DataFrame:
    for col in raw.columns:
        if col not in exclude:
            raw[col] = pd.to_numeric(raw[col], errors='coerce')
    return raw


# ============================================================
# 新物体识别 NOR
# ============================================================

def _compute_nor_derived(df: pd.DataFrame) -> pd.DataFrame:
    for region, suffix in [('鼻尖', '鼻尖'), ('中心', '中心')]:
        old_c = _find_col(df, '物体1', region, '区域时间')
        new_c = _find_col(df, '新物体1', region, '区域时间')
        if old_c and new_c:
            total = pd.to_numeric(df[old_c], errors='coerce') + pd.to_numeric(df[new_c], errors='coerce')
            df[f'DI_{suffix}'] = (pd.to_numeric(df[new_c], errors='coerce')
                                  - pd.to_numeric(df[old_c], errors='coerce')) / total.replace(0, pd.NA)
            df[f'总探索时间_{suffix}(s)'] = total
    return df


def _load_nor(xlsx_path, sheet_name, group_map) -> pd.DataFrame:
    df = _read_multilevel(xlsx_path, sheet_name, header_rows=3)
    df = _normalize_id(df)
    df = _numeric_cols(df)
    df = _add_group(df, group_map)
    df = _compute_nor_derived(df)
    return df


# ============================================================
# 旷场 Open Field (代谢组: 单行表头; flux 后续补)
# ============================================================

def _compute_of_derived(df: pd.DataFrame) -> pd.DataFrame:
    time_c = _find_col(df, '区域时间')
    if time_c:
        df['中心区域时间占比(%)'] = pd.to_numeric(df[time_c], errors='coerce') / 300 * 100
    return df


def _load_open_field(xlsx_path, sheet_name, group_map) -> pd.DataFrame:
    raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=0)
    raw.columns = [_zh_part(c) for c in raw.columns]
    raw = _normalize_id(raw)
    raw = _numeric_cols(raw)
    raw = _add_group(raw, group_map)
    raw = _compute_of_derived(raw)
    return raw


# ============================================================
# Y迷宫 Y-Maze
# ============================================================

def _compute_ym_derived(df: pd.DataFrame) -> pd.DataFrame:
    for arm in ['A', 'B', 'C', 'O']:
        col = _find_col(df, f'{arm}', '鼻尖', '区域时间') if arm != 'O' else _find_col(df, 'O', '鼻尖', '区域时间')
        if col and not _find_col(df, 'Zone', '鼻尖', '区域时间') == col:
            df[f'{arm}臂探索时间_鼻尖(s)'] = pd.to_numeric(df[col], errors='coerce')
    for region in ['鼻尖', '中心', '尾根']:
        col = _find_col(df, 'Zone_1', region, '交替百分比')
        if col:
            df[f'交替百分比_{region}(%)'] = pd.to_numeric(df[col], errors='coerce')
    for metric, label in [('轮转总次数', '轮转总次数_鼻尖'),
                          ('最大连续次数', '最大连续_鼻尖'),
                          ('直接逗留次数', '直接逗留_鼻尖'),
                          ('间接逗留次数', '间接逗留_鼻尖'),
                          ('最大次数', '最大次数_鼻尖')]:
        col = _find_col(df, 'Zone_1', '鼻尖', metric)
        if col:
            df[label] = pd.to_numeric(df[col], errors='coerce')
    return df


def _load_y_maze(xlsx_path, sheet_name, group_map) -> pd.DataFrame:
    df = _read_multilevel(xlsx_path, sheet_name, header_rows=3)
    df = _normalize_id(df)
    df = _numeric_cols(df)
    df = _add_group(df, group_map)
    df = _compute_ym_derived(df)
    return df


# ============================================================
# 水迷宫 Morris Water Maze (代谢组: gb18030 csv; flux 后续补)
# ============================================================

def _load_water_maze_csv(csv_path) -> pd.DataFrame:
    raw = pd.read_csv(csv_path, encoding='gb18030')
    raw.drop(columns=[c for c in raw.columns if c.startswith('Unnamed')], inplace=True)
    raw.rename(columns={raw.columns[-1]: 'MetaboID'}, inplace=True)
    raw['MetaboID'] = raw['MetaboID'].astype(str).str.strip().str.upper()
    raw['Daily_Mean'] = pd.to_numeric(raw['Daily_Mean'], errors='coerce')
    raw['Day_Num'] = pd.to_numeric(raw['Day_Num'], errors='coerce')
    return raw


def _compute_mwm_derived(dist_df: pd.DataFrame, latency_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for mid, g in dist_df.groupby('MetaboID'):
        g = g.sort_values('Day_Num').dropna(subset=['Day_Num', 'Daily_Mean'])
        g = g[g['Daily_Mean'] > 0]
        slope = (float(np.polyfit(g['Day_Num'].values, np.log(g['Daily_Mean'].values), 1)[0])
                 if len(g) >= 2 else np.nan)
        lat_mean = pd.to_numeric(latency_df[latency_df['MetaboID'] == mid]['Daily_Mean'],
                                 errors='coerce').mean()
        dist_mean = g['Daily_Mean'].mean()
        swim_speed = float(dist_mean / lat_mean) if (lat_mean and lat_mean > 0) else np.nan
        records.append({'MetaboID': mid, 'learn_rate': -slope, 'swim_speed': swim_speed})
    return pd.DataFrame(records).set_index('MetaboID')


# ============================================================
# flux 分组映射 (从编号对应读: control→ad, drug→sema)
# ============================================================

def _load_flux_group_map() -> dict[str, str]:
    p = _FLUX_DIR / '代谢流小鼠编号对应_Mouse_ID_Mapping.xlsx'
    df = pd.read_excel(p)
    df.columns = [_zh_part(c) for c in df.columns]
    id_col = next((c for c in df.columns if '编号' in c and '组织' not in c), None)
    grp_col = next((c for c in df.columns if 'GROUP' in c or 'Group' in c), None)
    gm = {}
    for _, row in df.iterrows():
        mid = str(row[id_col]).strip().upper()
        grp = str(row[grp_col]).strip().lower()
        gm[mid] = Keyword.ad if grp == 'control' else Keyword.sema
    return gm


def _flux_sheet_to_mouse_map() -> dict[str, str]:
    """flux 通量 sheet → flux 小鼠 ID (经器官编号中转, 真实映射链)。

    链: 小鼠ID→器官编号 (代谢流小鼠编号对应_Mouse_ID_Mapping.xlsx) +
        器官编号→(condition,index) (specific_data_parameters.py:172-177:
        ctrl_condition_list A-F→ctrl 1-6, sema_condition_list G-L→sema 1-6) →
        glucose_ad_{condition}_{index} sheet。
    """
    df = pd.read_excel(_FLUX_DIR / '代谢流小鼠编号对应_Mouse_ID_Mapping.xlsx')
    df.columns = [_zh_part(c) for c in df.columns]
    id_col = next((c for c in df.columns if '编号' in c and '组织' not in c), None)
    organ_col = next((c for c in df.columns if '组织编号' in c), None)
    # organ → sheet (cite specific_data_parameters.py ctrl/sema_condition_list)
    ctrl_organs = ['A', 'B', 'C', 'D', 'E', 'F']
    sema_organs = ['G', 'H', 'I', 'J', 'K', 'L']
    organ_to_sheet = {o: f'glucose_ad_ctrl_{i+1}' for i, o in enumerate(ctrl_organs)}
    organ_to_sheet.update({o: f'glucose_ad_semaglutide_{i+1}' for i, o in enumerate(sema_organs)})
    sheet_to_mouse = {}
    for _, row in df.iterrows():
        mouse = str(row[id_col]).strip().upper()
        organ = str(row[organ_col]).strip()
        sheet = organ_to_sheet.get(organ)
        if sheet:
            sheet_to_mouse[sheet] = mouse
    return sheet_to_mouse


# ============================================================
# 数据集类
# ============================================================

class BehavioralDataset:
    """行为学数据集（NOR + 旷场 + Y迷宫 + 水迷宫; 任意可为 None）。

    group_map: {mouse_id: Keyword.ad/sema} (代谢组或 flux 各自)。
    """

    def __init__(self, nor=None, open_field=None, y_maze=None, mwm=None, group_map=None):
        self.nor = nor
        self.open_field = open_field
        self.y_maze = y_maze
        self.mwm = mwm
        self.group_map = group_map or {}

    @property
    def animals(self) -> list[str]:
        for df in [self.nor, self.y_maze, self.open_field]:
            if df is not None and 'ID' in df.columns:
                return df['ID'].tolist()
        return []

    def get_group(self, animal_id: str) -> str:
        return self.group_map.get(animal_id.upper(), 'unknown')

    def list_groups(self) -> list[str]:
        return sorted(set(self.group_map.values()))

    def get_nor_by_group(self, group: str) -> pd.DataFrame:
        return self.nor[self.nor['group'] == group].reset_index(drop=True)

    def get_open_field_by_group(self, group: str) -> pd.DataFrame:
        return self.open_field[self.open_field['group'] == group].reset_index(drop=True)

    def get_y_maze_by_group(self, group: str) -> pd.DataFrame:
        return self.y_maze[self.y_maze['group'] == group].reset_index(drop=True)

    def summary(self) -> str:
        lines = [f'动物数: {len(self.animals)}']
        for g in self.list_groups():
            n = sum(1 for a in self.animals if self.group_map.get(a) == g)
            lines.append(f'  {g}: {n} 只')
        if self.nor is not None:
            lines.append(f'NOR 指标数: {len(self.nor.columns)}')
        if self.open_field is not None:
            lines.append(f'旷场 指标数: {len(self.open_field.columns)}')
        if self.y_maze is not None:
            lines.append(f'Y迷宫 指标数: {len(self.y_maze.columns)}')
        if self.mwm is not None:
            lines.append(f'水迷宫 动物数: {len(self.mwm)}')
        return '\n'.join(lines)


# ---- 加载函数 (分流) ----

def load_behavioral_data(data_source: str) -> BehavioralDataset:
    """加载行为学数据。data_source: 'metabolomics' | 'flux'。"""
    if data_source == 'metabolomics':
        return _load_metabolomics()
    if data_source == 'flux':
        return _load_flux()
    raise ValueError(f"未知 data_source: {data_source} (期望 'metabolomics' | 'flux')")


def _load_metabolomics() -> BehavioralDataset:
    """代谢组小鼠: 单 xlsx (NOR/旷场/Y迷宫) + 水迷宫 csv。全 derived。"""
    xlsx = _METABOLOMICS_DIR / '代谢组学小鼠行为学结果.xlsx'
    gm = METABOLOMICS_GROUP_MAP
    nor = _load_nor(xlsx, '新物体识别', gm)
    open_field = _load_open_field(xlsx, '旷场', gm)
    y_maze = _load_y_maze(xlsx, 'Y迷宫', gm)
    dist_csv = _METABOLOMICS_DIR / 'clean_raw_distance_subset.csv'
    lat_csv = _METABOLOMICS_DIR / 'clean_raw_lantency_subset.csv'
    if dist_csv.exists() and lat_csv.exists():
        mwm = _compute_mwm_derived(_load_water_maze_csv(dist_csv), _load_water_maze_csv(lat_csv))
    else:
        print(f'未找到水迷宫数据 ({dist_csv.name} / {lat_csv.name}), 跳过 MWM')
        mwm = None
    return BehavioralDataset(nor=nor, open_field=open_field, y_maze=y_maze, mwm=mwm, group_map=gm)


def _load_flux_open_field(xlsx_path, sheet_name, group_map) -> pd.DataFrame:
    """flux 旷场 (多级表头: 中心 movement + 中心区/边缘区 区域时间)。"""
    df = _read_multilevel(xlsx_path, sheet_name, header_rows=3)
    df = _normalize_id(df)
    df = _numeric_cols(df)
    df = _add_group(df, group_map)
    # 中心区/边缘区时间 (焦虑指标) + 中心区占比
    for zone_key, zone_label in [('中心区', '中心区时间'), ('边缘区', '边缘区时间')]:
        col = _find_col(df, zone_key, '区域时间')
        if col:
            df[f'{zone_label}(s)'] = pd.to_numeric(df[col], errors='coerce')
    center_col = _find_col(df, '中心区', '区域时间')
    if center_col:
        df['中心区域时间占比(%)'] = pd.to_numeric(df[center_col], errors='coerce') / 300 * 100
    return df


def _load_flux_water_maze_csv(csv_path) -> pd.DataFrame:
    """flux 水迷宫 csv (utf-8, 有 ID 列); rename ID→MetaboID 供 _compute_mwm_derived。"""
    raw = pd.read_csv(csv_path)  # utf-8, 列 Group/MouseID/Day_Num/Daily_Mean/ID
    raw['ID'] = raw['ID'].astype(str).str.strip().str.upper()
    raw['Daily_Mean'] = pd.to_numeric(raw['Daily_Mean'], errors='coerce')
    raw['Day_Num'] = pd.to_numeric(raw['Day_Num'], errors='coerce')
    return raw.rename(columns={'ID': 'MetaboID'})


def _load_flux() -> BehavioralDataset:
    """代谢流小鼠: Y迷宫 + 新物体 + 旷场 + 水迷宫 (全 derived)。"""
    gm = _load_flux_group_map()
    y_maze = _load_y_maze(_FLUX_DIR / 'Y迷宫数据_Y_Maze.xlsx', 'Sheet1', gm)
    nor = _load_nor(_FLUX_DIR / '新物体识别数据_NOR.xlsx', 'Sheet2', gm)
    open_field = _load_flux_open_field(_FLUX_DIR / '旷场数据_Open_Field.xlsx', 'Sheet1', gm)
    dist_csv = _FLUX_DIR / '水迷宫数据' / 'path_lenth_subset.csv'   # path length = distance
    lat_csv = _FLUX_DIR / '水迷宫数据' / 'latency_subset.csv'
    if dist_csv.exists() and lat_csv.exists():
        mwm = _compute_mwm_derived(_load_flux_water_maze_csv(dist_csv),
                                   _load_flux_water_maze_csv(lat_csv))
    else:
        print(f'flux: 未找到水迷宫 csv, 跳过 MWM')
        mwm = None
    return BehavioralDataset(nor=nor, open_field=open_field, y_maze=y_maze, mwm=mwm, group_map=gm)
