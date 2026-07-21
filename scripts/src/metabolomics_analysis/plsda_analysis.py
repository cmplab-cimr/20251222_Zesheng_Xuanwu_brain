# -*- coding: utf-8 -*-
"""PLS-DA（偏最小二乘判别分析）

对两组代谢组学数据进行判别分析，输出 VIP 分数排序和得分散点图。
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler


def _impute_half_min(X: np.ndarray) -> np.ndarray:
    """用每列最小值的 0.5 倍填充 NaN"""
    X = X.copy()
    for j in range(X.shape[1]):
        col = X[:, j]
        nan_mask = np.isnan(col)
        if nan_mask.any():
            fill_val = col[~nan_mask].min() * 0.5 if (~nan_mask).any() else 0
            col[nan_mask] = fill_val
    return X


def _calculate_vip(pls: PLSRegression, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """计算 VIP (Variable Importance in Projection) 分数"""
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_
    p, h = w.shape

    # 每个成分解释的 Y 方差
    ss = np.array([(q[0, c] ** 2) * (t[:, c] @ t[:, c]) for c in range(h)])
    total_ss = ss.sum()
    if total_ss == 0:
        return np.zeros(p)

    w_sq_norm = (w ** 2).sum(axis=0)  # 每个成分的权重平方和
    vip = np.zeros(p)
    for j in range(p):
        vip[j] = np.sqrt(p * np.sum(ss * w[j] ** 2 / w_sq_norm) / total_ss)
    return vip


def _confidence_ellipse(x, y, ax, n_std=1.96, **kwargs):
    """绘制 95% 置信椭圆"""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width = 2 * n_std * np.sqrt(max(vals[0], 0))
    height = 2 * n_std * np.sqrt(max(vals[1], 0))
    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                  width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ell)


class PLSDAResult:
    """PLS-DA 分析结果"""

    def __init__(self):
        self.vip_scores: pd.Series | None = None
        self.coefficients: pd.Series | None = None
        self.scores: pd.DataFrame | None = None
        self.r2x: float = 0
        self.r2y: float = 0
        self.n_components: int = 2
        self.metabolite_info: pd.DataFrame | None = None
        self.group_names: tuple[str, str] = ('Group1', 'Group2')

    def get_ranked_metabolites(self) -> pd.DataFrame:
        """按 |VIP| 降序返回代谢物"""
        vip = self.vip_scores
        order = vip.abs().sort_values(ascending=False).index
        df = pd.DataFrame({'VIP': vip[order], '|VIP|': vip.abs()[order]})
        if self.coefficients is not None:
            df['coefficient'] = self.coefficients[order]
        if self.metabolite_info is not None:
            info_cols = [c for c in ['Compounds', '物质', 'Class I', 'Class II', 'kegg_map']
                         if c in self.metabolite_info.columns]
            df = df.join(self.metabolite_info[info_cols].loc[order], how='left')
        return df


def run_plsda(
        data: pd.DataFrame,
        group1_cols: list[str], group2_cols: list[str],
        group1_name: str = 'Group1', group2_name: str = 'Group2',
        n_components: int = 2,
        min_valid_ratio: float = 0.5,
) -> PLSDAResult:
    """运行 PLS-DA

    Args:
        data: 代谢物(Index) × 样本 丰度矩阵
        group1_cols: 组1 列名
        group2_cols: 组2 列名
        n_components: PLS 成分数
        min_valid_ratio: 每组中有效值最少比例

    Returns:
        PLSDAResult
    """
    all_cols = group1_cols + group2_cols

    # 过滤缺失值过多的代谢物
    valid_mask = (
        (data[group1_cols].notna().mean(axis=1) >= min_valid_ratio)
        & (data[group2_cols].notna().mean(axis=1) >= min_valid_ratio)
    )
    data_filtered = data.loc[valid_mask, all_cols]

    # 转置：样本(行) × 代谢物(列)，再做 log2 + 填充 + 标准化
    X_raw = data_filtered.T.to_numpy(dtype=float)
    X = StandardScaler().fit_transform(_impute_half_min(np.log2(X_raw + 1)))

    # Y 编码
    Y = np.array([0] * len(group1_cols) + [1] * len(group2_cols), dtype=float).reshape(-1, 1)

    n_comp = min(n_components, X.shape[0] - 1, X.shape[1])
    pls = PLSRegression(n_components=n_comp)
    pls.fit(X, Y)

    # 结果
    result = PLSDAResult()
    result.n_components = n_comp
    metabolite_index = data_filtered.index
    result.vip_scores = pd.Series(
        _calculate_vip(pls, X, Y), index=metabolite_index, name='VIP')
    result.coefficients = pd.Series(
        pls.coef_.flatten(), index=metabolite_index, name='coefficient')

    scores_df = pd.DataFrame(
        pls.x_scores_, index=all_cols,
        columns=[f'Component {i + 1}' for i in range(n_comp)])
    scores_df['group'] = [group1_name] * len(group1_cols) + [group2_name] * len(group2_cols)
    result.scores = scores_df

    # R²X
    X_hat = pls.x_scores_ @ pls.x_loadings_.T
    result.r2x = 1 - np.sum((X - X_hat) ** 2) / np.sum((X - X.mean(axis=0)) ** 2)
    # R²Y
    result.r2y = float(pls.score(X, Y))
    result.group_names = (group1_name, group2_name)

    return result


def plot_plsda(
        result: PLSDAResult,
        output_path: str | Path,
        colors: tuple[str, str] = ('#E64B35', '#4DBBD5'),
):
    """绘制 PLSDA 得分散点图 + 95% 置信椭圆"""
    scores = result.scores
    if scores is None or result.n_components < 2:
        print('成分数不足 2，无法绘图')
        return

    fig, ax = plt.subplots(figsize=(6, 5.5))

    g1_name, g2_name = result.group_names
    for name, color in zip([g1_name, g2_name], colors):
        mask = scores['group'] == name
        x = scores.loc[mask, 'Component 1'].values
        y = scores.loc[mask, 'Component 2'].values
        ax.scatter(x, y, c=color, s=60, edgecolors='white', linewidths=0.8,
                   label=name, zorder=3)
        _confidence_ellipse(x, y, ax, n_std=1.96,
                            facecolor=color, alpha=0.15, edgecolor=color,
                            linewidth=1.5, linestyle='--')

    ax.set_xlabel('Component 1', fontsize=12)
    ax.set_ylabel('Component 2', fontsize=12)
    ax.set_title('PLS-DA', fontsize=13, fontweight='bold')
    ax.legend(frameon=True, fontsize=10, loc='best')

    stats_text = f'R²X = {result.r2x:.3f}\nR²Y = {result.r2y:.3f}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存图片: {path}')


def export_plsda(result: PLSDAResult, output_path: str | Path):
    """导出 VIP 排序结果到 xlsx"""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ranked = result.get_ranked_metabolites()
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        ranked.to_excel(writer, sheet_name='VIP排名')
        if result.scores is not None:
            result.scores.to_excel(writer, sheet_name='样本得分')
    print(f'已导出: {path}')


# ---- OPLS-DA ----

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

TISSUE_NAMES_ZH = {
    'cortex': '皮层', 'hippocampus': '海马', 'hypothalamus': '下丘脑',
    'liver': '肝脏', 'kidney': '肾脏', 'fat': '脂肪',
    'muscle': '肌肉', 'ileum': '回肠', 'serum': '血浆',
    'gastrocnemius': '腓肠肌', 'quadriceps': '股四头肌',
}


class OPLSDAResult:
    """OPLS-DA 分析结果"""
    def __init__(self):
        self.scores: pd.DataFrame | None = None
        self.r2x: float = 0
        self.r2y: float = 0
        self.q2: float = 0
        self.vip_scores: pd.Series | None = None
        self.tissue: str = ''
        self.group_names: tuple[str, str] = ('Group1', 'Group2')
        self.metabolite_info: pd.DataFrame | None = None
        # S-plot 数据
        self.p_loadings: pd.Series | None = None  # p[1]: 预测载荷 (covariance)
        self.p_corr: pd.Series | None = None      # p(corr)[1]: 相关系数载荷

    def get_ranked_metabolites(self) -> pd.DataFrame:
        """按 |VIP| 降序返回代谢物"""
        vip = self.vip_scores
        if vip is None:
            return pd.DataFrame()
        order = vip.abs().sort_values(ascending=False).index
        df = pd.DataFrame({'VIP': vip[order], '|VIP|': vip.abs()[order]})
        if self.metabolite_info is not None:
            info_cols = [c for c in ['Compounds', '物质', 'Class I', 'Class II', 'kegg_map']
                         if c in self.metabolite_info.columns]
            df = df.join(self.metabolite_info[info_cols].loc[order], how='left')
        return df


def _oplsda_core(X: np.ndarray, y: np.ndarray):
    """OPLS-DA 核心算法: 1 预测成分 + 1 正交成分"""
    n, p = X.shape
    y = y.ravel()

    # 预测权重 (与 Y 最大协变方向)
    w_pred = X.T @ y
    w_norm = np.linalg.norm(w_pred)
    if w_norm < 1e-10:
        return np.zeros(n), np.zeros((n, 1)), 0., 0., np.zeros(p), np.zeros(p)
    w_pred = w_pred / w_norm

    # 预测得分
    t_pred = X @ w_pred
    t_ss = t_pred @ t_pred
    if t_ss < 1e-10:
        return np.zeros(n), np.zeros((n, 1)), 0., 0., w_pred, np.zeros(p)

    # 预测载荷
    p_pred = (X.T @ t_pred) / t_ss

    # 正交权重 (Gram-Schmidt: 去除与预测权重平行的分量)
    w_ortho = p_pred - w_pred * (w_pred @ p_pred)
    w_ortho_norm = np.linalg.norm(w_ortho)
    has_ortho = w_ortho_norm > 1e-10

    if has_ortho:
        w_ortho /= w_ortho_norm
        t_ortho = (X @ w_ortho).reshape(-1, 1)
        p_ortho = (X.T @ t_ortho.ravel()) / (t_ortho.ravel() @ t_ortho.ravel() + 1e-12)
    else:
        t_ortho = np.zeros((n, 1))
        p_ortho = np.zeros(p)

    # R2Y
    b = (y @ t_pred) / t_ss
    Y_hat = t_pred * b
    ss_y = y @ y
    R2Y = 1 - np.sum((y - Y_hat) ** 2) / ss_y if ss_y > 0 else 0.

    # R2X
    X_hat = np.outer(t_pred, p_pred)
    if has_ortho:
        X_hat += np.outer(t_ortho.ravel(), p_ortho)
    ss_x = np.sum(X ** 2)
    R2X = 1 - np.sum((X - X_hat) ** 2) / ss_x if ss_x > 0 else 0.

    return t_pred, t_ortho, R2X, R2Y, w_pred, p_pred


def _oplsda_q2(X: np.ndarray, y: np.ndarray) -> float:
    """交叉验证 Q² (PLS1 近似预测成分)"""
    from sklearn.model_selection import KFold, LeaveOneOut, cross_val_predict

    n = X.shape[0]
    n_splits = min(5, n // 2) if n >= 10 else n
    cv = LeaveOneOut() if n < 10 else KFold(n_splits=n_splits, shuffle=True, random_state=42)

    pls = PLSRegression(n_components=1)
    Y_pred = cross_val_predict(pls, X, y.reshape(-1, 1), cv=cv)

    ss_res = np.sum((y - Y_pred.ravel()) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.


def run_oplsda(
        data: pd.DataFrame,
        group1_cols: list[str], group2_cols: list[str],
        group1_name: str = 'Group1', group2_name: str = 'Group2',
        min_valid_ratio: float = 0.5,
) -> OPLSDAResult:
    """运行 OPLS-DA"""
    all_cols = group1_cols + group2_cols

    valid_mask = (
        (data[group1_cols].notna().mean(axis=1) >= min_valid_ratio)
        & (data[group2_cols].notna().mean(axis=1) >= min_valid_ratio)
    )
    data_filtered = data.loc[valid_mask, all_cols]

    X_raw = data_filtered.T.to_numpy(dtype=float)
    X = StandardScaler().fit_transform(_impute_half_min(np.log2(X_raw + 1)))
    Y = np.array([0.] * len(group1_cols) + [1.] * len(group2_cols))
    y_c = Y - Y.mean()

    t_pred, t_ortho, r2x, r2y, w_pred, p_pred = _oplsda_core(X, y_c)
    q2 = _oplsda_q2(X, y_c)

    # VIP (单预测成分: VIP_j = sqrt(p * w_j² / Σw_k²))
    p_feat = X.shape[1]
    w_sq = w_pred @ w_pred
    vip = np.sqrt(p_feat * w_pred ** 2 / w_sq) if w_sq > 0 else np.zeros(p_feat)

    # S-plot: p(corr)[1] = loading * std(t) / std(x_j)
    t_std = t_pred.std(ddof=1) if t_pred.std(ddof=1) > 0 else 1.0
    x_std = X.std(axis=0, ddof=1)
    x_std[x_std == 0] = 1.0
    p_corr = p_pred * t_std / x_std

    result = OPLSDAResult()
    result.r2x, result.r2y, result.q2 = r2x, r2y, q2
    result.group_names = (group1_name, group2_name)
    result.vip_scores = pd.Series(vip, index=data_filtered.index, name='VIP')
    result.p_loadings = pd.Series(p_pred, index=data_filtered.index, name='p[1]')
    result.p_corr = pd.Series(p_corr, index=data_filtered.index, name='p(corr)[1]')
    result.scores = pd.DataFrame({
        'Predictive': t_pred,
        'Orthogonal': t_ortho.ravel(),
        'group': [group1_name] * len(group1_cols) + [group2_name] * len(group2_cols),
    }, index=all_cols)

    return result


def permutation_test_oplsda(
        data: pd.DataFrame,
        group1_cols: list[str], group2_cols: list[str],
        n_permutations: int = 200,
        min_valid_ratio: float = 0.5,
) -> dict:
    """OPLS-DA 置换检验"""
    all_cols = group1_cols + group2_cols

    valid_mask = (
        (data[group1_cols].notna().mean(axis=1) >= min_valid_ratio)
        & (data[group2_cols].notna().mean(axis=1) >= min_valid_ratio)
    )
    data_filtered = data.loc[valid_mask, all_cols]

    X_raw = data_filtered.T.to_numpy(dtype=float)
    X = StandardScaler().fit_transform(_impute_half_min(np.log2(X_raw + 1)))
    Y = np.array([0.] * len(group1_cols) + [1.] * len(group2_cols))
    y_c = Y - Y.mean()

    _, _, _, orig_r2y, _, _ = _oplsda_core(X, y_c)
    orig_q2 = _oplsda_q2(X, y_c)

    rng = np.random.default_rng(42)
    perm_r2y, perm_q2 = np.zeros(n_permutations), np.zeros(n_permutations)

    for i in range(n_permutations):
        y_perm = rng.permutation(y_c)
        _, _, _, r2y, _, _ = _oplsda_core(X, y_perm)
        perm_r2y[i] = r2y
        perm_q2[i] = _oplsda_q2(X, y_perm)

    p_r2y = (np.sum(perm_r2y >= orig_r2y) + 1) / (n_permutations + 1)
    p_q2 = (np.sum(perm_q2 >= orig_q2) + 1) / (n_permutations + 1)

    return {
        'original_r2y': orig_r2y, 'original_q2': orig_q2,
        'perm_r2y': perm_r2y, 'perm_q2': perm_q2,
        'p_r2y': p_r2y, 'p_q2': p_q2,
        'n_permutations': n_permutations,
    }


def plot_oplsda(
        result: OPLSDAResult,
        output_path: str | Path,
        colors: tuple[str, str] = ('#E64B35', '#4DBBD5'),
):
    """单个 OPLS-DA 得分散点图"""
    scores = result.scores
    if scores is None:
        return

    fig, ax = plt.subplots(figsize=(6, 5.5))
    g1_name, g2_name = result.group_names

    for name, color in zip([g1_name, g2_name], colors):
        mask = scores['group'] == name
        x = scores.loc[mask, 'Predictive'].values
        y = scores.loc[mask, 'Orthogonal'].values
        ax.scatter(x, y, c=color, s=60, edgecolors='white', linewidths=0.8,
                   label=name, zorder=3)
        _confidence_ellipse(x, y, ax, n_std=1.96,
                            facecolor=color, alpha=0.15, edgecolor=color,
                            linewidth=1.5, linestyle='--')

    ax.axhline(0, color='gray', lw=0.5, alpha=0.5)
    ax.axvline(0, color='gray', lw=0.5, alpha=0.5)
    ax.set_xlabel('Predictive component (组间)', fontsize=12)
    ax.set_ylabel('Orthogonal component (组内)', fontsize=12)
    tissue_label = TISSUE_NAMES_ZH.get(result.tissue, result.tissue)
    ax.set_title(f'OPLS-DA — {tissue_label}', fontsize=13, fontweight='bold')
    ax.legend(frameon=True, fontsize=10, loc='best')

    stats_text = f'R²Y = {result.r2y:.3f}\nQ² = {result.q2:.3f}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存图片: {path}')


def plot_oplsda_grid(
        results: list[OPLSDAResult],
        output_path: str | Path,
        ncols: int = 3,
        colors: tuple[str, str] = ('#E64B35', '#4DBBD5'),
):
    """OPLS-DA small multiples 网格图"""
    n = len(results)
    if n == 0:
        return

    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows))
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    axes_flat = axes.flatten()

    for idx, result in enumerate(results):
        ax = axes_flat[idx]
        scores = result.scores
        if scores is None:
            ax.set_visible(False)
            continue

        g1_name, g2_name = result.group_names
        for name, color in zip([g1_name, g2_name], colors):
            mask = scores['group'] == name
            x = scores.loc[mask, 'Predictive'].values
            y = scores.loc[mask, 'Orthogonal'].values
            ax.scatter(x, y, c=color, s=40, edgecolors='white', linewidths=0.5,
                       label=name, zorder=3)
            _confidence_ellipse(x, y, ax, n_std=1.96,
                                facecolor=color, alpha=0.15, edgecolor=color,
                                linewidth=1.0, linestyle='--')

        ax.axhline(0, color='gray', lw=0.5, alpha=0.5)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.5)

        tissue_label = TISSUE_NAMES_ZH.get(result.tissue, result.tissue)
        ax.set_title(tissue_label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Predictive', fontsize=9)
        ax.set_ylabel('Orthogonal', fontsize=9)

        stats_text = f'R²Y={result.r2y:.2f} Q²={result.q2:.2f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                fontsize=7, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.5))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=7)

    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=10,
               frameon=True, bbox_to_anchor=(0.5, 0.01))

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存网格图: {path}')


def plot_permutation_test(
        perm_result: dict,
        tissue: str = '',
        output_path: str | Path | None = None,
):
    """置换检验图 (R²Y vs Q² 散点 + 分布直方图)"""
    orig_r2y = perm_result['original_r2y']
    orig_q2 = perm_result['original_q2']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # R²Y vs Q² 散点 (SIMCA 风格)
    ax = axes[0]
    ax.scatter(perm_result['perm_r2y'], perm_result['perm_q2'],
               c='gray', s=15, alpha=0.5, label='Permuted')
    ax.scatter([orig_r2y], [orig_q2], c='#E64B35', s=80, marker='*',
               zorder=5, label='Original')
    if len(perm_result['perm_r2y']) > 1:
        z = np.polyfit(perm_result['perm_r2y'], perm_result['perm_q2'], 1)
        x_line = np.linspace(perm_result['perm_r2y'].min(), orig_r2y * 1.1, 100)
        ax.plot(x_line, np.poly1d(z)(x_line), 'k--', lw=1, alpha=0.7)
        intercept_text = f'Intercept Q²={z[1]:.3f}'
    else:
        intercept_text = ''

    ax.set_xlabel('R²Y', fontsize=10)
    ax.set_ylabel('Q²', fontsize=10)
    title = f'{TISSUE_NAMES_ZH.get(tissue, tissue)} — Permutation Test' if tissue else 'Permutation Test'
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 分布直方图
    ax = axes[1]
    ax.hist(perm_result['perm_r2y'], bins=20, color='#4DBBD5', alpha=0.6, label='R²Y')
    ax.hist(perm_result['perm_q2'], bins=20, color='#E64B35', alpha=0.6, label='Q²')
    ax.axvline(orig_r2y, color='#4DBBD5', lw=2, ls='--')
    ax.axvline(orig_q2, color='#E64B35', lw=2, ls='--')

    info = f'p(R²Y)={perm_result["p_r2y"]:.4f}\np(Q²)={perm_result["p_q2"]:.4f}'
    if intercept_text:
        info += f'\n{intercept_text}'
    ax.text(0.98, 0.98, info, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))
    ax.set_xlabel('Value', fontsize=10)
    ax.set_ylabel('Count', fontsize=10)
    ax.set_title(f'Permutation Distribution ({perm_result["n_permutations"]}x)',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if output_path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f'已保存置换检验图: {path}')
    plt.close(fig)


def plot_permutation_grid(
        perm_results: dict[str, dict],
        output_path: str | Path,
        ncols: int = 3,
):
    """置换检验网格图（补充材料）"""
    n = len(perm_results)
    if n == 0:
        return

    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    axes_flat = axes.flatten()

    for idx, (tissue, pr) in enumerate(perm_results.items()):
        ax = axes_flat[idx]
        ax.scatter(pr['perm_r2y'], pr['perm_q2'], c='gray', s=10, alpha=0.4)
        ax.scatter([pr['original_r2y']], [pr['original_q2']],
                   c='#E64B35', s=60, marker='*', zorder=5)
        if len(pr['perm_r2y']) > 1:
            z = np.polyfit(pr['perm_r2y'], pr['perm_q2'], 1)
            x_line = np.linspace(pr['perm_r2y'].min(), pr['original_r2y'] * 1.1, 50)
            ax.plot(x_line, np.poly1d(z)(x_line), 'k--', lw=0.8, alpha=0.7)

        ax.set_title(TISSUE_NAMES_ZH.get(tissue, tissue), fontsize=10, fontweight='bold')
        ax.set_xlabel('R²Y', fontsize=8)
        ax.set_ylabel('Q²', fontsize=8)
        ax.text(0.02, 0.98, f"p={pr['p_q2']:.3f}", transform=ax.transAxes, fontsize=7,
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=7)

    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存置换检验网格图: {path}')


def run_combined_pca(
        tissue_data: dict[str, tuple[pd.DataFrame, list[str], list[str]]],
        group1_name: str, group2_name: str,
        output_path: str | Path,
):
    """全组织合并 PCA"""
    from sklearn.decomposition import PCA

    if len(tissue_data) < 2:
        print('有效组织不足 2，跳过合并 PCA')
        return

    # 共有代谢物
    common_idx = None
    for data, _, _ in tissue_data.values():
        common_idx = data.index if common_idx is None else common_idx.intersection(data.index)
    if len(common_idx) < 5:
        print(f'共有代谢物仅 {len(common_idx)} 个，跳过合并 PCA')
        return

    # 构建矩阵: 样本(行) × 代谢物(列)
    X_list, meta = [], []
    for t, (data, g1_cols, g2_cols) in tissue_data.items():
        data_c = data.loc[common_idx]
        for col in g1_cols:
            if col in data_c.columns:
                X_list.append(data_c[col].values)
                meta.append({'tissue': t, 'group': group1_name})
        for col in g2_cols:
            if col in data_c.columns:
                X_list.append(data_c[col].values)
                meta.append({'tissue': t, 'group': group2_name})

    X = np.array(X_list, dtype=float)
    X = StandardScaler().fit_transform(_impute_half_min(np.log2(X + 1)))

    pca = PCA(n_components=2)
    scores = pca.fit_transform(X)

    meta_df = pd.DataFrame(meta)
    meta_df['PC1'], meta_df['PC2'] = scores[:, 0], scores[:, 1]

    # 绘图: 颜色=组织, 形状=分组
    tissue_list = sorted(meta_df['tissue'].unique())
    n_tissues = len(tissue_list)
    tissue_colors = plt.cm.tab10(np.linspace(0, 1, max(n_tissues, 10)))[:n_tissues]
    t2c = dict(zip(tissue_list, tissue_colors))

    fig, ax = plt.subplots(figsize=(8, 7))
    for t in tissue_list:
        mask_t = meta_df['tissue'] == t
        for grp, marker in zip([group1_name, group2_name], ['o', '^']):
            mask = mask_t & (meta_df['group'] == grp)
            if not mask.any():
                continue
            label = f'{TISSUE_NAMES_ZH.get(t, t)} ({grp})'
            ax.scatter(meta_df.loc[mask, 'PC1'], meta_df.loc[mask, 'PC2'],
                       c=[t2c[t]], marker=marker, s=50,
                       edgecolors='white', linewidths=0.5,
                       label=label, zorder=3)

    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=12)
    ax.set_title('Combined PCA — All Tissues', fontsize=13, fontweight='bold')
    ax.legend(fontsize=7, ncol=2, loc='best', frameon=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存合并 PCA 图: {path}')


def export_oplsda(result: OPLSDAResult, output_path: str | Path):
    """导出 OPLS-DA 结果"""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ranked = result.get_ranked_metabolites()
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        ranked.to_excel(writer, sheet_name='VIP排名')
        if result.scores is not None:
            result.scores.to_excel(writer, sheet_name='样本得分')
    print(f'已导出: {path}')
