"""队列参照基线：留一法（leave-one-out）Percentile 与 Z-score。

对样本 i 的基因 g：
  percentile = 100 * (count(x_other <= x_i) + 1) / N
  Z          = (x_i - mean_others) / std_others   (std 用 ddof=1)

缓存为 outputs/cache/baseline.npz，避免重复计算。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg


def compute_baseline(data: dict, force: bool = False) -> dict:
    cache_path = cfg.CACHE_DIR / "baseline.npz"
    cfg.ensure_output_dirs()
    samples = data["samples"]
    X = data["log2_tpm"].values  # (genes, samples)
    genes = data["log2_tpm"].index

    if cache_path.exists() and not force:
        cached = np.load(cache_path)
        cached_samples = [str(s) for s in cached["samples"]]
        if (cached_samples == samples
                and cached["percentile"].shape == X.shape):
            return {
                "percentile": pd.DataFrame(
                    cached["percentile"], index=genes, columns=samples),
                "zscore": pd.DataFrame(
                    cached["zscore"], index=genes, columns=samples),
            }

    n = X.shape[1]
    pct = np.zeros_like(X, dtype=np.float32)
    zsc = np.zeros_like(X, dtype=np.float32)
    for i in range(n):
        xi = X[:, i:i + 1]
        others = np.concatenate([X[:, :i], X[:, i + 1:]], axis=1)
        count = np.count_nonzero(others <= xi, axis=1)
        pct[:, i] = 100.0 * (count + 1.0) / n
        mean = others.mean(axis=1)
        std = others.std(axis=1, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (X[:, i] - mean) / std
        z[~np.isfinite(z)] = 0.0
        zsc[:, i] = z

    # 样本名用定长 Unicode 数组保存，避免 object 数组需 allow_pickle 才能读回
    np.savez_compressed(
        cache_path,
        percentile=pct,
        zscore=zsc,
        samples=np.array(samples, dtype="U64"),
    )
    return {
        "percentile": pd.DataFrame(pct, index=genes, columns=samples),
        "zscore": pd.DataFrame(zsc, index=genes, columns=samples),
    }


def direction(percentile: float, zscore: float | None = None) -> str:
    """纯百分位驱动的方向判定：≥P80 为 high，≤P20 为 low，其余 normal。

    Z-score 不再参与判定（仅保留在 JSON 供审计），避免偏态分布下 Z 与 P 的直觉冲突。
    """
    if percentile >= cfg.P_HIGH:
        return "high"
    if percentile <= cfg.P_LOW:
        return "low"
    return "normal"
