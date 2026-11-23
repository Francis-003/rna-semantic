"""ssGSEA 通路打分（自实现，不依赖 gseapy）。

对每个样本：按表达降序排列基因，对目标基因集计算运行累加得分（ES），
权重采用 |rank statistic|^p（p=1）标准化；再将该 ES 在队列中取百分位。
过滤：仅 protein_coding，且排除 MT- / RPL / RPS / 假基因等非编码背景。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg
from .load import name_to_ids, pick_representative_id


def read_gmt(path) -> dict[str, dict]:
    """读取 GMT：name, description, genes... -> {name: {desc, genes}}"""
    out = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            out[parts[0]] = {"description": parts[1], "genes": parts[2:]}
    return out


def read_gene_list(path) -> list[str]:
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            out.append(s)
    return out


def _usable_gene_mask(gene_info: pd.DataFrame) -> np.ndarray:
    """保留 protein_coding 且非线粒体/核糖体/假基因的基因。"""
    gtype = gene_info["gene_type"].astype(str)
    gname = gene_info["gene_name"].astype(str)
    mask = (gtype == "protein_coding")
    mask &= ~gname.str.startswith("MT-")
    mask &= ~gname.str.match(r"^RPL|^RPS")
    return mask.values


def ssGSEA_scores(data: dict) -> pd.DataFrame:
    """返回 DataFrame：index=通路名，columns=样本，值=该样本通路 ES。"""
    gmt = read_gmt(cfg.GENE_SETS_DIR / "escc_pathways.gmt")
    log2_tpm = data["log2_tpm"]
    gene_info = data["gene_info"]
    samples = data["samples"]
    name_map = name_to_ids(gene_info)
    mask = _usable_gene_mask(gene_info)
    usable = gene_info.index[mask]
    values = log2_tpm.loc[usable, samples]

    # 每个通路映射到代表 gene_id 集合（每个 gene_name 取其在基因集中的最高表达 ID）
    set_ids = {}
    for pname, info in gmt.items():
        ids = set()
        for gname in info["genes"]:
            gid = pick_representative_id(gname, name_map, log2_tpm, samples[0])
            if gid is not None and gid in set(usable):
                ids.add(gid)
        set_ids[pname] = ids

    scores = pd.DataFrame(index=list(gmt.keys()), columns=samples, dtype=float)
    # 预计算每个样本的排序与秩统计量
    for s in samples:
        v = values[s].values
        order = np.argsort(-v, kind="stable")  # 降序
        sorted_v = v[order]
        sorted_ids = usable.values[order]
        N = len(sorted_ids)
        absv = np.abs(sorted_v)
        total_abs = absv.sum()
        for pname, ids in set_ids.items():
            k = len(ids)
            if k == 0 or k >= N:
                scores.loc[pname, s] = 0.0
                continue
            hit = np.isin(sorted_ids, list(ids))
            # 命中基因权重 = |表达|^p 归一化（正向累加）；
            # 未命中基因统一施加负向惩罚 -1/(N-k)，使累积曲线随命中/未命中上下波动，
            # 这样 ES 才能反映通路在排序顶端/底端的富集程度。
            w_hit = absv[hit]
            sw = w_hit.sum()
            if sw <= 0:
                w_hit = np.ones(k)
                sw = k
            increment = np.zeros(N)
            increment[hit] = w_hit / sw
            increment[~hit] = -1.0 / (N - k)
            running = np.cumsum(increment)
            mx = running.max()
            mn = running.min()
            # 富集在顶端取最大正值，富集在底端取最大负值
            es = mx if mx >= -mn else mn
            scores.loc[pname, s] = float(es)
    return scores.astype(float)
