"""数据加载：89 个 GDC STAR 定量 TSV → 统一矩阵。

输出：
- log2_tpm : DataFrame，index=gene_id（带版本），columns=样本目录 UUID，值=log2(TPM+1)
- raw_counts: DataFrame，同形状，值=unstranded 原始计数（仅 QC 用）
- gene_info : DataFrame，index=gene_id，列 gene_name / gene_type
- summaries : dict[sample] -> {N_unmapped: ..., ...}
- sample_meta: dict[sample] -> {file_uuid, tsv_path, dir_path}
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg


def discover_tsvs() -> list[dict]:
    """返回每个样本的信息字典，按目录 UUID 排序，保证顺序稳定。"""
    samples = []
    for sub in sorted(cfg.RNA_DIR.iterdir()):
        if not sub.is_dir():
            continue
        tsvs = sorted(sub.glob("*.rna_seq.augmented_star_gene_counts.tsv"))
        if not tsvs:
            continue
        tsv = tsvs[0]
        samples.append({
            "sample": sub.name,
            "file_uuid": tsv.name.split(".rna_seq.")[0],
            "tsv_path": str(tsv),
            "dir_path": str(sub),
        })
    return samples


def _read_one(path: str):
    """读单个 TSV，跳过首行注释，分离 N_* 汇总行与基因行。"""
    df = pd.read_csv(path, sep="\t", skiprows=1, dtype={
        "gene_id": str, "gene_name": str, "gene_type": str,
    })
    mask_summary = df["gene_id"].astype(str).str.startswith(cfg.SUMMARY_PREFIX)
    summary_rows = df.loc[mask_summary, ["gene_id"] + [
        cfg.COL_COUNT, "stranded_first", "stranded_second"
    ]].copy()
    summary = {}
    for _, row in summary_rows.iterrows():
        summary[str(row["gene_id"])] = int(float(row[cfg.COL_COUNT]))
    genes = df.loc[~mask_summary].copy()
    genes = genes.set_index("gene_id")
    return genes, summary


def load_all() -> dict:
    samples = discover_tsvs()
    if not samples:
        raise FileNotFoundError(f"未在 {cfg.RNA_DIR} 找到 STAR 定量 TSV")

    gene_order = None
    gene_info = None
    log2_frames = []
    count_frames = []
    summaries = {}
    sample_meta = {}

    for s in samples:
        genes, summary = _read_one(s["tsv_path"])
        summaries[s["sample"]] = summary
        sample_meta[s["sample"]] = {
            "file_uuid": s["file_uuid"],
            "tsv_path": s["tsv_path"],
            "dir_path": s["dir_path"],
        }
        if gene_order is None:
            gene_order = genes.index.tolist()
            gene_info = genes[["gene_name", "gene_type"]].copy()
            gene_info.index.name = "gene_id"
        else:
            if genes.index.tolist() != gene_order:
                genes = genes.reindex(gene_order)
        tpm = pd.to_numeric(genes[cfg.COL_TPM], errors="coerce").fillna(0.0).astype(float)
        cnt = pd.to_numeric(genes[cfg.COL_COUNT], errors="coerce").fillna(0.0).astype(float)
        log2_frames.append(np.log2(tpm + 1.0).rename(s["sample"]))
        count_frames.append(cnt.rename(s["sample"]))

    log2_tpm = pd.concat(log2_frames, axis=1)
    raw_counts = pd.concat(count_frames, axis=1)
    log2_tpm.index.name = "gene_id"
    raw_counts.index.name = "gene_id"

    return {
        "log2_tpm": log2_tpm,
        "raw_counts": raw_counts,
        "gene_info": gene_info,
        "summaries": summaries,
        "sample_meta": sample_meta,
        "samples": [s["sample"] for s in samples],
    }


def name_to_ids(gene_info: pd.DataFrame) -> dict[str, list[str]]:
    """gene_name -> 该名称对应的所有 gene_id（处理重复 gene_name）。"""
    mapping: dict[str, list[str]] = {}
    for gid, row in gene_info.iterrows():
        name = row["gene_name"]
        if not isinstance(name, str) or not name:
            continue
        mapping.setdefault(name, []).append(str(gid))
    return mapping


def pick_representative_id(
    name: str, name_map: dict[str, list[str]], log2_tpm: pd.DataFrame, sample: str
) -> str | None:
    """一个 gene_name 多 ID 时，取该样本中表达最高的 ID 作为展示值。"""
    ids = name_map.get(name)
    if not ids:
        return None
    if len(ids) == 1:
        return ids[0]
    vals = log2_tpm.loc[[i for i in ids if i in log2_tpm.index], sample]
    if vals.empty:
        return None
    return str(vals.idxmax())
