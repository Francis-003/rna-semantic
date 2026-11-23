"""质控：线粒体 reads 占比、比对汇总、低质量标注。"""
from __future__ import annotations

from pathlib import Path

from . import config as cfg
from .load import name_to_ids


def compute_qc(data: dict) -> dict[str, dict]:
    log2_tpm = data["log2_tpm"]
    raw_counts = data["raw_counts"]
    gene_info = data["gene_info"]
    summaries = data["summaries"]
    sample_meta = data["sample_meta"]
    samples = data["samples"]

    name_map = name_to_ids(gene_info)
    mt_ids = [
        gid for name, ids in name_map.items()
        if name.startswith("MT-") for gid in ids
    ]
    # 以防不同 GENCODE 版本命名差异，再按 gene_name 前缀直接兜底
    mt_ids = sorted(set(mt_ids) | set(
        gene_info.index[gene_info["gene_name"].astype(str).str.startswith("MT-")]
    ))
    mt_ids = [g for g in mt_ids if g in raw_counts.index]

    qc = {}
    for s in samples:
        total = float(raw_counts[s].sum())
        mt_reads = float(raw_counts.loc[mt_ids, s].sum())
        mt_pct = 100.0 * mt_reads / total if total > 0 else 0.0
        sm = summaries.get(s, {})
        n_unmap = sm.get("N_unmapped", 0)
        n_multi = sm.get("N_multimapping", 0)
        n_nofeature = sm.get("N_noFeature", 0)
        n_ambig = sm.get("N_ambiguous", 0)
        denominator = n_unmap + n_multi + n_nofeature + n_ambig + total
        unmapped_pct = 100.0 * n_unmap / denominator if denominator > 0 else 0.0
        low_quality = bool(
            (Path(sample_meta[s]["dir_path"]) / "annotations.txt").exists()
        )
        verdict = "需谨慎解读（GDC 低质量标注）" if low_quality else "数据可信"
        if mt_pct > 20:
            verdict += "；MT% 偏高，提示 RNA 降解可能"
        qc[s] = {
            "sample": s,
            "mt_pct": round(mt_pct, 2),
            "mt_reads": int(mt_reads),
            "total_reads": int(total),
            "unmapped_pct": round(unmapped_pct, 2),
            "N_unmapped": n_unmap,
            "N_multimapping": n_multi,
            "N_noFeature": n_nofeature,
            "N_ambiguous": n_ambig,
            "low_quality": low_quality,
            "qc_verdict": verdict,
        }
    return qc
