"""特征聚合：把 QC / 队列分位 / 通路打分 / 基因集 组织为结构化事实 JSON。

维度对应方案 5 段（放弃配对维度）：
  ① 样本质控与数据可信度
  ② 免疫微环境与免疫治疗提示
  ③ 增殖与核心驱动通路
  ④ 潜在靶向 / ADC 靶点线索
  ⑤ 临床解读边界（固定）
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg
from .baseline import direction
from .load import name_to_ids, pick_representative_id
from .pathways import read_gene_list, ssGSEA_scores

# 驱动基因按功能类别（与 gene_sets/drivers.txt 对应）
DRIVER_GROUPS = {
    "细胞周期与增殖": ["CCND1", "CDK4", "CDK6", "MKI67", "E2F1"],
    "NRF2 氧化应激": ["NFE2L2", "NQO1", "GCLC", "GCLM"],
    "鳞状分化与干性": ["SOX2", "TP63", "KRT17"],
    "EMT 与侵袭": ["VIM", "SNAI1", "TGFB1"],
}

PATHWAY_LABELS = {
    "G2M_E2F": "G2M/E2F 细胞周期",
    "NRF2_OXIDATIVE": "NRF2 氧化应激",
    "SQUAMOUS_STEMNESS": "鳞状分化与干性",
    "EMT_INVASION": "EMT 与侵袭",
    "IFNG_GEP": "IFN-γ / T 细胞炎症",
    "STROMAL": "基质/成纤维细胞",
}

# 仅对拷贝数驱动基因给"提示扩增/过表达"的审慎措辞
CNV_AMPLIFIABLE = {"CCND1", "SOX2", "EGFR", "ERBB2"}


def _gene_record(
    gene_name: str, sample: str, name_map: dict, log2_tpm: pd.DataFrame,
    pct: pd.DataFrame, zsc: pd.DataFrame,
) -> dict | None:
    gid = pick_representative_id(gene_name, name_map, log2_tpm, sample)
    if gid is None or gid not in pct.index:
        return None
    p = float(pct.loc[gid, sample])
    z = float(zsc.loc[gid, sample])
    rec = {
        "gene": gene_name,
        "gene_id": gid,
        "log2_tpm": round(float(log2_tpm.loc[gid, sample]), 3),
        "percentile": round(p, 1),
        "zscore": round(z, 2),
        "direction": direction(p, z),
    }
    if gene_name in CNV_AMPLIFIABLE and rec["direction"] == "high":
        rec["amp_hint"] = True  # 只允许"提示可能存在扩增/过表达"
    return rec


def _score_percentile(score_row: pd.Series, sample: str) -> float:
    """对通路打分矩阵的一行（全样本）做留一法百分位。"""
    vals = score_row.dropna().astype(float)
    x = float(score_row.get(sample, np.nan))
    if not np.isfinite(x) or vals.empty:
        return 50.0
    others = vals.drop(index=sample) if sample in vals.index else vals
    count = int(np.count_nonzero(others.values <= x))
    return round(100.0 * (count + 1.0) / len(vals), 1)


def build_facts(
    data: dict, baseline: dict, qc: dict, sample: str,
    ss: pd.DataFrame | None = None,
) -> dict:
    log2_tpm = data["log2_tpm"]
    gene_info = data["gene_info"]
    pct = baseline["percentile"]
    zsc = baseline["zscore"]
    samples = data["samples"]
    name_map = name_to_ids(gene_info)

    def gene_group(names: list[str]) -> list[dict]:
        out = []
        for g in names:
            rec = _gene_record(g, sample, name_map, log2_tpm, pct, zsc)
            if rec is not None:
                out.append(rec)
        return out

    # ── ① 质控 ──
    q = qc[sample]
    qc_fact = {
        "mt_pct": q["mt_pct"],
        "unmapped_pct": q["unmapped_pct"],
        "low_quality": q["low_quality"],
        "verdict": q["qc_verdict"],
        "note": "基于线粒体转录本占比与比对汇总评估数据可信度。",
    }

    # ── 通路打分与分位 ──
    if ss is None:
        ss = ssGSEA_scores(data)
    pathway_facts = []
    for pname, info in ss.iterrows():
        pathway_facts.append({
            "name": pname,
            "label": PATHWAY_LABELS.get(pname, pname),
            "es": round(float(info.get(sample, 0.0)), 4),
            "percentile": _score_percentile(info, sample),
        })
    pathway_map = {p["name"]: p for p in pathway_facts}

    # ── ② 免疫 ──
    checkpoints = read_gene_list(cfg.GENE_SETS_DIR / "checkpoints.txt")
    immune = {
        "checkpoints": gene_group(checkpoints),
        "ifng_gep": pathway_map.get("IFNG_GEP"),
        "stromal": pathway_map.get("STROMAL"),
    }
    # 免疫表型：IFN-γ 炎症 × CD8A × 耗竭检查点聚合 × 基质，细粒度判定
    ifng_p = immune["ifng_gep"]["percentile"] if immune["ifng_gep"] else 50.0
    stromal_p = immune["stromal"]["percentile"] if immune["stromal"] else 50.0
    cd8 = _gene_record("CD8A", sample, name_map, log2_tpm, pct, zsc)
    immune["cd8a"] = cd8
    cd8_p = cd8["percentile"] if cd8 else 50.0

    # 耗竭/共抑制检查点聚合分位 = 7 个检查点分位的中位数（排除 CD274）
    exhaust = [
        r["percentile"] for r in immune["checkpoints"]
        if r["gene"] in cfg.EXHAUSTION_CHECKPOINTS
    ]
    exhaust_p = float(np.median(exhaust)) if exhaust else 50.0
    immune["exhaustion_percentile"] = round(exhaust_p, 1)

    # 炎症轴：IFN-γ 与 CD8A 取较高者判定热倾向
    inflamed = max(ifng_p, cd8_p) >= cfg.P_HIGH
    cold = (ifng_p <= cfg.P_LOW or cd8_p <= cfg.P_LOW) and exhaust_p <= cfg.P_LOW
    if inflamed:
        base_label = "免疫炎性热"
    elif cold:
        base_label = "免疫冷排斥"
    else:
        base_label = "免疫中间"

    # 附加修饰：低耗竭 + 低基质 → "中间偏冷 / 低耗竭、低间质"
    modifiers = []
    if exhaust_p <= cfg.P_LOW:
        modifiers.append("低耗竭")
    elif exhaust_p >= cfg.P_HIGH:
        modifiers.append("高耗竭")
    if stromal_p <= cfg.P_LOW:
        modifiers.append("低间质")
    elif stromal_p >= cfg.P_HIGH:
        modifiers.append("高间质")
    if base_label == "免疫中间" and modifiers:
        immune["phenotype"] = f"中间偏冷 / {'、'.join(modifiers)}"
    elif modifiers:
        immune["phenotype"] = f"{base_label} / {'、'.join(modifiers)}"
    else:
        immune["phenotype"] = base_label
    # 触发了耗竭标签时，记录 Top2–3 代表性检查点供正文引用（标签与证据不脱节）
    if exhaust_p <= cfg.P_LOW:
        ex_sorted = sorted(
            (r for r in immune["checkpoints"] if r["gene"] in cfg.EXHAUSTION_CHECKPOINTS),
            key=lambda r: r["percentile"],
        )[:3]
        immune["exhaustion_evidence"] = [
            {"gene": r["gene"], "percentile": r["percentile"]} for r in ex_sorted
        ]
        immune["exhaustion_direction"] = "low"
    elif exhaust_p >= cfg.P_HIGH:
        ex_sorted = sorted(
            (r for r in immune["checkpoints"] if r["gene"] in cfg.EXHAUSTION_CHECKPOINTS),
            key=lambda r: r["percentile"],
            reverse=True,
        )[:3]
        immune["exhaustion_evidence"] = [
            {"gene": r["gene"], "percentile": r["percentile"]} for r in ex_sorted
        ]
        immune["exhaustion_direction"] = "high"
    else:
        immune["exhaustion_evidence"] = []
        immune["exhaustion_direction"] = "normal"
    immune["phenotype_basis"] = (
        f"IFN-γ P{int(round(ifng_p))}、CD8A P{int(round(cd8_p))}、"
        f"耗竭检查点聚合 P{int(round(exhaust_p))}、基质 P{int(round(stromal_p))} 综合判定。"
    )

    # ── ③ 增殖与核心通路 ──
    pathways = []
    for pname in ["G2M_E2F", "NRF2_OXIDATIVE", "SQUAMOUS_STEMNESS", "EMT_INVASION"]:
        pathways.append(pathway_map.get(pname))
    pathways = [p for p in pathways if p is not None]
    drivers = {}
    for label, genes in DRIVER_GROUPS.items():
        drivers[label] = gene_group(genes)
    proliferation = {"pathways": pathways, "drivers": drivers}

    # ── ④ 靶点 ──
    targets = gene_group(read_gene_list(cfg.GENE_SETS_DIR / "targets_adc.txt"))

    return {
        "sample": sample,
        "reference": f"TCGA-ESCC 队列内参照（n={len(samples)}，留一法）",
        "qc": qc_fact,
        "immune": immune,
        "proliferation": proliferation,
        "targets": targets,
        "boundary": (
            "仅反映 bulk RNA-seq 转录层面相对趋势，不能替代 NGS/FISH 或 IHC 病理确诊。"
        ),
    }
