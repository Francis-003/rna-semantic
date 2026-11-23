"""CLI 入口：TCGA-ESCC bulk RNA-seq 临床语义化管线。

用法：
  python -m pipeline.main --samples all             # 全量，规则版文本
  python -m pipeline.main --samples all --polish    # 全量 + 可选 LLM 润色
  python -m pipeline.main --samples <uuid前缀>      # 单样本调试
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from . import config as cfg
from . import load, qc, baseline, features, render, llm_polish, pathways, id_map


def _resolve_samples(all_samples: list[str], query: str) -> list[str]:
    if query.lower() == "all":
        return all_samples
    q = query.lower()
    hits = [s for s in all_samples if s.lower().startswith(q)]
    if not hits:
        raise SystemExit(f"未找到匹配样本：{query}（共 {len(all_samples)} 个样本）")
    if len(hits) > 1:
        print(f"[main] 前缀匹配到 {len(hits)} 个样本，将全部处理：{hits}")
    return hits


def write_cohort_summary(rows: list[dict]) -> None:
    cfg.ensure_output_dirs()
    fields = [
        "sample", "mt_pct", "unmapped_pct", "low_quality",
        "ifng_gep_percentile", "stromal_percentile", "immune_phenotype",
        "nrf2_percentile", "g2m_percentile", "n_high_targets",
    ]
    # COHORT_SUMMARY 是由 config.py 中代码常量派生的固定路径（OUTPUT_DIR/字面量文件名），
    # 不含任何外部输入；用 Path.open() 写死在受控目录内。
    with cfg.COHORT_SUMMARY.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fields})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", default="all", help="all 或样本 UUID 前缀")
    parser.add_argument("--polish", action="store_true", help="尝试用 LLM 润色文本（无 key 自动降级）")
    parser.add_argument("--force-baseline", action="store_true", help="强制重算队列基线缓存")
    args = parser.parse_args()

    cfg.ensure_output_dirs()
    print("[1/6] 加载 89 个 STAR 定量 TSV ...")
    data = load.load_all()
    print(f"      矩阵形状：{data['log2_tpm'].shape[0]} 基因 × {len(data['samples'])} 样本")

    print("[2/6] 计算质控指标 ...")
    qc_results = qc.compute_qc(data)

    print("[3/6] 计算队列参照基线（留一法 Percentile / Z-score）...")
    bl = baseline.compute_baseline(data, force=args.force_baseline)

    print("      预计算 ssGSEA 通路打分（全队列一次）...")
    ss = pathways.ssGSEA_scores(data)

    print("      加载 File ID → Case ID 映射表 ...")
    case_map = id_map.load_case_map()
    print(f"      映射条目：{len(case_map)}")

    targets = _resolve_samples(data["samples"], args.samples)
    print(f"[4/6] 构建结构化特征并渲染文本（{len(targets)} 个样本）...")
    summary_rows = []
    patient_descriptions: dict[str, dict] = {}
    for idx, s in enumerate(targets, 1):
        case_id = id_map.to_case_id(s, case_map)
        facts = features.build_facts(data, bl, qc_results, s, ss=ss)
        paths = render.write_outputs(facts, case_id=case_id)
        patient_descriptions[case_id] = {"description": paths["body"]}
        status = "rule"
        if args.polish:
            md_path = Path(paths["markdown"])
            with md_path.open("r", encoding="utf-8") as fh:
                md = fh.read()
            polished, status = llm_polish.polish_markdown(md, facts)
            if status == "polished":
                # 保留规则版供审计（同一受控目录内，仅改后缀）
                rule_path = md_path.with_suffix(".rule.md")
                with rule_path.open("w", encoding="utf-8") as fh:
                    fh.write(md)
                with md_path.open("w", encoding="utf-8") as fh:
                    fh.write(polished)
        im = facts["immune"]
        pr = facts["proliferation"]
        pmap = {p["name"]: p["percentile"] for p in pr["pathways"]}
        summary_rows.append({
            "sample": s,
            "mt_pct": facts["qc"]["mt_pct"],
            "unmapped_pct": facts["qc"]["unmapped_pct"],
            "low_quality": facts["qc"]["low_quality"],
            "ifng_gep_percentile": im["ifng_gep"]["percentile"] if im["ifng_gep"] else "",
            "stromal_percentile": im["stromal"]["percentile"] if im["stromal"] else "",
            "immune_phenotype": im["phenotype"],
            "nrf2_percentile": pmap.get("NRF2_OXIDATIVE", ""),
            "g2m_percentile": pmap.get("G2M_E2F", ""),
            "n_high_targets": sum(1 for t in facts["targets"] if t["direction"] == "high"),
        })
        if idx <= 3 or idx % 20 == 0 or idx == len(targets):
            print(f"      [{idx}/{len(targets)}] {s[:8]}... 文本={status}")

    # 汇总表与聚合 JSON 仅全量运行时写入，避免单样本调试覆盖完整队列产物
    if len(targets) == len(data["samples"]):
        print("[5/6] 写入队列汇总表与患者描述聚合 JSON（全量）...")
        write_cohort_summary(summary_rows)
        with cfg.PATIENT_DESCRIPTIONS.open("w", encoding="utf-8") as fh:
            json.dump(patient_descriptions, fh, ensure_ascii=False, indent=2)
    else:
        print(f"[5/6] 跳过汇总表与聚合 JSON 写入（本次仅 {len(targets)} 个样本，保留已有全量产物）")

    print("[6/6] 完成。")
    print(f"      特征 JSON：{cfg.FEATURES_DIR}")
    print(f"      临床文本  ：{cfg.REPORTS_DIR}")
    if len(targets) == len(data["samples"]):
        print(f"      队列汇总  ：{cfg.COHORT_SUMMARY}")
        print(f"      患者描述  ：{cfg.PATIENT_DESCRIPTIONS}")


if __name__ == "__main__":
    main()
