"""规则模板：结构化事实 JSON → 医生可读 Markdown 摘要。

原则：每个数值/百分位都必须来自事实 JSON；禁止无锚点的自由表述。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from . import config as cfg


def _p(rec: dict) -> str:
    """格式化百分位锚点，如 'P96'。文本层只保留百分位。

    显示取整与分档判定同口径：精确值处于中等区间（20 < p < 80）时，
    显示值截断在 21–79，避免出现“中等分位（P80）”（实际 79.8）式自相矛盾。
    """
    p = float(rec["percentile"])
    v = int(round(p))
    if cfg.P_LOW < p < cfg.P_HIGH:
        v = min(max(v, int(cfg.P_LOW) + 1), int(cfg.P_HIGH) - 1)
    return f"P{v}"


def _gene(rec: dict) -> str:
    """如 `CCND1` P96（文本层不含 Z / ES）。"""
    return f"`{rec['gene']}` {_p(rec)}"


def _gene_alias(name: str) -> str:
    """临床常用别名：TACSTD2 → TACSTD2/TROP2。"""
    return "TACSTD2/TROP2" if name == "TACSTD2" else name


def render_body(facts: dict) -> str:
    """100% 纯客观正文（四段）：测序质控 / 免疫 / 核心通路 / 常见靶点。

    只陈述可追溯到矩阵的统计分位与质控数值；不含标题、患者 ID 行、
    诊断标签、因果推断、临床建议或免责声明。供聚合 JSON 的 description 使用。
    """
    q = facts["qc"]
    tier = cfg.percentile_tier

    body: list[str] = []

    # ── 测序质控指标 ──
    qc_text = f"测序质控指标：线粒体转录本占比 {q['mt_pct']}%，未比对率 {q['unmapped_pct']}%。"
    if q["low_quality"]:
        qc_text += "GDC 标注：Item Flagged Low Quality。"
    body.append(qc_text)

    # ── 免疫相关特征与基因（客观功能名词 + 分位）──
    im = facts["immune"]
    ifng = im.get("ifng_gep")
    stromal = im.get("stromal")
    cd8 = im.get("cd8a")
    cps = im.get("checkpoints", [])
    cd274 = next((r for r in cps if r["gene"] == "CD274"), None)

    immune_parts = []
    if ifng:
        immune_parts.append(
            f"IFN-γ 炎症特征评分位于队列{tier(ifng['percentile'])}（{_p(ifng)}）"
        )
    # CD274/PD-L1 与 CD8A 合并为一句"位于<tier>"（统一带反引号）
    anchor_genes = []
    if cd274:
        anchor_genes.append(f"`CD274/PD-L1`（{_p(cd274)}）")
    if cd8:
        anchor_genes.append(f"`CD8A`（{_p(cd8)}）")
    if anchor_genes:
        # 若两者分档相同合并陈述，否则分别标注
        tiers = {tier(r["percentile"]) for r in (cd274, cd8) if r}
        if len(tiers) == 1:
            immune_parts.append(f"{'、'.join(anchor_genes)}位于{tiers.pop()}")
        else:
            if cd274:
                immune_parts.append(f"`CD274/PD-L1`（{_p(cd274)}）位于{tier(cd274['percentile'])}")
            if cd8:
                immune_parts.append(f"`CD8A`（{_p(cd8)}）位于{tier(cd8['percentile'])}")

    # 免疫检查点基因：按分档分组客观列出（不依赖 exhaustion_* 推断字段）。
    # CD274 已由上方锚点句单独播报，从汇总列表排除以免重复。
    hi_cp = sorted((r for r in cps if r["gene"] != "CD274" and tier(r["percentile"]) == "高分位"),
                   key=lambda r: r["percentile"], reverse=True)
    lo_cp = sorted((r for r in cps if r["gene"] != "CD274" and tier(r["percentile"]) == "低分位"),
                   key=lambda r: r["percentile"])
    if hi_cp:
        if len(hi_cp) == 1:
            immune_parts.append(f"免疫检查点基因 `{hi_cp[0]['gene']}`（{_p(hi_cp[0])}）位于高分位")
        else:
            shown = "、".join(f"`{r['gene']}` {_p(r)}" for r in hi_cp[:3])
            immune_parts.append(f"多个免疫检查点基因位于高分位（{shown}）")
    if lo_cp:
        if len(lo_cp) == 1:
            immune_parts.append(f"免疫检查点基因 `{lo_cp[0]['gene']}`（{_p(lo_cp[0])}）位于低分位")
        else:
            shown = "、".join(f"`{r['gene']}` {_p(r)}" for r in lo_cp[:3])
            immune_parts.append(f"多个免疫检查点基因位于低分位（{shown}）")

    if stromal:
        immune_parts.append(f"基质评分位于队列{tier(stromal['percentile'])}（{_p(stromal)}）")
    body.append("免疫相关特征与基因：" + "；".join(immune_parts) + "。")

    # ── 核心通路及代表基因 ──
    pr = facts["proliferation"]
    low_items: list[str] = []
    high_items: list[str] = []
    mid_anchor: list[str] = []
    med_pathway_extremes: list[str] = []  # 通路评分中等但成员基因处于极端分位（MKI67 除外）

    DRIVER_TO_PATHWAY = {
        "细胞周期与增殖": "G2M_E2F",
        "NRF2 氧化应激": "NRF2_OXIDATIVE",
        "鳞状分化与干性": "SQUAMOUS_STEMNESS",
        "EMT 与侵袭": "EMT_INVASION",
    }
    path_map = {p["name"]: {"p": p, "ev": [], "rev": []} for p in pr["pathways"]}
    for grp_label, recs in pr["drivers"].items():
        node = path_map.get(DRIVER_TO_PATHWAY.get(grp_label))
        if node is None:
            continue
        path_tier = tier(node["p"]["percentile"])
        # 按方向收纳：同向极端基因进“代表基因”，反向极端基因单独标注，
        # 避免低分位通路却列出高分位“代表基因”的语义矛盾（Bug 修复）。
        lo = sorted((r for r in recs if tier(r["percentile"]) == "低分位"),
                    key=lambda r: r["percentile"])
        hi = sorted((r for r in recs if tier(r["percentile"]) == "高分位"),
                    key=lambda r: r["percentile"])
        if path_tier == "低分位":
            node["ev"].extend(f"`{r['gene']}` {_p(r)}" for r in lo)
            node["rev"].extend(f"`{r['gene']}` {_p(r)}" for r in hi)
        elif path_tier == "高分位":
            node["ev"].extend(f"`{r['gene']}` {_p(r)}" for r in hi)
            node["rev"].extend(f"`{r['gene']}` {_p(r)}" for r in lo)
        else:
            # 通路评分中等 → 条目不输出；极端成员基因补报，避免静默丢失
            for r in lo + hi:
                if r["gene"] != "MKI67":
                    med_pathway_extremes.append(f"`{r['gene']}`（{_p(r)}）")
        # 必报锚点中处于中等分位者单列
        for r in recs:
            if r["gene"] in cfg.MUST_REPORT_GENES and tier(r["percentile"]) == "中等分位":
                mid_anchor.append(f"`{r['gene']}`（{_p(r)}）")

    # MKI67 全覆盖：G2M 通路处于中等分位（条目不输出）时，MKI67 无论自身分档必报，
    # 避免其极端表达证据随“中等通路不输出”规则静默丢失。
    g2m = path_map.get("G2M_E2F")
    mki67 = next((r for r in pr["drivers"].get("细胞周期与增殖", [])
                  if r["gene"] == "MKI67"), None)
    mki_anchor = None
    if g2m is not None and mki67 is not None:
        if (tier(g2m["p"]["percentile"]) == "中等分位"
                and tier(mki67["percentile"]) != "中等分位"):
            mki_anchor = f"`MKI67`（{_p(mki67)}）"

    for node in path_map.values():
        p = node["p"]
        t = tier(p["percentile"])
        if t not in ("低分位", "高分位"):
            continue
        text = f"{p['label']}通路评分（{_p(p)}"
        if node["ev"]:
            text += f"；代表基因 {'、'.join(node['ev'])}"
        if node["rev"]:
            dir_word = "反向低表达" if t == "高分位" else "反向高表达"
            text += f"；成员基因{dir_word} {'、'.join(node['rev'])}"
        text += "）"
        if t == "低分位":
            low_items.append(text)
        elif t == "高分位":
            high_items.append(text)

    # 血管生成相关基因 VEGFA：仅功能名词 + 分位（无推断）
    vegfa = next((r for r in facts["targets"] if r["gene"] in cfg.ANGIOGENESIS_TARGETS), None)
    if vegfa and tier(vegfa["percentile"]) == "高分位":
        high_items.append(f"血管生成相关基因 `VEGFA`（{_p(vegfa)}）")

    core_lines = ["核心通路及代表基因："]
    if low_items:
        core_lines.append(f"  低分位：{'；'.join(low_items)}。")
    if high_items:
        core_lines.append(f"  高分位：{'；'.join(high_items)}。")
    if mid_anchor:
        core_lines.append(f"  中等分位锚点基因：{'、'.join(mid_anchor)}。")
    if mki_anchor:
        core_lines.append(f"  增殖锚点基因：{mki_anchor}。")
    if med_pathway_extremes:
        core_lines.append(f"  成员基因极端表达（通路评分中等）：{'、'.join(med_pathway_extremes)}。")
    if len(core_lines) == 1:
        core_lines.append("  各核心通路评分均位于队列中等分位。")
    body.append("\n".join(core_lines))

    # ── 常见靶点基因表达（仅分位，无建议/价值判断）──
    # VEGFA 作为血管生成相关基因已在核心通路段呈现，靶点段不重复列出
    tg = [r for r in facts["targets"] if r["gene"] not in cfg.ANGIOGENESIS_TARGETS]
    hi_tg = sorted((r for r in tg if tier(r["percentile"]) == "高分位"),
                   key=lambda r: r["percentile"], reverse=True)
    lo_tg = sorted((r for r in tg if tier(r["percentile"]) == "低分位"),
                   key=lambda r: r["percentile"])
    mid_tg = [r for r in tg if tier(r["percentile"]) == "中等分位"]
    target_parts = []
    if hi_tg:
        target_parts.append("、".join(f"`{_gene_alias(r['gene'])}`（{_p(r)}）" for r in hi_tg) + "位于高分位")
    if lo_tg:
        target_parts.append("、".join(f"`{_gene_alias(r['gene'])}`（{_p(r)}）" for r in lo_tg) + "位于低分位")
    if mid_tg:
        target_parts.append("、".join(f"`{_gene_alias(r['gene'])}`" for r in mid_tg) + "位于中等分位")
    body.append("常见靶点基因表达：" + "；".join(target_parts) + "。")

    return "\n".join(body)


def render_markdown(facts: dict, case_id: str | None = None) -> str:
    """完整报告 = 标题 + 患者 ID 行 + 纯客观正文。

    case_id 为 TCGA Case ID；未提供时回退到 facts['sample']（file id）。
    """
    header = f"## 【患者转录组表达特征客观描述｜{facts['reference']}】\n"
    shown_id = case_id or facts["sample"]
    header += f"患者 ID：{shown_id}\n"
    return header + render_body(facts)


def _safe_output_path(directory: Path, filename: str) -> Path:
    """在写入前规范化并校验路径：禁止 ../、绝对路径越界，必须位于 directory 内。"""
    import os
    base = os.path.realpath(str(directory))
    candidate = os.path.realpath(os.path.join(base, filename))
    # 禁止路径穿越片段
    if ".." in filename.split("/") or ".." in filename.split("\\"):
        raise ValueError(f"非法文件名（含 ..）: {filename!r}")
    if os.path.commonpath([base, candidate]) != base:
        raise ValueError(f"输出路径越界: {candidate} 不在 {base} 内")
    return Path(candidate)


def write_outputs(
    facts: dict, out_dir=None, with_markdown: bool = True,
    case_id: str | None = None,
) -> dict:
    out = out_dir or cfg.OUTPUT_DIR
    feat_dir = out / "features"
    report_dir = out / "reports"
    feat_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    # 样本名仅作数据内容，不直接进路径：文件名由纯数字哈希派生，杜绝路径穿越。
    # 哈希仍基于 file id（facts['sample']），保持审计链与既有路径不变。
    import hashlib
    digest = hashlib.sha256(str(facts["sample"]).encode("utf-8")).hexdigest()[:12]
    safe_name = "s_" + digest

    json_path = _safe_output_path(feat_dir, safe_name + ".json")
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(facts, fh, ensure_ascii=False, indent=2)

    body = render_body(facts)
    md_path = None
    if with_markdown:
        md_path = _safe_output_path(report_dir, safe_name + ".md")
        with md_path.open("w", encoding="utf-8") as fh:
            fh.write(render_markdown(facts, case_id=case_id))
    return {
        "json": str(json_path),
        "markdown": str(md_path) if md_path else None,
        "body": body,
    }
