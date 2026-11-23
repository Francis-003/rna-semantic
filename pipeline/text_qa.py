"""队列级文本质检。

用法：python -m pipeline.text_qa
"""
from __future__ import annotations

import json
import re
import sys

from . import config as cfg
from .render import _p, render_body


def _band_of_display(d):
    """显示整数 → 分档标签（口径同 config.percentile_tier：≥80 高、≤20 低）。"""
    if d >= 80:
        return "高分位"
    if d <= 20:
        return "低分位"
    return "中等分位"


def _case_feature_map():
    """经 reports/*.md 的患者 ID 行建立 case_id → 特征文件映射。"""
    mapping = {}
    for md in sorted(cfg.REPORTS_DIR.glob("s_*.md")):
        if ".rule." in md.name:
            continue
        m = re.search(r"患者 ID：(\S+)", md.read_text(encoding="utf-8"))
        if m:
            mapping[m.group(1)] = cfg.FEATURES_DIR / f"{md.stem}.json"
    return mapping


def _summary(stats, failures):
    print(json.dumps(stats, ensure_ascii=False, indent=2))
    if failures:
        print(f"QA 未通过：{len(failures)} 项失败")
        for f in failures[:30]:
            print("  -", f)
    else:
        print("QA 全部通过。")


def run():
    descs = json.loads(cfg.PATIENT_DESCRIPTIONS.read_text(encoding="utf-8"))
    mapping = _case_feature_map()
    failures = []
    stats = {
        "n_cases": len(descs),
        "n_rev_annotation": 0,
        "n_mki_anchor": 0,
        "n_med_extreme": 0,
        "n_single_cp": 0,
        "n_percentile_values_checked": 0,
    }
    for case, feat_path in mapping.items():
        if case not in descs:
            failures.append(f"{case}: 不在 patient_descriptions.json 中")
            continue
        text = descs[case]["description"]
        facts = json.load(open(feat_path, encoding="utf-8"))
        _check_case(text, facts, case, failures, stats)
    _summary(stats, failures)
    return not failures


def _check_genes(genes_seg, expect, case, failures):
    """校验基因段显示方向：expect 为 “低”/“高”。基因形如 `TP63` P3。"""
    pairs = re.findall(r"`([A-Z0-9/]+)` P(\d+)", genes_seg)
    for gname, gdisp in pairs:
        g = int(gdisp)
        if expect == "低" and g > 20:
            failures.append(f"{case}: 低方向基因 {gname} P{g} 超出低分位")
        if expect == "高" and g < 80:
            failures.append(f"{case}: 高方向基因 {gname} P{g} 低于高分位")


def _items(text, band):
    """提取核心通路段指定分档行（如“低分位：A；B。”）的条目列表。"""
    m = re.search(rf"{band}：([^。]*)。", text)
    seg = m.group(1) if m else ""
    if not seg:
        return []
    return seg.split("；")


def _check_case(text, facts, case, failures, stats):
    """对单例描述做全部文本层检查。"""
    # ── 1. MKI67 恰好一次 ──
    n = text.count("MKI67")
    if n != 1:
        failures.append(f"{case}: MKI67 出现 {n} 次（期望 1）")

    # ── 2. 通路条目方向一致 ──
    for band in ("低分位", "高分位"):
        for item in _items(text, band):
            pm = re.search(r"通路评分（P(\d+)；", item)
            if pm is None:
                pm = re.search(r"通路评分（P(\d+)）", item)
            if pm is None:
                continue
            d = int(pm.group(1))
            if _band_of_display(d) != band:
                failures.append(f"{case}: 通路显示 P{d} 与 {band} 行矛盾")
            seg = item.split("；成员基因反向")[0]
            if "代表基因" in seg:
                genes_seg = seg.split("代表基因 ")[1]
                _check_genes(genes_seg, band[0], case, failures)
            rm = re.search(r"；成员基因反向(高|低)表达 ([^）]*)）", item)
            if rm:
                direction = rm.group(1)
                _check_genes(rm.group(2), direction, case, failures)
                if band == "低分位" and direction == "低":
                    failures.append(f"{case}: 低分位行出现反向低表达标注")
                if band == "高分位" and direction == "高":
                    failures.append(f"{case}: 高分位行出现反向高表达标注")

    # ── 3. 检查点播报 ──
    for band in ("高", "低"):
        m = re.search(rf"多个免疫检查点基因位于{band}分位（([^）]*)）", text)
        if m:
            inner = m.group(1)
            if "、" not in inner:
                failures.append(f"{case}: “多个”检查点括号内仅 1 个基因：{inner}")
            if "CD274" in inner:
                failures.append(f"{case}: CD274 重复出现在检查点汇总括号内")
        sm = re.search(rf"免疫检查点基因 `([A-Z0-9/]+)`（P(\d+)）位于{band}分位", text)
        if sm:
            stats["n_single_cp"] += 1

    # ── 4. 格式 ──
    if "） " in text:
        failures.append(f"{case}: 中文右括号后半角空格")
    if "。 GDC" in text:
        failures.append(f"{case}: 句号后半角空格（GDC 行）")
    if "`CD274/PD-L1`（P" not in text:
        failures.append(f"{case}: CD274/PD-L1 锚点缺少反引号")
    if "`CD8A`（P" not in text:
        failures.append(f"{case}: CD8A 锚点缺少反引号")

    # ── 5. 显示与分档一致（遍历特征 JSON 全部百分位值）──
    for value in _iter_percentiles(facts):
        shown = int(_p({"percentile": value})[1:])
        stats["n_percentile_values_checked"] += 1
        if _band_of_display(shown) != cfg.percentile_tier(value):
            failures.append(f"{case}: 显示 P{shown} 跨档（精确值 {value}）")

    # ── 6. 产物同步：描述与当前 render 代码重渲染一致 ──
    fresh = render_body(facts)
    if fresh != text:
        failures.append(f"{case}: 描述与重渲染结果不一致（陈旧产物）")

    # ── 统计 ──
    if "成员基因反向" in text:
        stats["n_rev_annotation"] += 1
    if "增殖锚点基因" in text:
        stats["n_mki_anchor"] += 1
    if "成员基因极端表达" in text:
        stats["n_med_extreme"] += 1
    return True


def _iter_percentiles(node):
    """递归收集 facts 中所有 percentile 数值。"""
    out = []
    if isinstance(node, dict):
        if "percentile" in node:
            out.append(node["percentile"])
        for k in node:
            out.extend(_iter_percentiles(node[k]))
    elif isinstance(node, list):
        for sub in node:
            out.extend(_iter_percentiles(sub))
    return out


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 1)
