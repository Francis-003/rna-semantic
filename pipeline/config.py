"""全局配置：路径、阈值、LLM 设置。"""
from __future__ import annotations

import os
from pathlib import Path

# 项目根目录（本文件位于 <root>/pipeline/config.py）
ROOT = Path(__file__).resolve().parent.parent
RNA_DIR = ROOT / "RNA"
GENE_SETS_DIR = ROOT / "gene_sets"
OUTPUT_DIR = ROOT / "outputs"
FEATURES_DIR = OUTPUT_DIR / "features"
REPORTS_DIR = OUTPUT_DIR / "reports"
CACHE_DIR = OUTPUT_DIR / "cache"
COHORT_SUMMARY = OUTPUT_DIR / "cohort_summary.csv"
SAMPLE_SHEET = ROOT / "tcga_metadata_sample_sheet_ESCC_RNA_STAR_counts.tsv"
PATIENT_DESCRIPTIONS = OUTPUT_DIR / "patient_descriptions.json"

# ── 五级分位阈值（文本层统一口径，纯百分位驱动；Z 仅保留在 JSON 供审计）──
P_EXTREME_LOW = 10.0   # ≤P10 极低
P_LOW = 20.0           # ≤P20 偏低
P_HIGH = 80.0          # ≥P80 偏高
P_EXTREME_HIGH = 90.0  # ≥P90 显著高

# 单基因点名阈值：只点名 ≥P85 或 ≤P15 的极端基因，中间态折叠
GENE_EXTREME_LOW = 15.0
GENE_EXTREME_HIGH = 85.0

# Z 阈值保留供后台审计，不再参与文本/方向判定
Z_HIGH = 1.5
Z_LOW = -1.5

# 靶点分组话术（临床语义）
MEMBRANE_AMP_TARGETS = {"ERBB2", "EGFR", "FGFR1", "FGFR2", "MET", "TACSTD2", "NECTIN4"}
ANGIOGENESIS_TARGETS = {"VEGFA"}
# 耗竭/共抑制检查点（checkpoints.txt 去掉 CD274 这类配体）
EXHAUSTION_CHECKPOINTS = {
    "PDCD1", "LAG3", "TIGIT", "CTLA4", "PDCD1LG2", "HAVCR2", "BTLA",
}
# 必报锚点：即使处于中间态也在文本中点名
MUST_REPORT_GENES = {"CD274", "MKI67"}


def percentile_band(p: float) -> str:
    """返回五级定性标签：极低/偏低/中等/偏高/显著高。"""
    if p <= P_EXTREME_LOW:
        return "极低"
    if p <= P_LOW:
        return "偏低"
    if p >= P_EXTREME_HIGH:
        return "显著高"
    if p >= P_HIGH:
        return "偏高"
    return "中等"


def percentile_tier(p: float) -> str:
    """纯客观三档分位标签（严格对应固定阈值），供客观语义化文本使用。

    ≥P80 → 高分位；≤P20 → 低分位；其余 → 中等分位。
    """
    if p >= P_HIGH:
        return "高分位"
    if p <= P_LOW:
        return "低分位"
    return "中等分位"

# 定量列（GDC STAR augmented 输出）
COL_TPM = "tpm_unstranded"
COL_FPKM = "fpkm_unstranded"
COL_FPKM_UQ = "fpkm_uq_unstranded"
COL_COUNT = "unstranded"

SUMMARY_PREFIX = "N_"  # N_unmapped 等汇总行前缀


def ensure_output_dirs() -> None:
    for d in (OUTPUT_DIR, FEATURES_DIR, REPORTS_DIR, CACHE_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_path(base: Path, *parts: str) -> Path:
    """安全拼接路径：禁止 ..、路径分隔符、盘符与空字节，且规范化后必须仍在 base 内。

    用于所有输出文件写入，防止样本名等片段造成路径穿越。
    """
    root = base.resolve()
    for part in parts:
        if not isinstance(part, str) or not part:
            raise ValueError(f"非法路径片段: {part!r}")
        if part in (".", "..") or any(
            ch in part for ch in ("/", "\\", ":", "\x00")
        ):
            raise ValueError(f"非法路径片段: {part!r}")
    target = base.joinpath(*parts).resolve()
    try:
        target.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"路径越界: {target} 不在 {root} 内") from exc
    return target


def llm_config() -> dict:
    """从环境变量读取 OpenAI 兼容接口配置；无 key 返回 enabled=False。"""
    api_key = os.environ.get("LLM_API_KEY", "").strip()
    return {
        "enabled": bool(api_key),
        "api_key": api_key,
        "base_url": os.environ.get(
            "LLM_BASE_URL", "https://api.openai.com/v1"
        ).strip().rstrip("/"),
        "model": os.environ.get("LLM_MODEL", "gpt-4o-mini").strip(),
        "timeout": float(os.environ.get("LLM_TIMEOUT", "60")),
    }
