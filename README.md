# TCGA-ESCC bulk RNA-seq 表达特征客观语义化管线

把 GDC STAR 定量的 89 个食管鳞癌样本，从"样本内 Top20 高表达基因"转换为**队列参照的相对偏离度 + 通路聚合 + 靶点表达分位**的 100% 纯客观描述文本，作为下游结合随访结局做诊断/预后推理的干净输入。

> **输出定位**：报告只做客观数据描述——每句话都可追溯到 TSV 矩阵的数学计算结果，不含诊断标签、因果推断、临床建议或免责声明。主观/推断字段（phenotype、exhaustion_direction、amp_hint、boundary 等）仍保留在特征 JSON 中供审计，但不进入客观文本。

## 核心范式

| 旧做法 | 本管线 |
|---|---|
| 样本内绝对表达排名（全是 MT-/RPL/RPS） | 过滤线粒体/核糖体/非编码，仅保留蛋白编码基因做特征 |
| TPM=180 无参照系 | 留一法队列百分位 Pxx；客观三档 `≤P20 低分位 / 21–79 中等分位 / ≥P80 高分位`（严格固定阈值，Z 仅存 JSON 供审计） |
| 单基因列表 | ssGSEA 通路聚合 + 代表基因括号佐证 + 精选呈现（极端分位点名、中间态折叠）|

## 数据边界（务必先读）

- **输入**：`RNA/` 下 89 个 `*.rna_seq.augmented_star_gene_counts.tsv`（GENCODE v36，每样本一文件，60660 基因 × 6 定量列）。
- **无 TCGA barcode、无临床/分组/配对信息** → 无法区分肿瘤(01A)/正常(11A)，**放弃"肿瘤 vs 正常 Log2FC"维度**（已与需求方确认）。
- 8 个样本带 GDC `Item Flagged Low Quality` 注释 → 在报告中标注"需谨慎解读"，不静默剔除。
- `gene_name` 有 111 个重复 → 统一以 `gene_id`（带版本）为主键；一名多 ID 时取该样本表达最高者作展示值。

## 运行

环境：项目 `.venv`（已用 `--system-site-packages` 继承系统 pandas/numpy/scipy）。所有文本文件 UTF-8。

```bash
# 全量 89 样本（规则版文本）
.venv/Scripts/python.exe -m pipeline.main --samples all

# 单样本调试（UUID 前缀即可）
.venv/Scripts/python.exe -m pipeline.main --samples 0518caa3

# 可选：LLM 润色（无 key 自动降级为规则版，管线不中断）
#   PowerShell:  $env:LLM_API_KEY="sk-..."; $env:LLM_BASE_URL="https://api.openai.com"
#   Git Bash:    export LLM_API_KEY=sk-... LLM_BASE_URL=https://api.openai.com
.venv/Scripts/python.exe -m pipeline.main --samples all --polish
```

LLM 润色安全约束：仅允许 HTTPS + 域名白名单（`api.openai.com` / `api.anthropic.com` / `openai.googleapis.com`），解析 IP 必须为公网、连接固定到已校验 IP 防 DNS rebinding、禁代理与重定向；输入仅限结构化事实 JSON，system prompt 强制"只润色语言、保留全部数值与 P/Z 锚点、不得新增事实或临床结论"。

## 输出

- `outputs/features/s_<hash>.json`：结构化事实（QC、免疫检查点/IFN-γ GEP、通路、驱动基因、靶点的 percentile/zscore/direction、扩增提示标志、边界声明）。`sample` 字段仍为 file id（目录 UUID），保持审计链不变。
- `outputs/reports/s_<hash>.md`：100% 纯客观描述文本，**逐点分行列举**（患者 ID / 测序质控指标 / 免疫相关特征与基因 / 核心通路及代表基因「低分位 + 高分位 + 中等分位锚点」/ 常见靶点基因表达），只陈述可追溯到矩阵的统计分位与质控数值。
  - 文件名 `s_<hash>` 是 file id 的 SHA-256 前 12 位，仅作安全文件名（避免路径穿越）；报告的"患者 ID"行已按映射表显示为 TCGA Case ID（如 `TCGA-LN-A49N`）。
- `outputs/patient_descriptions.json`：**患者描述聚合文件**，结构 `{"<case_id>": {"description": "<纯客观正文>"}, ...}`，共 89 个患者；description 为去掉标题与 ID 行后的四段正文，供下游结合随访结局直接消费。**仅在 `--samples all` 全量运行时写入**。
- `outputs/cohort_summary.csv`：89 样本关键指标一览（MT%、低质量标注、免疫表型、通路分位、高表达靶点数）。**仅在 `--samples all` 全量运行时写入**，单样本调试不会覆盖。
- `outputs/cache/baseline.npz`：留一法基线缓存（按样本列表与矩阵形状自动失效）；需强制重算加 `--force-baseline`。

### ID 映射

根目录 `tcga_metadata_sample_sheet_ESCC_RNA_STAR_counts.tsv` 记录 File ID（= `RNA/` 目录 UUID = 管线内部 sample id）→ Case ID（TCGA barcode）的 1:1 映射（89 条，无重复）。`pipeline/id_map.py` 读取该表：报告的"患者 ID"行与聚合 JSON 的 key 使用 Case ID，特征 JSON 的文件名与 `sample` 字段仍用 File ID。映射缺失时回退原 ID 并告警。

## 方法细节

- **定量**：`tpm_unstranded` → `log2(TPM+1)`；`unstranded` 仅算 MT% 与比对汇总。
- **基线（留一法）**：样本 i 的基因分位 `P = 100*(count(others<=x_i)+1)/89`，`Z=(x_i-mean_others)/std_others(ddof=1)`；避免自身参与造成分位收缩。
- **客观分档**：文本层用 `percentile_tier()` 三档——`≥P80 高分位 / ≤P20 低分位 / 其余 中等分位`（严格固定阈值，中性措辞）；基因/代表基因按同一阈值点名或折叠。`MUST_REPORT_GENES = {CD274, MKI67}` 仍单列；VEGFA 作为"血管生成相关基因"只在核心通路段呈现，靶点段不重复。`direction()`/五级 `percentile_band` 仍保留供 JSON 与汇总表使用。
- **ssGSEA**：自实现秩次算法（`pipeline/pathways.py`），命中基因按 `|log2(TPM+1)|` 加权正向累加、未命中基因负向惩罚，取累积曲线极值；通路含 G2M/E2F、NRF2、鳞状分化、EMT、IFN-γ GEP、基质。**ES/Z 仅存 JSON，文本层不输出**。
- **客观语义化红线**：文本只保留基因/通路公认功能名词（客观本体知识）+ `Pxx` 统计分位；**不输出**诊断标签（phenotype/耗竭分型）、因果推断（"提示…可能存在"、"显著高表达"）、临床建议（"可结合 IHC/FISH"、"未见优势靶点"）、免责声明（解读边界）。这些主观字段仅存 JSON 供审计/下游自行取用。
- **基因功能注释 vs 状态推断**：保留 `G2M/E2F 细胞周期通路`、`免疫检查点基因`、`血管生成相关基因 VEGFA` 等数据库固有名词；剔除任何"提示/可能存在/倾向于/建议验证"等针对当前患者的推断句。

## 第二阶段（仅留接口，未实现）

`pipeline/features.py` 中为 ESTIMATE（纯度/基质）、免疫反卷积、TIDE、亚型映射预留了接入点。跨院泛化时按方案避坑点 3：用样本内秩次算法并重新校准分位数基线，不把外部医院 TPM 硬套 TCGA 的 Z-score。

## 验证结果（本次运行）

- 矩阵 `60660×89`，`log2(TPM+1)` 与原始 TSV 手算一致。
- 留一法 percentile/Z 独立重算与输出一致。
- MT% 范围 1.84%–25.08%；8 个低质量样本标注齐全。
- 89 份纯客观报告、89 份 JSON；严格禁用词表（诊断标签/推断/建议/免责声明）命中 0；所有 `Pxx` 可溯源到事实 JSON；`cohort_summary.csv` 89 行。
- 无 API key 时 `--polish` 正确降级为规则版（status=skipped），管线无中断。
- 客观描述基准样本 `5693ee07`：IFN-γ P97、CD274 P92、CD8A P88、CTLA4 P100/HAVCR2 P98/PDCD1 P96、基质 P62、G2M P4、EMT P91、VEGFA P100、MKI67 P21、MET P92、TACSTD2/TROP2 P100、FGFR2 P16 均与预期一致，且无任何主观措辞。
