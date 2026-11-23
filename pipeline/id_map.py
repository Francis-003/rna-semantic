"""File ID（目录 UUID）→ Case ID（TCGA barcode）映射读取。"""
from __future__ import annotations

import csv

from . import config as cfg


def load_case_map() -> dict[str, str]:
    """读取样本表，返回 File ID -> Case ID。

    列按表头名取值，避免硬编码列序；文件不存在或列缺失时给出清晰错误。
    """
    if not cfg.SAMPLE_SHEET.exists():
        raise FileNotFoundError(f"未找到样本映射表：{cfg.SAMPLE_SHEET}")
    mapping: dict[str, str] = {}
    with cfg.SAMPLE_SHEET.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "File ID" not in reader.fieldnames \
                or "Case ID" not in reader.fieldnames:
            raise ValueError(f"映射表缺少 File ID / Case ID 列：{reader.fieldnames}")
        for row in reader:
            fid = (row.get("File ID") or "").strip()
            cid = (row.get("Case ID") or "").strip()
            if fid and cid:
                mapping[fid] = cid
    return mapping


def to_case_id(file_id: str, case_map: dict[str, str]) -> str:
    """把 file id 转成 case id；查不到则回退原值并告警，保证管线不中断。"""
    case_id = case_map.get(file_id)
    if case_id is None:
        print(f"[id_map] 警告：file id {file_id} 无对应 Case ID，回退使用原 ID")
        return file_id
    return case_id
