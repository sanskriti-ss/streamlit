"""Shared I/O helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import yaml


def normalize_bacid(value: Any) -> str:
    """Normalize BacID for joins across CSVs (158570 vs 158570.0)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    s = str(value).strip()
    if s.lower() in ("nan", "none", ""):
        return ""
    if s.endswith(".0"):
        s = s[:-2]
    return s


def species_slug(name: str) -> str:
    s = re.sub(r"[^\w\s-]", "", str(name).strip().lower())
    s = re.sub(r"[-\s]+", "_", s)
    return s[:80] or "unknown_species"


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_selected_species(path: Path) -> Dict[str, Any]:
    cfg = load_yaml(path)
    species: List[str] = []
    bacdive_ids: List[str] = []

    if isinstance(cfg.get("species"), list):
        species.extend(str(s) for s in cfg["species"])
    if isinstance(cfg.get("bacdive_ids"), list):
        bacdive_ids.extend(str(x) for x in cfg["bacdive_ids"])

    for block in cfg.get("categories") or []:
        if isinstance(block, dict) and isinstance(block.get("species"), list):
            species.extend(str(s) for s in block["species"])

    return {
        "species": list(dict.fromkeys(species)),
        "bacdive_ids": list(dict.fromkeys(bacdive_ids)),
        "enable_antismash": bool(cfg.get("enable_antismash", False)),
        "raw": cfg,
    }
