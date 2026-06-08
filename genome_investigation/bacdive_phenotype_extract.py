"""
Extract BacDive utilization / production / resistance metabolite lists per strain.

Uses API cache when available; can fetch live with network.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import pandas as pd

from genome_investigation.api_cache import DEFAULT_CACHE_DIR, fetch_json

BACDIVE_API_BASE = "https://api.bacdive.dsmz.de"

PHENOTYPE_COLUMNS = [
    "BacID",
    "species",
    "strain",
    "activity",
    "metabolite",
    "chebi_id",
    "result",
    "detail",
    "source_ref",
]


def _norm_list(block: Any) -> List[dict]:
    if isinstance(block, list):
        return [x for x in block if isinstance(x, dict)]
    if isinstance(block, dict):
        return [block]
    return []


def _chebi(row: dict) -> str:
    for key in ("Chebi-ID", "Chebi ID", "ChEBI", "ChEBI-ID", "chebi_id"):
        if row.get(key) not in (None, ""):
            return str(row[key]).strip()
    return ""


def _activity_rows(strain: dict, bacid: str, species: str, strain_name: str) -> List[dict]:
    phys = strain.get("Physiology and metabolism") or strain.get("physiology_and_metabolism") or {}
    if not isinstance(phys, dict):
        return []

    rows: List[dict] = []

    for item in _norm_list(phys.get("metabolite utilization")):
        activity = str(item.get("utilization activity") or item.get("activity") or "").strip()
        rows.append(
            {
                "BacID": bacid,
                "species": species,
                "strain": strain_name,
                "activity": "utilization",
                "metabolite": str(item.get("metabolite") or "").strip(),
                "chebi_id": _chebi(item),
                "result": activity or "unknown",
                "detail": str(item.get("kind of utilization tested") or ""),
                "source_ref": str(item.get("@ref") or ""),
            }
        )

    for item in _norm_list(phys.get("metabolite production")):
        prod = str(item.get("production") or item.get("activity") or "").strip().lower()
        rows.append(
            {
                "BacID": bacid,
                "species": species,
                "strain": strain_name,
                "activity": "production",
                "metabolite": str(item.get("metabolite") or "").strip(),
                "chebi_id": _chebi(item),
                "result": "yes" if prod in ("yes", "+", "positive") else prod or "unknown",
                "detail": "",
                "source_ref": str(item.get("@ref") or ""),
            }
        )

    for item in _norm_list(phys.get("antibiotic resistance")):
        res = "resistant" if str(item.get("is resistant") or "").lower() in ("yes", "+") else "unknown"
        rows.append(
            {
                "BacID": bacid,
                "species": species,
                "strain": strain_name,
                "activity": "resistance",
                "metabolite": str(item.get("metabolite") or "").strip(),
                "chebi_id": _chebi(item),
                "result": res,
                "detail": str(item.get("resistance conc.") or item.get("resistance conc") or ""),
                "source_ref": str(item.get("@ref") or ""),
            }
        )

    for item in _norm_list(phys.get("antibiotic sensitivity")):
        sens = str(item.get("is sensitive") or "").lower()
        rows.append(
            {
                "BacID": bacid,
                "species": species,
                "strain": strain_name,
                "activity": "sensitivity",
                "metabolite": str(item.get("metabolite") or "").strip(),
                "chebi_id": _chebi(item),
                "result": "sensitive" if sens in ("yes", "+") else sens or "unknown",
                "detail": str(item.get("sensitivity conc.") or item.get("sensitivity conc") or ""),
                "source_ref": str(item.get("@ref") or ""),
            }
        )

    return [r for r in rows if r["metabolite"]]


def fetch_strain_payload(bacid: str, *, cache_dir: Path, force_refresh: bool = False) -> Optional[dict]:
    bacid_clean = re.sub(r"\.0$", "", str(bacid).strip())
    if not bacid_clean:
        return None
    url = f"{BACDIVE_API_BASE}/v2/fetch/{quote(bacid_clean)}"
    body, _ = fetch_json(url, cache_dir=cache_dir, force_refresh=force_refresh)
    if not isinstance(body, dict) or body.get("error"):
        return None
    results = body.get("results")
    if not isinstance(results, dict):
        return None
    strain = results.get(str(bacid_clean)) or results.get(bacid_clean)
    return strain if isinstance(strain, dict) else None


def extract_phenotype_rows(
    bacid: str,
    species: str = "",
    strain: str = "",
    *,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    force_refresh: bool = False,
) -> List[dict]:
    strain_data = fetch_strain_payload(bacid, cache_dir=cache_dir, force_refresh=force_refresh)
    if not strain_data:
        return []

    name_tax = strain_data.get("Name and taxonomic classification") or {}
    sp = species or (name_tax.get("species") if isinstance(name_tax, dict) else "") or ""
    st = strain or (name_tax.get("strain designation") if isinstance(name_tax, dict) else "") or ""
    return _activity_rows(strain_data, str(bacid), str(sp), str(st))


def build_phenotype_table(
    strains: pd.DataFrame,
    *,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    force_refresh: bool = False,
) -> pd.DataFrame:
    """strains: columns BacID, species (optional strain)."""
    rows: List[dict] = []
    for _, r in strains.iterrows():
        rows.extend(
            extract_phenotype_rows(
                str(r["BacID"]),
                str(r.get("species", "")),
                str(r.get("strain", "") or ""),
                cache_dir=cache_dir,
                force_refresh=force_refresh,
            )
        )
    if not rows:
        return pd.DataFrame(columns=PHENOTYPE_COLUMNS)
    return pd.DataFrame(rows, columns=PHENOTYPE_COLUMNS)


def positive_metabolites(df: pd.DataFrame, activity: str) -> pd.DataFrame:
    """Filter to positive utilization / production / resistance calls."""
    sub = df[df["activity"] == activity].copy()
    if activity == "utilization":
        return sub[sub["result"].isin(["+", "positive", "yes"])]
    if activity == "production":
        return sub[sub["result"].isin(["yes", "+", "positive"])]
    if activity == "resistance":
        return sub[sub["result"] == "resistant"]
    return sub
