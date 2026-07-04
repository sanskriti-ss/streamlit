"""
Load and merge cross-kingdom phenotype layers for the symbiosis network app.
"""

from __future__ import annotations

import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

LAYER_LABELS = {
    "bacteria_experimental": "Experimental bacteria",
    "bacteria_predicted": "Predicted bacteria",
    "fungi_experimental": "Experimental fungi",
    "fungi_predicted": "Predicted fungi",
}

UNIFIED_COLUMNS = [
    "entity_key",
    "kingdom",
    "species",
    "genus",
    "activity",
    "metabolite",
    "confidence",
    "observed",
    "source_layer",
]

METADATA_COLS = {
    "BacID",
    "species",
    "genus",
    "order",
    "type_strain",
    "is_strain",
    "species_with_id",
    "strain",
    "entity_id",
    "kingdom",
    "entity_key",
}


@dataclass
class LayerStatus:
    name: str
    label: str
    available: bool
    path: Optional[Path]
    row_count: int
    message: str = ""


def load_paths_config(config_path: Optional[Path] = None) -> dict:
    path = config_path or (REPO_ROOT / "symbiosis_data_paths.yaml")
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_metabolite_aliases(path: Optional[Path]) -> Dict[str, str]:
    """Map alternate metabolite names to canonical BacDive-style names."""
    if not path or not Path(path).exists():
        return {}
    df = pd.read_csv(path)
    if df.empty or "alias" not in df.columns or "canonical" not in df.columns:
        return {}
    return {
        str(a).strip().lower(): str(c).strip()
        for a, c in zip(df["alias"], df["canonical"])
        if pd.notna(a) and pd.notna(c)
    }


def normalize_metabolite(name: str, aliases: Dict[str, str]) -> str:
    key = str(name).strip()
    return aliases.get(key.lower(), key)


DEGRADATION_SUFFIX = "+degradation"


def load_metabolite_bridges(path: Optional[Path]) -> List[Tuple[str, str, float]]:
    """Load polymer→monomer degradation bridges.

    Each row means: a species that *utilizes* (extracellularly degrades) ``polymer``
    releases ``monomer`` into the shared environment — i.e. effectively *produces* it
    for cross-feeding partners. ``release_confidence`` scales the utilization
    confidence when deriving the inferred production row.
    """
    if not path or not Path(path).exists():
        return []
    df = pd.read_csv(path)
    if df.empty or "polymer" not in df.columns or "monomer" not in df.columns:
        return []
    out: List[Tuple[str, str, float]] = []
    for _, r in df.iterrows():
        polymer = str(r.get("polymer", "")).strip()
        monomer = str(r.get("monomer", "")).strip()
        if not polymer or not monomer:
            continue
        try:
            factor = float(r.get("release_confidence", 0.7))
        except (TypeError, ValueError):
            factor = 0.7
        factor = min(max(factor, 0.0), 1.0)
        out.append((polymer, monomer, factor))
    return out


def apply_degradation_bridges(
    df: pd.DataFrame,
    bridges: List[Tuple[str, str, float]],
    aliases: Dict[str, str],
) -> pd.DataFrame:
    """Derive production rows for monomers released by polymer degradation.

    Bridges the production/utilization ontology gap: BGC-based production and
    carbon-source utilization otherwise share no vocabulary, so synergy is always
    zero. Extracellular breakdown of a polymer (a utilization capability) yields
    monomers partners can consume, so we emit inferred production rows in the same
    carbon-source vocabulary. Original rows are preserved.
    """
    if df.empty or not bridges:
        return df

    bridge_map: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for polymer, monomer, factor in bridges:
        bridge_map[polymer.strip().lower()].append((monomer, factor))

    util = df[df["activity"] == "utilization"]
    if util.empty:
        return df

    new_rows: List[dict] = []
    for _, r in util.iterrows():
        met_key = str(r.get("metabolite", "")).strip().lower()
        targets = bridge_map.get(met_key)
        if not targets:
            continue
        base_conf = float(r.get("confidence", 1.0))
        for monomer, factor in targets:
            new_rows.append(
                {
                    "entity_key": r.get("entity_key"),
                    "kingdom": r.get("kingdom"),
                    "species": r.get("species"),
                    "genus": r.get("genus", ""),
                    "activity": "production",
                    "metabolite": normalize_metabolite(monomer, aliases),
                    "confidence": round(base_conf * factor, 3),
                    "observed": False,
                    "source_layer": f"{r.get('source_layer', '')}{DEGRADATION_SUFFIX}",
                }
            )

    if not new_rows:
        return df

    derived = pd.DataFrame(new_rows, columns=UNIFIED_COLUMNS)
    derived = derived.sort_values("confidence", ascending=False).drop_duplicates(
        subset=["entity_key", "activity", "metabolite", "source_layer"], keep="first"
    )
    return pd.concat([df, derived], ignore_index=True)


def _read_csv_or_zip(path: Path, *, nrows: Optional[int] = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            csvs = [n for n in zf.namelist() if n.lower().endswith(".csv")]
            if not csvs:
                raise ValueError(f"No CSV in zip: {path}")
            with zf.open(csvs[0]) as f:
                return pd.read_csv(f, low_memory=False, nrows=nrows)
    return pd.read_csv(path, low_memory=False, nrows=nrows)


def _metabolite_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c not in METADATA_COLS]


def _entity_key_from_row(row: pd.Series, kingdom: str) -> str:
    if pd.notna(row.get("entity_key")) and str(row["entity_key"]).strip():
        return str(row["entity_key"]).strip()
    if kingdom == "bacteria" and pd.notna(row.get("BacID")):
        sp = str(row.get("species", "") or "")
        return f"{row['BacID']} | {sp}" if sp else str(row["BacID"])
    if pd.notna(row.get("entity_id")):
        sp = str(row.get("species", "") or "")
        return f"{row['entity_id']} | {sp}" if sp else str(row["entity_id"])
    return str(row.get("species", "unknown"))


def _wide_to_long_positive(
    df: pd.DataFrame,
    activity: str,
    kingdom: str,
    source_layer: str,
    aliases: Dict[str, str],
) -> pd.DataFrame:
    """Convert wide matrix to long format, keeping only positive phenotypes."""
    mets = _metabolite_columns(df)
    if not mets:
        return pd.DataFrame(columns=UNIFIED_COLUMNS)

    rows: List[dict] = []
    for _, row in df.iterrows():
        entity_key = _entity_key_from_row(row, kingdom)
        species = str(row.get("species", "") or entity_key)
        genus = str(row.get("genus", "") or "")
        for met in mets:
            val = row[met]
            if pd.isna(val):
                continue
            try:
                v = float(val)
            except (TypeError, ValueError):
                continue
            if activity in ("utilization", "production"):
                if v <= 0:
                    continue
                confidence = 1.0
                observed = True
            else:
                if v != 1:
                    continue
                confidence = 1.0
                observed = True
            rows.append(
                {
                    "entity_key": entity_key,
                    "kingdom": kingdom,
                    "species": species,
                    "genus": genus,
                    "activity": activity,
                    "metabolite": normalize_metabolite(met, aliases),
                    "confidence": confidence,
                    "observed": observed,
                    "source_layer": source_layer,
                }
            )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def aggregate_by_genus(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate strain-level wide tables to genus × metabolite (any-positive)."""
    mets = _metabolite_columns(df)
    if "genus" not in df.columns or not mets:
        return df
    counts = df.groupby("genus")["species"].count().rename("species_count")
    met_max = df.groupby("genus")[mets].max()
    out = met_max.join(counts).reset_index()
    return out


def _apply_entity_keys_genus(df: pd.DataFrame, kingdom: str) -> pd.DataFrame:
    out = df.copy()
    out["entity_key"] = out["genus"].astype(str) + f" | ({kingdom} genus)"
    out["species"] = out["genus"].astype(str)
    return out


def _limit_entities_wide(
    prod_df: pd.DataFrame,
    util_df: pd.DataFrame,
    *,
    max_entities: int,
    aggregation: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if aggregation == "genus":
        prod_df = _apply_entity_keys_genus(aggregate_by_genus(prod_df), "bacteria")
        util_df = _apply_entity_keys_genus(aggregate_by_genus(util_df), "bacteria")
        common = sorted(set(prod_df["entity_key"]) & set(util_df["entity_key"]))
        prod_df = prod_df[prod_df["entity_key"].isin(common)]
        util_df = util_df[util_df["entity_key"].isin(common)]
        return prod_df, util_df

    mets = sorted(set(_metabolite_columns(prod_df)) & set(_metabolite_columns(util_df)))
    keys_prod = prod_df.apply(lambda r: _entity_key_from_row(r, "bacteria"), axis=1)
    keys_util = util_df.apply(lambda r: _entity_key_from_row(r, "bacteria"), axis=1)
    prod_df = prod_df.copy()
    util_df = util_df.copy()
    prod_df["_key"] = keys_prod
    util_df["_key"] = keys_util
    prod_sub = prod_df.drop_duplicates("_key").set_index("_key")
    util_sub = util_df.drop_duplicates("_key").set_index("_key")
    common = sorted(set(prod_sub.index) & set(util_sub.index))
    if not common:
        return prod_df.head(0), util_df.head(0)

    prod_pos = (prod_sub.loc[common, mets].apply(pd.to_numeric, errors="coerce").fillna(0) > 0).sum(axis=1)
    util_pos = (util_sub.loc[common, mets].apply(pd.to_numeric, errors="coerce").fillna(0) > 0).sum(axis=1)
    breadth = prod_pos + util_pos
    keep = breadth.sort_values(ascending=False).head(max_entities).index.tolist()
    prod_df = prod_df[prod_df["_key"].isin(keep)].drop(columns=["_key"])
    util_df = util_df[util_df["_key"].isin(keep)].drop(columns=["_key"])
    return prod_df, util_df


def load_bacteria_experimental(
    config: dict,
    aliases: Dict[str, str],
    *,
    max_entities: int = 500,
    aggregation: str = "species",
) -> Tuple[pd.DataFrame, LayerStatus]:
    label = LAYER_LABELS["bacteria_experimental"]
    paths = (config.get("bacteria_experimental") or {})
    prod_path = REPO_ROOT / paths.get("production", "")
    util_path = REPO_ROOT / paths.get("utilization", "")
    if not prod_path.exists() or not util_path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "bacteria_experimental",
            label,
            False,
            None,
            0,
            "BacDive Step3 production/utilization files not found",
        )
    try:
        prod_df = _read_csv_or_zip(prod_path)
        util_df = _read_csv_or_zip(util_path)
        prod_df, util_df = _limit_entities_wide(
            prod_df, util_df, max_entities=max_entities, aggregation=aggregation
        )
        prod_long = _wide_to_long_positive(
            prod_df, "production", "bacteria", "bacteria_experimental", aliases
        )
        util_long = _wide_to_long_positive(
            util_df, "utilization", "bacteria", "bacteria_experimental", aliases
        )
        out = pd.concat([prod_long, util_long], ignore_index=True)
        return out, LayerStatus(
            "bacteria_experimental",
            label,
            True,
            prod_path,
            len(out),
            f"Loaded {len(out):,} positive phenotype rows ({aggregation} level)",
        )
    except Exception as exc:
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "bacteria_experimental",
            label,
            False,
            prod_path,
            0,
            str(exc),
        )


def load_long_phenotype_file(
    path: Path,
    *,
    kingdom: str,
    source_layer: str,
    aliases: Dict[str, str],
    default_activity: Optional[str] = None,
) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLUMNS)

    rows: List[dict] = []
    for _, r in df.iterrows():
        activity = str(r.get("activity", default_activity or "")).strip()
        if not activity:
            continue
        met = normalize_metabolite(str(r.get("metabolite", "")), aliases)
        if not met:
            continue
        conf = float(r.get("confidence", 1.0))
        if conf <= 0:
            continue
        observed = bool(r.get("observed", conf >= 1.0))
        kingdom_val = str(r.get("kingdom", kingdom))
        species = str(r.get("species", "") or "")
        genus = str(r.get("genus", "") or "")
        entity_key = str(r.get("entity_key", "") or "")
        if not entity_key:
            if kingdom_val == "bacteria" and pd.notna(r.get("BacID")):
                entity_key = f"{r['BacID']} | {species}" if species else str(r["BacID"])
            elif pd.notna(r.get("entity_id")):
                entity_key = f"{r['entity_id']} | {species}" if species else str(r["entity_id"])
            else:
                entity_key = f"{species} | {source_layer}" if species else source_layer
        rows.append(
            {
                "entity_key": entity_key,
                "kingdom": kingdom_val,
                "species": species or entity_key,
                "genus": genus,
                "activity": activity,
                "metabolite": met,
                "confidence": conf,
                "observed": observed,
                "source_layer": source_layer,
            }
        )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def load_bacteria_predicted(
    config: dict, aliases: Dict[str, str]
) -> Tuple[pd.DataFrame, LayerStatus]:
    label = LAYER_LABELS["bacteria_predicted"]
    rel = config.get("bacteria_predicted", "")
    path = REPO_ROOT / rel if rel else None
    if not path or not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "bacteria_predicted",
            label,
            False,
            path,
            0,
            "phenotype_confidence.csv not found",
        )
    try:
        raw = pd.read_csv(path)
        rows: List[dict] = []
        for _, r in raw.iterrows():
            conf = float(r.get("confidence", 0))
            if conf <= 0:
                continue
            activity = str(r.get("activity", ""))
            if activity not in ("production", "utilization"):
                continue
            sp = str(r.get("species", ""))
            rows.append(
                {
                    "entity_key": f"{sp} | predicted",
                    "kingdom": "bacteria",
                    "species": sp,
                    "genus": sp.split()[0] if sp else "",
                    "activity": activity,
                    "metabolite": normalize_metabolite(str(r["metabolite"]), aliases),
                    "confidence": conf,
                    "observed": bool(r.get("observed", False)),
                    "source_layer": "bacteria_predicted",
                }
            )
        out = pd.DataFrame(rows, columns=UNIFIED_COLUMNS)
        return out, LayerStatus(
            "bacteria_predicted",
            label,
            True,
            path,
            len(out),
            f"Loaded {len(out):,} predicted phenotype rows",
        )
    except Exception as exc:
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "bacteria_predicted", label, False, path, 0, str(exc)
        )


def load_fungi_experimental(
    config: dict, aliases: Dict[str, str]
) -> Tuple[pd.DataFrame, LayerStatus]:
    label = LAYER_LABELS["fungi_experimental"]
    rel = config.get("fungi_experimental", "")
    path = REPO_ROOT / rel if rel else None
    if not path or not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "fungi_experimental",
            label,
            False,
            path,
            0,
            "FUNG-GROWTH data not found — see fungi_investigation/README.md",
        )
    try:
        if path.suffix.lower() == ".csv" and "wide" in path.name:
            df = pd.read_csv(path)
            out = _wide_to_long_positive(
                df, "utilization", "fungi", "fungi_experimental", aliases
            )
        else:
            out = load_long_phenotype_file(
                path, kingdom="fungi", source_layer="fungi_experimental", aliases=aliases
            )
        return out, LayerStatus(
            "fungi_experimental",
            label,
            True,
            path,
            len(out),
            f"Loaded {len(out):,} experimental fungi phenotype rows",
        )
    except Exception as exc:
        return pd.DataFrame(columns=UNIFIED_COLUMNS), LayerStatus(
            "fungi_experimental", label, False, path, 0, str(exc)
        )


def load_fungi_predicted(
    config: dict, aliases: Dict[str, str]
) -> Tuple[pd.DataFrame, List[LayerStatus]]:
    statuses: List[LayerStatus] = []
    parts: List[pd.DataFrame] = []

    for key, activity in (
        ("fungi_predicted_production", "production"),
        ("fungi_predicted_utilization", "utilization"),
    ):
        rel = config.get(key, "")
        path = REPO_ROOT / rel if rel else None
        label = LAYER_LABELS["fungi_predicted"]
        if not path or not path.exists():
            statuses.append(
                LayerStatus(
                    "fungi_predicted",
                    label,
                    False,
                    path,
                    0,
                    f"{path.name if path else key} not found",
                )
            )
            continue
        try:
            chunk = load_long_phenotype_file(
                path,
                kingdom="fungi",
                source_layer="fungi_predicted",
                aliases=aliases,
                default_activity=activity,
            )
            chunk["activity"] = activity
            parts.append(chunk)
            statuses.append(
                LayerStatus(
                    "fungi_predicted",
                    label,
                    True,
                    path,
                    len(chunk),
                    f"Loaded {len(chunk):,} rows from {path.name}",
                )
            )
        except Exception as exc:
            statuses.append(
                LayerStatus("fungi_predicted", label, False, path, 0, str(exc))
            )

    if not parts:
        return pd.DataFrame(columns=UNIFIED_COLUMNS), statuses
    return pd.concat(parts, ignore_index=True), statuses


def load_layers(
    selected_layers: Set[str],
    *,
    config_path: Optional[Path] = None,
    max_entities: int = 500,
    aggregation: str = "species",
    confidence_threshold: float = 0.5,
) -> Tuple[pd.DataFrame, List[LayerStatus]]:
    config = load_paths_config(config_path)
    alias_path = config.get("metabolite_aliases")
    aliases = load_metabolite_aliases(REPO_ROOT / alias_path if alias_path else None)

    frames: List[pd.DataFrame] = []
    statuses: List[LayerStatus] = []

    if "bacteria_experimental" in selected_layers:
        df, st = load_bacteria_experimental(
            config, aliases, max_entities=max_entities, aggregation=aggregation
        )
        statuses.append(st)
        if st.available and not df.empty:
            frames.append(df)

    if "bacteria_predicted" in selected_layers:
        df, st = load_bacteria_predicted(config, aliases)
        statuses.append(st)
        if st.available and not df.empty:
            frames.append(df)

    if "fungi_experimental" in selected_layers:
        df, st = load_fungi_experimental(config, aliases)
        statuses.append(st)
        if st.available and not df.empty:
            frames.append(df)

    if "fungi_predicted" in selected_layers:
        df, sts = load_fungi_predicted(config, aliases)
        statuses.extend(sts)
        if not df.empty:
            frames.append(df)

    if not frames:
        return pd.DataFrame(columns=UNIFIED_COLUMNS), statuses

    merged = pd.concat(frames, ignore_index=True)

    bridge_rel = config.get("metabolite_bridges")
    bridges = load_metabolite_bridges(REPO_ROOT / bridge_rel if bridge_rel else None)
    merged = apply_degradation_bridges(merged, bridges, aliases)

    merged = merged[merged["confidence"] >= confidence_threshold].copy()
    merged = merged.drop_duplicates(
        subset=["entity_key", "activity", "metabolite", "source_layer"], keep="first"
    )
    return merged, statuses


def list_metabolites(df: pd.DataFrame, query: str = "") -> List[str]:
    mets = sorted(df["metabolite"].dropna().unique())
    if not query:
        return mets
    q = query.strip().lower()
    return [m for m in mets if q in m.lower()]


def list_species(df: pd.DataFrame, query: str = "") -> List[str]:
    items = sorted(
        {
            f"{row['entity_key']} ({row['kingdom']})"
            for _, row in df[["entity_key", "kingdom"]].drop_duplicates().iterrows()
        }
    )
    if not query:
        return items
    q = query.strip().lower()
    return [s for s in items if q in s.lower()]
