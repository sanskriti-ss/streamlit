"""
Phase 1d: Merge BacDive outliers, genome metadata, and antiSMASH into ranked candidates.

Usage:
  python -m genome_investigation.species_prioritization \\
    --genome-enriched genome_investigation/results/Step2_5_genome_enriched.csv \\
    --species-data-dir species_data
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from genome_investigation.io_utils import load_yaml, normalize_bacid

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WEIGHTS = Path(__file__).resolve().parent / "config" / "prioritization_weights.yaml"
DEFAULT_OUT = Path(__file__).resolve().parent / "results" / "ranked_species_candidates.csv"
SPECIES_DATA = REPO_ROOT / "species_data"

OUTPUT_COLUMNS = [
    "species",
    "best_BacID",
    "category",
    "bacdive_breadth",
    "bacdive_z",
    "genome_accession",
    "genome_quality",
    "bgc_count_total",
    "bgc_types",
    "key_genome_features",
    "discordance_score",
    "novelty_score",
    "safety_flag",
    "priority_score",
    "rationale",
]


def load_weights(path: Path) -> dict:
    if not path.exists():
        return {
            "outlier_z_weight": 2.0,
            "genome_confidence_weight": 1.5,
            "bgc_count_weight": 0.4,
            "unusual_bgc_bonus": 1.0,
            "resistance_unexplained_penalty": 1.2,
            "low_confidence_genome_penalty": 1.0,
            "specialist_bgc_penalty": 0.0,
        }
    return load_yaml(path)


def _zscore(s: pd.Series) -> pd.Series:
    s = s.astype(float)
    sd = float(s.std(ddof=0)) or 1.0
    return (s - float(s.mean())) / sd


def compute_bacdive_breadth(species_data_dir: Path) -> pd.DataFrame:
    """Compute per-strain breadth and z-scores for util/prod/res."""
    import zipfile

    files = {
        "utilization": species_data_dir / "step3_met_util_exploded.csv.zip",
        "production": species_data_dir / "step3_met_prod_exploded.csv.zip",
        "resistance": species_data_dir / "step3_met_res_exploded.csv.zip",
    }
    meta = ["BacID", "species", "genus", "order", "type_strain"]
    rows = []
    for activity, path in files.items():
        if not path.exists():
            continue
        with zipfile.ZipFile(path, "r") as zf:
            csvs = [n for n in zf.namelist() if n.endswith(".csv")]
            with zf.open(csvs[0]) as f:
                df = pd.read_csv(f, low_memory=False)
        mets = [c for c in df.columns if c not in meta and c != "species_with_id"]
        x = df[mets].apply(pd.to_numeric, errors="coerce")
        if activity in ("resistance",):
            breadth = (x == 1).sum(axis=1)
        else:
            breadth = (x.replace(-1, 0).fillna(0) > 0).sum(axis=1)
        part = df[meta].copy()
        part["activity"] = activity
        part["breadth"] = breadth.astype(int)
        part["z"] = _zscore(breadth)
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def assign_categories(stats: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    """Tag strains into prioritization categories from z-scores."""
    tagged = []
    for activity, cat_hi, cat_lo in [
        ("utilization", "metabolic_generalist_util", "metabolic_specialist_util"),
        ("production", "metabolic_generalist_prod", "metabolic_specialist_prod"),
        ("resistance", "resistance_outlier", None),
    ]:
        sub = stats[stats["activity"] == activity].copy()
        if sub.empty:
            continue
        hi = sub.nlargest(top_n, "z")
        hi["category"] = cat_hi
        tagged.append(hi)
        if cat_lo:
            lo = sub.nsmallest(top_n, "z")
            lo["category"] = cat_lo
            tagged.append(lo)
    return pd.concat(tagged, ignore_index=True)


def genome_quality_label(row: dict) -> str:
    conf = float(row.get("match_confidence") or 0)
    level = str(row.get("assembly_level") or "")
    if conf >= 0.9 and level.lower() in ("complete", "chromosome"):
        return "high"
    if conf >= 0.7:
        return "medium"
    if conf > 0:
        return "low"
    return "missing"


def build_rationale(row: dict) -> str:
    parts = []
    if row.get("category"):
        parts.append(f"Category: {row['category']}.")
    if row.get("bacdive_breadth") is not None:
        parts.append(f"BacDive breadth={row['bacdive_breadth']} (z={row.get('bacdive_z', 'NA'):.2f}).")
    if row.get("genome_accession"):
        parts.append(
            f"Genome {row['genome_accession']} ({row.get('genome_quality', 'unknown')} confidence, "
            f"match={row.get('match_confidence', 'NA')})."
        )
    else:
        parts.append("No high-confidence genome accession.")
    bgc = int(row.get("bgc_count_total") or 0)
    if bgc:
        parts.append(f"antiSMASH reports {bgc} BGC(s): {row.get('bgc_types', '')}.")
    if row.get("safety_flag"):
        parts.append(f"Flag: {row['safety_flag']}.")
    if row.get("discordance_score", 0) > 0.5:
        parts.append("Broad resistance phenotype with limited genomic explanation.")
    return " ".join(parts)


def score_candidate(row: dict, weights: dict) -> tuple[float, float, float, str]:
    z = float(row.get("bacdive_z") or 0)
    conf = float(row.get("match_confidence") or 0)
    bgc = int(row.get("bgc_count_total") or 0)
    category = str(row.get("category") or "")

    priority = abs(z) * weights.get("outlier_z_weight", 2.0)
    priority += conf * weights.get("genome_confidence_weight", 1.5)

    discordance = 0.0
    novelty = 0.0
    safety = ""

    if "production" in category and bgc > 0:
        priority += bgc * weights.get("bgc_count_weight", 0.4)
        if bgc >= 5 or (row.get("bgc_types") and len(str(row["bgc_types"]).split(";")) >= 3):
            priority += weights.get("unusual_bgc_bonus", 1.0)
            novelty += 0.5

    if "resistance" in category and z > 1.5:
        if bgc < 2 and conf < 0.8:
            discordance = weights.get("resistance_unexplained_penalty", 1.2)
            priority += discordance
            safety = "resistance_outlier_unexplained"

    if conf < 0.6:
        priority -= weights.get("low_confidence_genome_penalty", 1.0)

    if "specialist" in category:
        penalty = weights.get("specialist_bgc_penalty", 0.0)
        if penalty and bgc == 0:
            priority -= penalty  # default 0 — do not penalize specialists

    if "cereus" in str(row.get("species", "")).lower():
        safety = safety or "biosafety_review_recommended"

    return round(priority, 3), round(discordance, 3), round(novelty, 3), safety


def prioritize(
    tagged: pd.DataFrame,
    genome_enriched: pd.DataFrame,
    antismash: Optional[pd.DataFrame],
    weights: dict,
) -> pd.DataFrame:
    genome_enriched = genome_enriched.copy()
    genome_enriched["_bacid_key"] = genome_enriched["BacID"].map(normalize_bacid)
    genome_by_bacid = genome_enriched.set_index("_bacid_key", drop=False)
    if antismash is not None and not antismash.empty:
        asm_index = antismash.set_index(antismash["BacID"].astype(str), drop=False)
    else:
        asm_index = pd.DataFrame()

    out_rows: List[dict] = []
    for _, t in tagged.iterrows():
        bacid = normalize_bacid(t["BacID"])
        g = genome_by_bacid.loc[bacid].to_dict() if bacid in genome_by_bacid.index else {}
        a = asm_index.loc[bacid].to_dict() if bacid in asm_index.index else {}

        row = {
            "species": t.get("species"),
            "best_BacID": bacid,
            "category": t.get("category"),
            "bacdive_breadth": int(t.get("breadth", 0)),
            "bacdive_z": float(t.get("z", 0)),
            "genome_accession": g.get("genome_accession", ""),
            "match_confidence": g.get("match_confidence", 0),
            "assembly_level": g.get("assembly_level", ""),
            "bgc_count_total": int(a.get("bgc_count_total") or 0),
            "bgc_types": a.get("bgc_types", ""),
        }
        row["genome_quality"] = genome_quality_label(row)
        row["key_genome_features"] = "; ".join(
            x
            for x in [
                f"assembly={row.get('assembly_level')}" if row.get("assembly_level") else "",
                f"gc={g.get('gc_percent')}" if g.get("gc_percent") not in (None, "") else "",
                f"size_bp={g.get('genome_size_bp')}" if g.get("genome_size_bp") not in (None, "") else "",
            ]
            if x
        )
        prio, disc, nov, safety = score_candidate(row, weights)
        row["priority_score"] = prio
        row["discordance_score"] = disc
        row["novelty_score"] = nov
        row["safety_flag"] = safety
        row["rationale"] = build_rationale(row)
        out_rows.append({c: row.get(c, "") for c in OUTPUT_COLUMNS})

    df = pd.DataFrame(out_rows, columns=OUTPUT_COLUMNS)
    return df.sort_values("priority_score", ascending=False)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Phase 1d: species prioritization")
    parser.add_argument("--genome-enriched", required=True)
    parser.add_argument("--antismash-summary", default=str(Path(__file__).resolve().parent / "results" / "antismash_summary.csv"))
    parser.add_argument("--species-data-dir", default=str(SPECIES_DATA))
    parser.add_argument("--weights", default=str(DEFAULT_WEIGHTS))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--top-outliers-per-category", type=int, default=50)
    args = parser.parse_args(argv)

    genome = pd.read_csv(args.genome_enriched)
    asm_path = Path(args.antismash_summary)
    antismash = pd.read_csv(asm_path) if asm_path.exists() else None
    weights = load_weights(Path(args.weights))

    print("[info] computing BacDive breadth z-scores …")
    stats = compute_bacdive_breadth(Path(args.species_data_dir))
    tagged = assign_categories(stats, top_n=args.top_outliers_per_category)

    ranked = prioritize(tagged, genome, antismash, weights)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(out, index=False)
    print(f"[done] wrote {out} ({len(ranked)} candidates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
