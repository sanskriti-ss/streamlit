"""
Per-metabolite phenotype confidence (resistance / production / utilization).

Combines BacDive-tested calls (pinned to 1.0 / 0.0) with genomic evidence for
*untested* metabolites:
  - resistance:   AMRFinder gene Class/Subclass + %identity*%coverage
  - production:   antiSMASH BGC product types (putative, capped)
  - utilization:  CAZyme / transporter annotations in protein FASTA (saturating)

These are heuristic, genomic-evidence confidences, NOT calibrated probabilities.
Only BacDive-observed phenotypes reach 1.0; untested production/utilization are
capped below 1.0 because genomic potential does not prove the phenotype.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from genome_investigation.genome_paths import find_protein_fasta
from genome_investigation.io_utils import species_slug

CONF_COLUMNS = ["species", "activity", "metabolite", "confidence", "evidence", "observed"]

PRODUCTION_PUTATIVE_CAP = 0.65
UTILIZATION_CAP = 0.9
RESISTANCE_CAP = 0.98

# antiSMASH product token -> displayed compound class
BGC_PRODUCT_LABELS: Dict[str, str] = {
    "nrps": "Nonribosomal peptide",
    "t1pks": "Polyketide (type I)",
    "t2pks": "Polyketide (type II)",
    "t3pks": "Polyketide (type III)",
    "pks": "Polyketide",
    "transatpks": "Polyketide (trans-AT)",
    "terpene": "Terpene",
    "terpene-precursor": "Terpene",
    "ripp-like": "RiPP",
    "lassopeptide": "Lasso peptide",
    "proteusin": "Proteusin (RiPP)",
    "lanthipeptide": "Lanthipeptide",
    "thiopeptide": "Thiopeptide",
    "sactipeptide": "Sactipeptide",
    "rre-containing": "RiPP (RRE)",
    "siderophore": "Siderophore",
    " nrp-metallophore": "Metallophore",
    "ectoine": "Ectoine",
    "betalactone": "Beta-lactone",
    "cyclic-lactone-autoinducer": "Acyl-homoserine-lactone",
    "redox-cofactor": "Redox cofactor",
    "butyrolactone": "Butyrolactone",
}

# utilization substrate -> annotation keywords (lowercase substrings)
UTILIZATION_ENZYMES: Dict[str, List[str]] = {
    "cellulose": ["cellulase", "endoglucanase", "cellobiohydrolase", "glycoside hydrolase"],
    "xylan": ["xylanase", "xylosidase", "xylan"],
    "starch": ["amylase", "pullulanase", "glucoamylase", "starch"],
    "chitin": ["chitinase", "chitin"],
    "pectin": ["pectinase", "pectate lyase", "polygalacturonase"],
    "lactose": ["beta-galactosidase", "lacz", "lactose"],
    "sucrose": ["sucrase", "invertase", "levanase", "sucrose"],
    "mannan": ["mannanase", "mannosidase", "mannan"],
    "arabinose": ["arabinofuranosidase", "arabinosidase", "arabinose"],
    "glucose": ["glucokinase", "glucose-6", "pts system glucose"],
    "fructose": ["fructokinase", "fructose", "pts system fructose"],
    "maltose": ["maltose", "maltodextrin", "malphosphorylase"],
}

# AMRFinder Class/Subclass token -> displayed antibiotic class
RESISTANCE_DRUGS: Dict[str, str] = {
    "beta-lactam": "Beta-lactams",
    "carbapenem": "Carbapenems",
    "cephalosporin": "Cephalosporins",
    "tetracycline": "Tetracycline",
    "macrolide": "Macrolides",
    "aminoglycoside": "Aminoglycosides",
    "phenicol": "Chloramphenicol",
    "glycopeptide": "Vancomycin",
    "vancomycin": "Vancomycin",
    "fosfomycin": "Fosfomycin",
    "streptogramin": "Streptogramins",
    "lincosamide": "Lincosamides",
    "fluoroquinolone": "Fluoroquinolones",
    "quinolone": "Fluoroquinolones",
    "sulfonamide": "Sulfonamides",
    "trimethoprim": "Trimethoprim",
    "rifamycin": "Rifampicin",
    "bleomycin": "Bleomycin",
    "nitroimidazole": "Nitroimidazoles",
}


def _amr_tsv_path(amr_root: Path, species: str, accession: str) -> Optional[Path]:
    p = amr_root / species_slug(species) / str(accession) / "amrfinder.tsv"
    return p if p.exists() else None


def resistance_confidence(amr_tsv: Optional[Path]) -> Dict[str, float]:
    """Antibiotic class label -> best confidence from matching AMR genes."""
    out: Dict[str, float] = {}
    if not amr_tsv or not amr_tsv.exists() or amr_tsv.stat().st_size == 0:
        return out
    df = pd.read_csv(amr_tsv, sep="\t", low_memory=False)
    if df.empty:
        return out

    ident = pd.to_numeric(df.get("% Identity to reference"), errors="coerce").fillna(80) / 100.0
    cov = pd.to_numeric(df.get("% Coverage of reference"), errors="coerce").fillna(80) / 100.0
    cls = df.get("Class", pd.Series([""] * len(df))).astype(str)
    sub = df.get("Subclass", pd.Series([""] * len(df))).astype(str)
    blob = (cls + " " + sub).str.lower()

    for i, text in blob.items():
        # 0.5 (gene present) .. ~1.0 (perfect identity+coverage)
        score = 0.5 + 0.5 * float(ident.get(i, 0.8)) * float(cov.get(i, 0.8))
        score = round(min(score, RESISTANCE_CAP), 3)
        for token, label in RESISTANCE_DRUGS.items():
            if token in text:
                out[label] = max(out.get(label, 0.0), score)
    return out


def production_confidence(bgc_types: str) -> Dict[str, float]:
    """Compound class label -> putative production confidence (capped)."""
    out: Dict[str, float] = {}
    for tok in str(bgc_types or "").split(";"):
        tok = tok.strip().lower()
        if not tok:
            continue
        label = BGC_PRODUCT_LABELS.get(tok)
        if label is None:
            for key, val in BGC_PRODUCT_LABELS.items():
                if key.strip() and key.strip() in tok:
                    label = val
                    break
        if label:
            out[label] = max(out.get(label, 0.0), PRODUCTION_PUTATIVE_CAP)
    return out


def utilization_confidence(faa: Optional[Path]) -> Dict[str, float]:
    """Substrate label -> saturating confidence from enzyme annotation counts."""
    out: Dict[str, float] = {}
    if not faa or not faa.exists():
        return out
    counts = {s: 0 for s in UTILIZATION_ENZYMES}
    with faa.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line.lower()
            for sub, kws in UTILIZATION_ENZYMES.items():
                if any(kw in header for kw in kws):
                    counts[sub] += 1
    for sub, n in counts.items():
        if n > 0:
            # n=1 -> 0.5, n=2 -> 0.75, n=3 -> 0.875, saturating at cap
            out[sub.capitalize()] = round(min(UTILIZATION_CAP, 1 - 0.5 ** n), 3)
    return out


def _observed_calls(pheno_species: pd.DataFrame) -> Dict[tuple, float]:
    """Map (activity, metabolite) -> 1.0/0.0 from BacDive observed results."""
    observed: Dict[tuple, float] = {}
    for _, p in pheno_species.iterrows():
        act = str(p["activity"])
        met = str(p["metabolite"]).strip()
        res = str(p["result"])
        if not met:
            continue
        if act == "utilization":
            observed[(act, met)] = 1.0 if res in ("+", "positive", "yes") else 0.0
        elif act == "production":
            observed[(act, met)] = 1.0 if res in ("yes", "+", "positive") else 0.0
        elif act == "resistance":
            observed[(act, met)] = 1.0 if res == "resistant" else 0.0
        elif act == "sensitivity":
            # sensitive => observed non-resistant
            key = ("resistance", met)
            val = 0.0 if res == "sensitive" else observed.get(key, 0.0)
            observed[key] = val
    return observed


def build_confidence_table(
    focal: pd.DataFrame,
    pheno: pd.DataFrame,
    *,
    amr_root: Path,
    genome_root: Path,
    antismash_summary: pd.DataFrame,
) -> pd.DataFrame:
    """One row per (species, activity, metabolite) with confidence in [0, 1]."""
    rows: List[dict] = []

    for _, r in focal.iterrows():
        sp = str(r["species"])
        acc = str(r.get("genome_accession", "") or "")
        gdir = genome_root / species_slug(sp) / acc

        bgc_types = ""
        if not antismash_summary.empty and "species" in antismash_summary.columns:
            match = antismash_summary[antismash_summary["species"].astype(str) == sp]
            if not match.empty:
                bgc_types = str(match.iloc[0].get("bgc_types", "") or "")

        preds = {
            "resistance": resistance_confidence(_amr_tsv_path(amr_root, sp, acc)),
            "production": production_confidence(bgc_types),
            "utilization": utilization_confidence(find_protein_fasta(gdir)),
        }

        pheno_sp = pheno[pheno["species"].astype(str) == sp] if not pheno.empty else pd.DataFrame()
        observed = _observed_calls(pheno_sp) if not pheno_sp.empty else {}

        for activity, pred_map in preds.items():
            mets = set(pred_map) | {m for (a, m) in observed if a == activity}
            for met in sorted(mets):
                if (activity, met) in observed:
                    rows.append(
                        {
                            "species": sp,
                            "activity": activity,
                            "metabolite": met,
                            "confidence": observed[(activity, met)],
                            "evidence": "BacDive (observed)",
                            "observed": True,
                        }
                    )
                else:
                    rows.append(
                        {
                            "species": sp,
                            "activity": activity,
                            "metabolite": met,
                            "confidence": pred_map[met],
                            "evidence": f"genomic ({activity})",
                            "observed": False,
                        }
                    )

    if not rows:
        return pd.DataFrame(columns=CONF_COLUMNS)
    return pd.DataFrame(rows, columns=CONF_COLUMNS)
