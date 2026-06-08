"""
Lightweight gene-name search linking BacDive metabolites to downloaded annotations.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import pandas as pd

from genome_investigation.genome_paths import find_gff3, find_protein_fasta

HIT_COLUMNS = [
    "BacID",
    "species",
    "genome_accession",
    "activity",
    "metabolite",
    "search_term",
    "hit_source",
    "gene_id",
    "product",
    "match_score",
]


# Metabolite keyword → extra gene/protein search terms
METABOLITE_ALIASES: Dict[str, List[str]] = {
    "penicillin": ["penicillin", "bla", "beta-lactam"],
    "ampicillin": ["ampicillin", "bla", "beta-lactam"],
    "tetracycline": ["tetracycline", "tet"],
    "erythromycin": ["erythromycin", "erm"],
    "streptomycin": ["streptomycin", "str", "aminoglycoside"],
    "chloramphenicol": ["chloramphenicol", "cat", "cml"],
    "vancomycin": ["vancomycin", "van"],
    "gentamicin": ["gentamicin", "aac", "aph"],
    "cellulose": ["cellulase", "cel", "cellulose", "glycoside hydrolase"],
    "xylan": ["xylanase", "xylan", "xyl"],
    "glucose": ["glucose", "glc", "pts", "transporter"],
    "lactate": ["lactate", "ldh", "lactate dehydrogenase"],
    "acetate": ["acetate", "ack", "acetate kinase"],
}


def metabolite_search_terms(metabolite: str) -> List[str]:
    name = str(metabolite).strip().lower()
    if not name:
        return []
    terms: Set[str] = set()
    for part in re.split(r"[\s/(),]+", name):
        if len(part) >= 4:
            terms.add(part)
    for key, aliases in METABOLITE_ALIASES.items():
        if key in name:
            terms.update(aliases)
    return sorted(terms)


def _score_match(text: str, term: str) -> float:
    t = text.lower()
    if term in t:
        return 1.0 if re.search(rf"\b{re.escape(term)}\b", t) else 0.7
    return 0.0


def search_protein_fasta(
    faa_path: Path,
    terms: Iterable[str],
    *,
    max_hits_per_term: int = 5,
) -> List[dict]:
    hits: List[dict] = []
    term_list = list(terms)
    if not term_list or not faa_path.exists():
        return hits

    header = ""
    seq_buf: List[str] = []
    with faa_path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if header:
                    hits.extend(_match_fasta_record(header, "".join(seq_buf), term_list, max_hits_per_term))
                header = line[1:].strip()
                seq_buf = []
            else:
                seq_buf.append(line.strip())
        if header:
            hits.extend(_match_fasta_record(header, "".join(seq_buf), term_list, max_hits_per_term))
    return hits


def _match_fasta_record(header: str, _seq: str, terms: List[str], cap: int) -> List[dict]:
    gene_id = header.split()[0] if header else ""
    product = header
    out: List[dict] = []
    for term in terms:
        sc = _score_match(product, term)
        if sc > 0:
            out.append(
                {
                    "hit_source": "protein_fasta",
                    "gene_id": gene_id,
                    "product": product[:200],
                    "search_term": term,
                    "match_score": sc,
                }
            )
        if len(out) >= cap:
            break
    return out[:cap]


def search_gff3(gff_path: Path, terms: Iterable[str], *, max_hits: int = 200) -> List[dict]:
    hits: List[dict] = []
    term_list = list(terms)
    if not term_list or not gff_path.exists():
        return hits

    with gff_path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in ("CDS", "gene", "mRNA"):
                continue
            attr = parts[8]
            gene_id = _gff_attr(attr, "ID") or _gff_attr(attr, "locus_tag") or parts[0]
            product = _gff_attr(attr, "product") or _gff_attr(attr, "Name") or attr
            blob = f"{gene_id} {product} {attr}".lower()
            for term in term_list:
                sc = _score_match(blob, term)
                if sc > 0:
                    hits.append(
                        {
                            "hit_source": "gff3",
                            "gene_id": gene_id,
                            "product": product[:200],
                            "search_term": term,
                            "match_score": sc,
                        }
                    )
                    break
            if len(hits) >= max_hits:
                break
    return hits


def _gff_attr(attr: str, key: str) -> str:
    m = re.search(rf"{key}=([^;]+)", attr)
    return m.group(1).strip() if m else ""


def run_targeted_search(
    job: dict,
    phenotype_df: pd.DataFrame,
    amr_summary_row: Optional[dict] = None,
) -> pd.DataFrame:
    """Search proteins/GFF for BacDive metabolites (+ AMR gene symbols)."""
    gdir = Path(job["genome_dir"])
    sp = str(job.get("species", ""))
    bacid = str(job.get("BacID", ""))
    acc = str(job.get("genome_accession", ""))

    sub = phenotype_df[phenotype_df["BacID"].astype(str) == bacid]
    if sub.empty and sp:
        sub = phenotype_df[phenotype_df["species"].astype(str) == sp]

    faa = find_protein_fasta(gdir)
    gff = find_gff3(gdir)

    rows: List[dict] = []
    for _, phen in sub.iterrows():
        activity = str(phen["activity"])
        if activity == "utilization" and phen["result"] not in ("+", "positive", "yes"):
            continue
        if activity == "production" and phen["result"] not in ("yes", "+", "positive"):
            continue
        if activity == "resistance" and phen["result"] != "resistant":
            continue
        if activity == "sensitivity":
            continue

        metabolite = str(phen["metabolite"])
        terms = metabolite_search_terms(metabolite)
        raw_hits: List[dict] = []
        if faa:
            raw_hits.extend(search_protein_fasta(faa, terms))
        if gff:
            raw_hits.extend(search_gff3(gff, terms))

        for h in raw_hits:
            rows.append(
                {
                    "BacID": bacid,
                    "species": sp,
                    "genome_accession": acc,
                    "activity": activity,
                    "metabolite": metabolite,
                    **h,
                }
            )

    amr_count = pd.to_numeric((amr_summary_row or {}).get("amr_gene_count"), errors="coerce")
    amr_count = 0 if pd.isna(amr_count) else int(amr_count)
    if amr_summary_row and amr_count > 0:
        for gene in str(amr_summary_row.get("amr_genes", "")).split(";"):
            gene = gene.strip()
            if not gene:
                continue
            rows.append(
                {
                    "BacID": bacid,
                    "species": sp,
                    "genome_accession": acc,
                    "activity": "resistance",
                    "metabolite": "(AMRFinder)",
                    "search_term": gene.lower(),
                    "hit_source": "amrfinder",
                    "gene_id": gene,
                    "product": f"AMRFinder hit: {gene}",
                    "match_score": 1.0,
                }
            )

    if not rows:
        return pd.DataFrame(columns=HIT_COLUMNS)
    return pd.DataFrame(rows, columns=HIT_COLUMNS)
