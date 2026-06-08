"""Match confidence scoring for genome metadata lookups."""

from __future__ import annotations

import re
from typing import Any, Dict, Optional, Tuple

RowLike = Dict[str, Any]


def _norm(s: Any) -> str:
    if s is None:
        return ""
    if isinstance(s, float) and s != s:  # NaN
        return ""
    t = str(s).strip().lower()
    if t in ("nan", "none", ""):
        return ""
    return re.sub(r"\s+", " ", t)


def score_bacdive_direct() -> Tuple[float, str]:
    return 1.0, "BacDive strain record lists genome assembly accession"


def score_ncbi_match(
    row: RowLike,
    candidate: Dict[str, Any],
    *,
    allow_species_only: bool,
) -> Tuple[float, str]:
    """
    Score NCBI assembly match against input row.
    candidate expects keys: organism_name, assembly_accession, strain_name (optional)
    """
    species = _norm(row.get("species"))
    genus = _norm(row.get("genus"))
    strain = _norm(row.get("strain"))
    type_strain = _norm(row.get("type_strain"))

    org = _norm(candidate.get("organism_name"))
    asm_name = _norm(candidate.get("assembly_name", ""))
    cand_strain = _norm(candidate.get("strain_name", ""))

    if not species and not genus:
        return 0.0, "missing species/genus in input row"

    species_match = species and (species in org or species in asm_name)
    genus_match = genus and genus in org

    if not species_match and not genus_match:
        return 0.0, f"organism '{candidate.get('organism_name')}' does not match input taxon"

    strain_tokens = [t for t in re.split(r"[\s_\-]+", strain) if len(t) > 2]
    strain_hit = bool(strain) and (
        strain in org or strain in asm_name or cand_strain == strain or any(t in asm_name for t in strain_tokens)
    )
    type_strain_yes = type_strain in ("yes", "y", "true", "1")

    if species_match and strain_hit:
        return 0.92, "NCBI assembly matches species and strain designation"
    if species_match and type_strain_yes and "type strain" in asm_name:
        return 0.88, "NCBI assembly likely type strain for species"
    if species_match and strain and not strain_hit:
        return 0.72, "NCBI species match; strain not confirmed in assembly metadata"
    if species_match and allow_species_only:
        return 0.58, "NCBI species-level match only (strain not verified)"
    if genus_match and allow_species_only:
        return 0.45, "NCBI genus-level match only"
    return 0.0, "species match too weak without --allow-species-only"


def pick_best_candidate(
    row: RowLike,
    candidates: list[Dict[str, Any]],
    *,
    allow_species_only: bool,
) -> Tuple[Optional[Dict[str, Any]], float, str]:
    best: Optional[Dict[str, Any]] = None
    best_score = -1.0
    best_note = ""
    for cand in candidates:
        score, note = score_ncbi_match(row, cand, allow_species_only=allow_species_only)
        if score > best_score:
            best_score = score
            best = cand
            best_note = note
    if best is None or best_score <= 0:
        return None, 0.0, "no acceptable NCBI candidate"
    return best, best_score, best_note
