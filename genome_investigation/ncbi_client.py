"""NCBI Datasets API client for genome metadata fallback."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from urllib.parse import quote

from genome_investigation.api_cache import fetch_json
from genome_investigation.match_scoring import pick_best_candidate

NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


def _report_to_candidate(report: dict) -> Dict[str, Any]:
    asm = report.get("assembly_info") or {}
    org = report.get("organism") or {}
    stats = asm.get("assembly_stats") or {}
    ann = report.get("annotation_info") or {}
    stats_ann = ann.get("stats") or {}

    return {
        "genome_accession": (
            report.get("accession")
            or report.get("current_accession")
            or asm.get("assembly_accession")
            or asm.get("accession")
        ),
        "assembly_level": asm.get("assembly_level") or asm.get("assembly_level_name"),
        "genome_size_bp": stats.get("total_sequence_length") or stats.get("genome_size"),
        "gc_percent": stats.get("gc_percent"),
        "gene_count": stats_ann.get("gene_count"),
        "cds_count": stats_ann.get("cds_count"),
        "ncbi_taxid": org.get("tax_id"),
        "organism_name": org.get("organism_name"),
        "assembly_name": asm.get("assembly_name"),
        "strain_name": org.get("infraspecific_names", {}).get("strain")
        if isinstance(org.get("infraspecific_names"), dict)
        else None,
        "source_database": "ncbi_datasets",
    }


def search_genomes_by_taxon(
    taxon: str,
    *,
    cache_dir,
    force_refresh: bool = False,
    limit: int = 10,
) -> List[Dict[str, Any]]:
    taxon_q = quote(taxon.strip())
    url = f"{NCBI_DATASETS_BASE}/genome/taxon/{taxon_q}/dataset_report?limit={limit}"
    body, _ = fetch_json(url, cache_dir=cache_dir, force_refresh=force_refresh)
    if not isinstance(body, dict):
        return []
    reports = body.get("reports") or []
    return [_report_to_candidate(r) for r in reports if isinstance(r, dict)]


def enrich_from_ncbi(
    row: dict,
    *,
    cache_dir,
    force_refresh: bool = False,
    allow_species_only: bool = False,
) -> Optional[Dict[str, Any]]:
    species = str(row.get("species") or "").strip()
    if not species:
        return None
    candidates = search_genomes_by_taxon(
        species, cache_dir=cache_dir, force_refresh=force_refresh, limit=12
    )
    candidates = [c for c in candidates if c.get("genome_accession")]
    best, score, note = pick_best_candidate(row, candidates, allow_species_only=allow_species_only)
    if not best or score <= 0 or not best.get("genome_accession"):
        return None
    out = dict(best)
    out["match_confidence"] = round(score, 3)
    out["match_notes"] = note
    return out


def fetch_assembly_details(
    accession: str,
    *,
    cache_dir,
    force_refresh: bool = False,
) -> Optional[Dict[str, Any]]:
    """Fetch richer stats for a known accession."""
    acc = str(accession).strip()
    if not acc:
        return None
    url = f"{NCBI_DATASETS_BASE}/genome/accession/{quote(acc)}/dataset_report"
    body, _ = fetch_json(url, cache_dir=cache_dir, force_refresh=force_refresh)
    if not isinstance(body, dict):
        return None
    reports = body.get("reports") or []
    if not reports:
        return None
    return _report_to_candidate(reports[0])
