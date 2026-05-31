"""BacDive API v2 client for genome metadata."""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional
from urllib.parse import quote

from genome_investigation.api_cache import fetch_json
from genome_investigation.match_scoring import score_bacdive_direct

BACDIVE_API_BASE = "https://api.bacdive.dsmz.de"


def _first_dict_value(d: Any) -> Optional[dict]:
    if isinstance(d, dict):
        if "INSDC accession" in d or "insdc_acc" in d or "accession" in d:
            return d
        for v in d.values():
            if isinstance(v, dict):
                found = _first_dict_value(v)
                if found:
                    return found
    elif isinstance(d, list):
        for item in d:
            found = _first_dict_value(item)
            if found:
                return found
    return None


def _collect_genome_records(seq_info: Any) -> List[dict]:
    records: List[dict] = []
    if not isinstance(seq_info, dict):
        return records
    genomes = seq_info.get("Genome sequences") or seq_info.get("genome sequences") or seq_info.get(
        "sequence_genomes"
    )
    if isinstance(genomes, dict):
        if any(k in genomes for k in ("INSDC accession", "insdc_acc", "assembly level")):
            records.append(genomes)
        else:
            for v in genomes.values():
                if isinstance(v, dict):
                    records.append(v)
    elif isinstance(genomes, list):
        records.extend(x for x in genomes if isinstance(x, dict))
    return records


def _gc_percent(seq_info: Any) -> Optional[float]:
    if not isinstance(seq_info, dict):
        return None
    gc = seq_info.get("GC content") or seq_info.get("GC_content")
    if isinstance(gc, list) and gc:
        gc = gc[0]
    if isinstance(gc, dict):
        val = gc.get("GC-content") or gc.get("GC_content")
        if val is not None:
            try:
                return float(str(val).replace("%", "").strip())
            except ValueError:
                return None
    return None


def parse_bacdive_strain(payload: dict, bacid: str) -> Optional[Dict[str, Any]]:
    """Extract genome fields from BacDive v2 fetch response for one strain."""
    results = payload.get("results") if isinstance(payload, dict) else None
    if not isinstance(results, dict):
        return None
    strain = results.get(str(bacid)) or results.get(bacid)
    if not isinstance(strain, dict):
        return None

    seq_info = strain.get("Sequence information") or strain.get("sequence_information")
    genomes = _collect_genome_records(seq_info)
    if not genomes:
        return None

    # Prefer highest score if present
    def _score(g: dict) -> float:
        try:
            return float(g.get("score", 0) or 0)
        except (TypeError, ValueError):
            return 0.0

    genomes.sort(key=_score, reverse=True)
    g = genomes[0]
    acc = g.get("INSDC accession") or g.get("insdc_acc") or g.get("accession")
    if not acc:
        return None

    conf, notes = score_bacdive_direct()
    tax_block = strain.get("General", {}).get("NCBI tax id") or strain.get("General", {}).get("NCBI tax ID")
    taxid = None
    if isinstance(tax_block, dict):
        taxid = tax_block.get("NCBI tax id") or tax_block.get("NCBI tax ID")

    name_tax = strain.get("Name and taxonomic classification", {})
    order = name_tax.get("order") if isinstance(name_tax, dict) else None
    strain_des = name_tax.get("strain designation") if isinstance(name_tax, dict) else None

    return {
        "genome_accession": str(acc).strip(),
        "assembly_level": str(g.get("assembly level") or g.get("assembly_lvl") or "").strip(),
        "genome_size_bp": None,
        "gc_percent": _gc_percent(seq_info),
        "gene_count": None,
        "cds_count": None,
        "ncbi_taxid": taxid,
        "source_database": "bacdive",
        "match_confidence": conf,
        "match_notes": notes,
        "strain": strain_des,
        "order": order,
    }


def fetch_strain(bacid: str, *, cache_dir, force_refresh: bool = False) -> Optional[dict]:
    bacid_clean = re.sub(r"\.0$", "", str(bacid).strip())
    if not bacid_clean or bacid_clean.lower() == "nan":
        return None
    url = f"{BACDIVE_API_BASE}/v2/fetch/{quote(bacid_clean)}"
    body, _ = fetch_json(url, cache_dir=cache_dir, force_refresh=force_refresh)
    if not isinstance(body, dict) or body.get("error"):
        return None
    return parse_bacdive_strain(body, bacid_clean)
