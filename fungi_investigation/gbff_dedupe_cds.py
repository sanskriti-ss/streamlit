"""
Clean fungal GBFF/GBK files so antiSMASH can parse them.

Fixes common RefSeq issues:
  - Multiple CDS features have the same location
  - multiple CDS features have the same name for mapping
  - location contains overlapping exons (CDS, mRNA, and gene)

Usage:
  python -m fungi_investigation.gbff_dedupe_cds \\
    --input path/to/genomic.gbff --output path/to/genomic.dedup.gbff
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple


# antiSMASH-style [start:end] and GenBank start..end
_LOC_RE = re.compile(r"\[(\d+):(\d+)\]|(\d+)\.\.(\d+)")


def _exons_overlap(loc: str) -> bool:
    """True if a join(...) location has overlapping exon intervals."""
    if "join" not in loc.lower():
        return False
    spans = []
    for a, b, c, d in _LOC_RE.findall(loc):
        if a and b:
            start, end = int(a), int(b)
        else:
            start, end = int(c), int(d)
        if start > end:
            start, end = end, start
        spans.append((start, end))
    if len(spans) < 2:
        return False
    spans.sort()
    for i in range(1, len(spans)):
        prev_a, prev_b = spans[i - 1]
        a, b = spans[i]
        # Inclusive interval overlap (excluding exact duplicate spans).
        if a <= prev_b and not (a == prev_a and b == prev_b):
            return True
    return False


def _feature_gene_name(block: List[str]) -> str:
    for bl in block:
        m = re.search(r'/(?:gene|locus_tag|protein_id)="([^"]+)"', bl)
        if m:
            return m.group(1)
    return ""


def _rewrite_gene_name(block: List[str], new_name: str) -> List[str]:
    out: List[str] = []
    replaced = False
    for bl in block:
        if not replaced and re.search(r'/(?:gene|locus_tag)="', bl):
            out.append(re.sub(r'/(gene|locus_tag)="[^"]+"', rf'/\1="{new_name}"', bl, count=1))
            replaced = True
        else:
            out.append(bl)
    return out


def _full_location(key_line_loc: str, block: List[str]) -> str:
    loc_full = key_line_loc
    for bl in block[1:]:
        stripped = bl.strip()
        if stripped.startswith("/"):
            break
        loc_full += stripped
    return loc_full.replace(" ", "")


def dedupe_cds_gbff(text: str) -> Tuple[str, int]:
    """Clean CDS features across all LOCUS records. Returns (text, n_changes)."""
    parts = re.split(r"(?=^LOCUS )", text, flags=re.M)
    out_parts: List[str] = []
    changes = 0
    for part in parts:
        if not part.strip():
            continue
        cleaned, n = _clean_record(part)
        out_parts.append(cleaned)
        changes += n
    return "".join(out_parts), changes


def _clean_record(record: str) -> Tuple[str, int]:
    lines = record.splitlines(keepends=True)
    if not lines:
        return record, 0

    out: List[str] = []
    i = 0
    seen_cds_locs: set[str] = set()
    seen_names: set[str] = set()
    changes = 0

    while i < len(lines):
        line = lines[i]
        m = re.match(r"^ {5}(\S+)\s+(\S.*)$", line)
        if not m:
            out.append(line)
            i += 1
            continue

        key, loc = m.group(1), m.group(2).strip()
        block = [line]
        i += 1
        while i < len(lines):
            nxt = lines[i]
            if re.match(r"^ {5}\S", nxt) or re.match(r"^[A-Z]{2,}", nxt):
                break
            block.append(nxt)
            i += 1

        loc_full = _full_location(loc, block)

        # antiSMASH also rejects overlapping exons on mRNA/gene features.
        if key in {"CDS", "mRNA", "gene"} and _exons_overlap(loc_full):
            changes += 1
            continue

        if key == "CDS":
            if loc_full in seen_cds_locs:
                changes += 1
                continue
            seen_cds_locs.add(loc_full)

            name = _feature_gene_name(block)
            if name:
                if name in seen_names:
                    n = 2
                    candidate = f"{name}_{n}"
                    while candidate in seen_names:
                        n += 1
                        candidate = f"{name}_{n}"
                    block = _rewrite_gene_name(block, candidate)
                    seen_names.add(candidate)
                    changes += 1
                else:
                    seen_names.add(name)

        out.extend(block)

    return "".join(out), changes


def write_deduped_gbff(input_path: Path, output_path: Path) -> int:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    cleaned, changes = dedupe_cds_gbff(text)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(cleaned, encoding="utf-8")
    return changes


def main() -> int:
    parser = argparse.ArgumentParser(description="Clean CDS features in a fungal GBFF file")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    n = write_deduped_gbff(args.input, args.output)
    print(f"[ok] Wrote {args.output} ({n} CDS fixes: drop/rename)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
