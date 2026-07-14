"""
Remove duplicate CDS features that share the same location from a GBFF/GBK file.

antiSMASH refuses some RefSeq fungal GBFF files with:
  Multiple CDS features have the same location: join{...}

Usage:
  python -m fungi_investigation.gbff_dedupe_cds \\
    --input path/to/genomic.gbff --output path/to/genomic.dedup.gbff
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def dedupe_cds_gbff(text: str) -> tuple[str, int]:
    """Drop CDS feature blocks that share the same location string within a record."""
    parts = re.split(r"(?=^LOCUS )", text, flags=re.M)
    out_parts: list[str] = []
    removed = 0
    for part in parts:
        if not part.strip():
            continue
        cleaned, n = _dedupe_record(part)
        out_parts.append(cleaned)
        removed += n
    return "".join(out_parts), removed


def _dedupe_record(record: str) -> tuple[str, int]:
    lines = record.splitlines(keepends=True)
    if not lines:
        return record, 0

    out: list[str] = []
    i = 0
    seen_cds_locs: set[str] = set()
    removed = 0
    while i < len(lines):
        line = lines[i]
        # Feature key lines start at column 5 (0-index: 5 spaces) per GenBank.
        m = re.match(r"^ {5}(\S+)\s+(\S.*)$", line)
        if not m:
            out.append(line)
            i += 1
            continue

        key, loc = m.group(1), m.group(2).strip()
        # Collect full feature block (qualifier lines start with 21 spaces).
        block = [line]
        i += 1
        while i < len(lines):
            nxt = lines[i]
            if re.match(r"^ {5}\S", nxt) or re.match(r"^[A-Z]{2,}", nxt):
                break
            block.append(nxt)
            i += 1

        if key == "CDS":
            # Continuation of location on following lines without '/'
            loc_full = loc
            for bl in block[1:]:
                stripped = bl.strip()
                if stripped.startswith("/"):
                    break
                loc_full += stripped
            if loc_full in seen_cds_locs:
                removed += 1
                continue
            seen_cds_locs.add(loc_full)

        out.extend(block)

    return "".join(out), removed


def write_deduped_gbff(input_path: Path, output_path: Path) -> int:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    cleaned, removed = dedupe_cds_gbff(text)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(cleaned, encoding="utf-8")
    return removed


def main() -> int:
    parser = argparse.ArgumentParser(description="Dedupe CDS features in a GBFF file")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    removed = write_deduped_gbff(args.input, args.output)
    print(f"[ok] Wrote {args.output} (removed {removed} duplicate CDS features)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
