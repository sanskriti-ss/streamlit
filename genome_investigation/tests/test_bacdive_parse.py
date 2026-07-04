"""Tests for BacDive response parsing."""

from genome_investigation.bacdive_client import parse_bacdive_strain

SAMPLE = {
    "count": 1,
    "results": {
        "158570": {
            "General": {"NCBI tax id": {"NCBI tax id": 2044859}},
            "Name and taxonomic classification": {
                "species": "Cohnella algarum",
                "strain designation": "Pch-40",
                "order": "Caryophanales",
            },
            "Sequence information": {
                "Genome sequences": {
                    "assembly level": "contig",
                    "INSDC accession": "GCA_016937515",
                    "score": 68.5,
                },
                "GC content": [{"GC-content": "55.6"}],
            },
        }
    },
}


def test_parse_bacdive_strain():
    parsed = parse_bacdive_strain(SAMPLE, "158570")
    assert parsed is not None
    assert parsed["genome_accession"] == "GCA_016937515"
    assert parsed["match_confidence"] == 1.0
    assert parsed["gc_percent"] == 55.6
