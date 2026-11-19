
import pytest

from project.pam_scanner import scan_pam, Candidate


def test_forward_pam_detection():
    sequence = "A" * 20 + "AGG"

    results = scan_pam("test_sequence", sequence)

    assert len(results) == 1

    hit = results[0]

    assert isinstance(hit, Candidate)
    assert hit.protospacer == "A" * 20
    assert hit.pam == "AGG"
    assert hit.start == 0
    assert hit.end == 20
    assert hit.pam_start == 20
    assert hit.pam_end == 23
    assert hit.strand == "+"


def test_reverse_pam_detection():
    sequence = "CCT" + "T" * 20

    results = scan_pam("test_sequence", sequence)

    reverse_hits = [
        hit for hit in results if hit.strand == "-"
    ]

    assert len(reverse_hits) == 1

    hit = reverse_hits[0]

    assert hit.pam == "AGG"
    assert hit.start == 0
    assert hit.end == 20
    assert hit.pam_start == 0
    assert hit.pam_end == 3
    assert hit.strand == "-"


def test_invalid_sequence():
    with pytest.raises(ValueError):
        scan_pam("test_sequence", "ATGCXYZ")


def test_empty_sequence():
    with pytest.raises(ValueError):
        scan_pam("test_sequence", "")
    