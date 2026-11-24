
from types import SimpleNamespace

import pytest

from src.project.pam_alignment import (
    count_mismatches,
    align_pam_candidates
)


def test_count_mismatches():
    assert count_mismatches(
        "ATCG",
        "ATGG"
    ) == 1


def test_unequal_sequence_length():
    with pytest.raises(ValueError):
        count_mismatches("ATCG", "ATC")


def test_pam_alignment():
    candidates = [
        SimpleNamespace(
            sequence_id="chr1",
            position=100,
            strand="+",
            protospacer="ATCG",
            pam="AGG"
        ),
        SimpleNamespace(
            sequence_id="chr1",
            position=200,
            strand="-",
            protospacer="ATGG",
            pam="TGG"
        ),
        SimpleNamespace(
            sequence_id="chr1",
            position=300,
            strand="+",
            protospacer="CCCC",
            pam="AGG"
        )
    ]

    results = align_pam_candidates(
        "ATCG",
        candidates,
        max_mismatches=1
    )

    assert len(results) == 2
    assert results[0].mismatches == 0
    assert results[1].mismatches == 1


def test_mismatch_threshold():
    candidates = [
        SimpleNamespace(
            sequence_id="chr1",
            position=100,
            strand="+",
            protospacer="CCCC",
            pam="AGG"
        )
    ]

    results = align_pam_candidates(
        "ATCG",
        candidates,
        max_mismatches=1
    )

    assert len(results) == 0
