
"""
Unit tests for GeneWeaver sequence alignment.
Tests the find_off_targets() CPU baseline implementation.
"""

import sys
from pathlib import Path

import pytest

# Add src directory to Python path
src_path = Path(__file__).resolve().parent.parent / "src"

if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from project.alignment import find_off_targets


class TestFindOffTargets:

    def test_exact_match(self):
        """Detect an exact DNA sequence match."""

        genome = "AACCGGTTAACCGG"
        guide = "CCGG"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert len(results) == 2
        assert results[0]["position"] == 2
        assert results[0]["sequence"] == "CCGG"
        assert results[0]["mismatches"] == 0

    def test_single_mismatch(self):
        """Detect a sequence containing one mismatch."""

        genome = "AACCGGTT"
        guide = "CCGA"

        results = find_off_targets(genome, guide, max_mismatches=1)

        assert len(results) == 1
        assert results[0]["position"] == 2
        assert results[0]["mismatches"] == 1

    def test_mismatch_threshold(self):
        """Reject matches exceeding the allowed mismatch count."""

        genome = "AACCGGTT"
        guide = "CCGA"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert results == []

    def test_rna_to_dna_conversion(self):
        """Verify U in guide RNA is converted into T."""

        genome = "AATCGA"
        guide = "AUCG"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert len(results) == 1
        assert results[0]["sequence"] == "ATCG"
        assert results[0]["mismatches"] == 0

    def test_case_insensitive_matching(self):
        """Verify lowercase and uppercase sequences match."""

        genome = "aatcga"
        guide = "ATCG"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert len(results) == 1
        assert results[0]["sequence"] == "ATCG"

    def test_overlapping_matches(self):
        """Verify overlapping sequence matches are detected."""

        genome = "AAAAA"
        guide = "AAA"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert len(results) == 3
        assert [r["position"] for r in results] == [0, 1, 2]

    def test_no_matching_sequence(self):
        """Return an empty list when no valid match exists."""

        genome = "AAAAAAAA"
        guide = "CCCC"

        results = find_off_targets(genome, guide, max_mismatches=0)

        assert results == []

    def test_result_structure(self):
        """Verify the returned dictionary structure."""

        results = find_off_targets(
            "AACCGGTT",
            "CCGG",
            max_mismatches=0
        )

        assert len(results) > 0

        for result in results:
            assert "position" in result
            assert "sequence" in result
            assert "mismatches" in result

            assert isinstance(result["position"], int)
            assert isinstance(result["sequence"], str)
            assert isinstance(result["mismatches"], int)
