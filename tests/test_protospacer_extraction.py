"""
Unit tests for 20-nt protospacer extraction on the forward strand (extract_forward_protospacers).
"""

from pathlib import Path
import sys
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.scanner import extract_forward_candidates, extract_forward_protospacers


class TestExtractForwardProtospacers:
    """Tests for forward-strand 20-nt protospacer extraction."""

    def test_valid_pam_with_20_preceding_nucleotides(self):
        protospacer = "ACGTACGTACGTACGTACGT"  # Exactly 20 nt
        pam = "AGG"
        seq = protospacer + pam  # Exactly 23 nt
        results = extract_forward_protospacers(seq)

        assert len(results) == 1
        expected = {
            "protospacer": protospacer,
            "pam": pam,
            "start": 0,
            "end": 20,
            "pam_position": 20,
            "strand": "+",
        }
        assert results[0] == expected

    def test_multiple_pam_candidates(self):
        protospacer1 = "AAAAAAAAAAAAAAAAAAAA"  # 20 A's
        pam1 = "AGG"
        protospacer2 = "CCCCCCCCCCCCCCCCCCCC"  # 20 C's
        pam2 = "CGG"
        seq = protospacer1 + pam1 + protospacer2 + pam2  # 20 + 3 + 20 + 3 = 46 nt

        results = extract_forward_protospacers(seq)
        assert len(results) == 2

        # First candidate at start 0, end 20, PAM at 20
        assert results[0] == {
            "protospacer": protospacer1,
            "pam": pam1,
            "start": 0,
            "end": 20,
            "pam_position": 20,
            "strand": "+",
        }

        # Second candidate starts at 23, end 43, PAM at 43
        assert results[1] == {
            "protospacer": protospacer2,
            "pam": pam2,
            "start": 23,
            "end": 43,
            "pam_position": 43,
            "strand": "+",
        }

    def test_pam_with_fewer_than_20_preceding_nucleotides_skipped(self):
        # PAM preceded by only 10 nucleotides
        seq = "AAAAAAAAAA" + "TGG"  # 10 + 3 = 13 nt
        results = extract_forward_protospacers(seq)
        assert results == []

    def test_pam_with_19_preceding_nucleotides_skipped(self):
        # Edge case: exactly 1 nucleotide short of a complete 20-nt protospacer
        seq = "A" * 19 + "TGG"
        results = extract_forward_protospacers(seq)
        assert results == []

    def test_pam_at_beginning_of_sequence_skipped(self):
        # PAM starting at index 0 has 0 preceding nucleotides
        seq = "TGG" + "A" * 30
        results = extract_forward_protospacers(seq)
        assert results == []

    def test_correct_start_and_end_coordinates_embedded(self):
        prefix = "TTTTTT"  # 6 nt
        protospacer = "GTCGATCGATCGATCGATCA"  # 20 nt
        pam = "TGG"  # 3 nt
        suffix = "AAAAAA"  # 6 nt
        seq = prefix + protospacer + pam + suffix

        results = extract_forward_protospacers(seq)
        assert len(results) == 1

        c = results[0]
        assert c["start"] == 6
        assert c["end"] == 26
        assert c["pam_position"] == 26
        assert c["strand"] == "+"
        assert c["protospacer"] == protospacer
        assert c["pam"] == pam
        # Coordinate slice parity check
        assert seq[c["start"] : c["end"]] == protospacer
        assert seq[c["pam_position"] : c["pam_position"] + 3] == pam

    def test_overlapping_pams_extraction(self):
        # Protospacer + AGGG (contains AGG at 20 and GGG at 21)
        protospacer = "ACGTACGTACGTACGTACGT"  # 20 nt
        seq = protospacer + "AGGG"  # 24 nt
        results = extract_forward_protospacers(seq)

        assert len(results) == 2
        # First PAM AGG at 20 -> protospacer seq[0:20]
        assert results[0]["start"] == 0
        assert results[0]["end"] == 20
        assert results[0]["pam_position"] == 20
        assert results[0]["pam"] == "AGG"
        assert results[0]["protospacer"] == protospacer

        # Second PAM GGG at 21 -> protospacer seq[1:21]
        assert results[1]["start"] == 1
        assert results[1]["end"] == 21
        assert results[1]["pam_position"] == 21
        assert results[1]["pam"] == "GGG"
        assert results[1]["protospacer"] == seq[1:21]

    def test_case_and_rna_uracil_normalization(self):
        # Lowercase and RNA uracil should be normalized
        protospacer_rna = "acguacguacguacguacgu"
        seq = protospacer_rna + "ugg"
        results = extract_forward_protospacers(seq)

        assert len(results) == 1
        assert results[0]["protospacer"] == "ACGTACGTACGTACGTACGT"
        assert results[0]["pam"] == "TGG"

    def test_invalid_input_rejection(self):
        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            extract_forward_protospacers("A" * 20 + "NGG")

        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            extract_forward_protospacers("")

        with pytest.raises(TypeError, match="Sequence must be a string"):
            extract_forward_protospacers(None)

        with pytest.raises(ValueError, match="protospacer_len must be a positive integer"):
            extract_forward_protospacers("A" * 20 + "TGG", protospacer_len=0)

    def test_alias_equality(self):
        seq = "A" * 20 + "TGG"
        assert extract_forward_protospacers(seq) == extract_forward_candidates(seq)
