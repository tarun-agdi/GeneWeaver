"""
Unit tests for forward-strand SpCas9 NGG PAM scanning (find_forward_pams).
"""

from pathlib import Path
import sys
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.scanner import find_forward_pams, scan_forward_pams


class TestFindForwardPams:
    """Tests for forward-strand SpCas9 PAM scanning."""

    def test_one_valid_pam(self):
        # Single TGG PAM at position 5
        seq = "ATCGATGGAACCT"
        results = find_forward_pams(seq)
        assert results == [{"pam": "TGG", "position": 5}]

    def test_example_sequence_all_matching_pams(self):
        # ATCGATGGAACCTAGG has TGG at position 5 and AGG at position 13
        seq = "ATCGATGGAACCTAGG"
        results = find_forward_pams(seq)
        assert results == [
            {"pam": "TGG", "position": 5},
            {"pam": "AGG", "position": 13},
        ]

    def test_multiple_pams_across_sequence(self):
        # TGG at 0, CGG at 4, AGG at 8
        seq = "TGGACGGTAGGT"
        results = find_forward_pams(seq)
        assert len(results) == 3
        assert results[0] == {"pam": "TGG", "position": 0}
        assert results[1] == {"pam": "CGG", "position": 4}
        assert results[2] == {"pam": "AGG", "position": 8}

    def test_overlapping_pams(self):
        # AGGG contains AGG at index 0 and GGG at index 1
        seq = "AGGG"
        results = find_forward_pams(seq)
        assert results == [
            {"pam": "AGG", "position": 0},
            {"pam": "GGG", "position": 1},
        ]

    def test_no_pam_present(self):
        # Sequence containing no NGG pattern
        assert find_forward_pams("ATATATATAT") == []
        assert find_forward_pams("ACACACACAC") == []
        assert find_forward_pams("TTTTTTTTTT") == []

    def test_sequence_shorter_than_pam(self):
        # Sequences shorter than 3 nt cannot contain a PAM
        assert find_forward_pams("A") == []
        assert find_forward_pams("AG") == []
        assert find_forward_pams("GG") == []

    def test_pam_at_beginning_of_sequence(self):
        # TGG at index 0
        seq = "TGGATCGATCG"
        results = find_forward_pams(seq)
        assert len(results) == 1
        assert results[0] == {"pam": "TGG", "position": 0}

    def test_pam_at_end_of_sequence(self):
        # AGG at the very end
        seq = "ATCGATCAGG"
        results = find_forward_pams(seq)
        assert len(results) == 1
        assert results[0] == {"pam": "AGG", "position": 7}

    def test_all_four_canonical_ngg_variants(self):
        for pam in ("AGG", "CGG", "GGG", "TGG"):
            results = find_forward_pams(f"{pam}AAAA")
            assert len(results) == 1
            assert results[0] == {"pam": pam, "position": 0}

    def test_non_ngg_patterns_ignored(self):
        for non_pam in ("AGA", "AGC", "AGT", "TGA", "TGC", "CGA", "CGC", "ACT"):
            results = find_forward_pams(f"AA{non_pam}AA")
            assert results == []

    def test_lowercase_and_rna_uracil_normalized(self):
        # Lowercase and RNA 'u' converted according to existing validation rules
        seq_rna_lower = "atcgauggaacct"
        results = find_forward_pams(seq_rna_lower)
        assert results == [{"pam": "TGG", "position": 5}]

    def test_invalid_characters_raise_value_error(self):
        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            find_forward_pams("ATCGATGGNAACCT")

        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            find_forward_pams("ATCG-TGG")

    def test_empty_sequence_raises_value_error(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            find_forward_pams("")

    def test_non_string_sequence_raises_type_error(self):
        with pytest.raises(TypeError, match="Sequence must be a string"):
            find_forward_pams(12345)

        with pytest.raises(TypeError, match="Sequence must be a string"):
            find_forward_pams(None)

    def test_scan_forward_pams_alias(self):
        seq = "ATCGATGGAACCT"
        assert scan_forward_pams(seq) == find_forward_pams(seq)
