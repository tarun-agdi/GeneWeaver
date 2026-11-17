"""
Unit tests for DNA sequence validation utilities.
"""

import sys
from pathlib import Path
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.sequence import validate_dna_sequence


class TestValidateDnaSequence:
    """Tests for validate_dna_sequence()."""

    def test_valid_uppercase_sequence(self):
        assert validate_dna_sequence("ACGT") == "ACGT"
        assert validate_dna_sequence("GATTACA") == "GATTACA"
        assert validate_dna_sequence("A") == "A"

    def test_lowercase_conversion(self):
        assert validate_dna_sequence("acgt") == "ACGT"
        assert validate_dna_sequence("aCgTa") == "ACGTA"
        assert validate_dna_sequence("gattaca") == "GATTACA"

    def test_rna_u_to_t_conversion(self):
        assert validate_dna_sequence("ACGU") == "ACGT"
        assert validate_dna_sequence("acgu") == "ACGT"
        assert validate_dna_sequence("UUU") == "TTT"
        assert validate_dna_sequence("u") == "T"

    def test_typical_crispr_guide_sequence(self):
        # 20-nt typical guide with RNA uracils
        guide_rna = "ggugaguaguugccccugac"
        expected_dna = "GGTGAGTAGTTGCCCCTGAC"
        assert validate_dna_sequence(guide_rna) == expected_dna

    def test_empty_sequence_raises_value_error(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            validate_dna_sequence("")

    def test_invalid_characters_raise_value_error_with_details(self):
        # Ambiguity code 'N'
        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence: 'N'"):
            validate_dna_sequence("ACGTN")

        # Numbers
        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence: '1', '2'"):
            validate_dna_sequence("ACG12T")

        # Whitespace
        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence: ' '"):
            validate_dna_sequence("AC GT")

        # Alignment gap characters
        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence: '-'"):
            validate_dna_sequence("AC-GT")

    def test_non_string_types_raise_type_error(self):
        with pytest.raises(TypeError, match="Sequence must be a string, got NoneType"):
            validate_dna_sequence(None)

        with pytest.raises(TypeError, match="Sequence must be a string, got int"):
            validate_dna_sequence(12345)

        with pytest.raises(TypeError, match="Sequence must be a string, got list"):
            validate_dna_sequence(["A", "C", "G", "T"])

        with pytest.raises(TypeError, match="Sequence must be a string, got bool"):
            validate_dna_sequence(True)

