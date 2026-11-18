"""
Unit tests for DNA sequence validation, reverse complement, and FASTA reading utilities.
"""

from pathlib import Path
import sys
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.sequence import (
    read_fasta,
    reverse_complement,
    validate_dna_sequence,
)


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


class TestReverseComplement:
    """Tests for reverse_complement()."""

    def test_canonical_bases_correctness(self):
        # A <-> T, C <-> G reversed
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("CCCC") == "GGGG"
        assert reverse_complement("GATTACA") == "TGTAATC"
        assert reverse_complement("AACCGGTT") == "AACCGGTT"

    def test_palindromic_sequences(self):
        # EcoRI restriction site (GAATTC)
        assert reverse_complement("GAATTC") == "GAATTC"
        # HindIII restriction site (AAGCTT)
        assert reverse_complement("AAGCTT") == "AAGCTT"
        # Small palindromes
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("GCGC") == "GCGC"
        assert reverse_complement("ATAT") == "ATAT"

    def test_lowercase_input(self):
        assert reverse_complement("gattaca") == "TGTAATC"
        assert reverse_complement("gaattc") == "GAATTC"
        assert reverse_complement("aCgT") == "ACGT"

    def test_rna_input(self):
        # Uracil should be treated as thymine: ACGU -> ACGT -> reverse complement is ACGT
        assert reverse_complement("ACGU") == "ACGT"
        assert reverse_complement("uuu") == "AAA"

    def test_empty_input_raises_value_error(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            reverse_complement("")

    def test_invalid_dna_input_raises_value_error(self):
        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence"):
            reverse_complement("ACGTN")

        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence"):
            reverse_complement("GATTACA-1")

    def test_non_string_type_raises_type_error(self):
        with pytest.raises(TypeError, match="Sequence must be a string"):
            reverse_complement(None)

        with pytest.raises(TypeError, match="Sequence must be a string"):
            reverse_complement(12345)


class TestReadFasta:
    """Tests for read_fasta()."""

    def test_single_sequence_fasta(self, tmp_path: Path):
        fasta_file = tmp_path / "single.fasta"
        fasta_file.write_text(">seq1\nACGTACGTACGT\n", encoding="utf-8")

        records = read_fasta(fasta_file)
        assert len(records) == 1
        assert records[0] == ("seq1", "ACGTACGTACGT")

    def test_multiple_sequence_fasta(self, tmp_path: Path):
        fasta_file = tmp_path / "multiple.fasta"
        content = (
            ">target_1\n"
            "ACGTACGT\n"
            ">target_2\n"
            "TGCATGCA\n"
            ">target_3\n"
            "GGCCGGCC\n"
        )
        fasta_file.write_text(content, encoding="utf-8")

        records = read_fasta(fasta_file)
        assert len(records) == 3
        assert records[0] == ("target_1", "ACGTACGT")
        assert records[1] == ("target_2", "TGCATGCA")
        assert records[2] == ("target_3", "GGCCGGCC")

    def test_multiline_sequence(self, tmp_path: Path):
        fasta_file = tmp_path / "multiline.fasta"
        content = (
            ">long_seq\n"
            "ACGT\n"
            "ACGT\n"
            "TGCA\n"
            "TGCA\n"
        )
        fasta_file.write_text(content, encoding="utf-8")

        records = read_fasta(fasta_file)
        assert len(records) == 1
        assert records[0] == ("long_seq", "ACGTACGTTGCATGCA")

    def test_preserves_sequence_identifiers_with_spaces_or_special_chars(self, tmp_path: Path):
        fasta_file = tmp_path / "annotated.fasta"
        content = (
            ">NC_000913.3 Escherichia coli str. K-12 substr. MG1655\n"
            "AGCT\n"
        )
        fasta_file.write_text(content, encoding="utf-8")

        records = read_fasta(fasta_file)
        assert len(records) == 1
        assert records[0][0] == "NC_000913.3 Escherichia coli str. K-12 substr. MG1655"
        assert records[0][1] == "AGCT"

    def test_normalizes_lowercase_and_rna(self, tmp_path: Path):
        fasta_file = tmp_path / "lowercase.fasta"
        fasta_file.write_text(">seq_rna\nacguacgu\n", encoding="utf-8")

        records = read_fasta(fasta_file)
        assert records[0] == ("seq_rna", "ACGTACGT")

    def test_accepts_string_path(self, tmp_path: Path):
        fasta_file = tmp_path / "str_path.fasta"
        fasta_file.write_text(">seq1\nACGT\n", encoding="utf-8")

        records = read_fasta(str(fasta_file))
        assert len(records) == 1
        assert records[0] == ("seq1", "ACGT")

    def test_empty_file_raises_value_error(self, tmp_path: Path):
        empty_file = tmp_path / "empty.fasta"
        empty_file.write_text("", encoding="utf-8")

        with pytest.raises(ValueError, match="FASTA file is empty"):
            read_fasta(empty_file)

    def test_whitespace_only_file_raises_value_error(self, tmp_path: Path):
        ws_file = tmp_path / "whitespace.fasta"
        ws_file.write_text("   \n\n\t  \n", encoding="utf-8")

        with pytest.raises(ValueError, match="FASTA file is empty"):
            read_fasta(ws_file)

    def test_missing_header_raises_value_error(self, tmp_path: Path):
        malformed = tmp_path / "no_header.fasta"
        malformed.write_text("ACGTACGT\n", encoding="utf-8")

        with pytest.raises(ValueError, match=r"Expected '>' header line"):
            read_fasta(malformed)

    def test_empty_header_identifier_raises_value_error(self, tmp_path: Path):
        malformed = tmp_path / "empty_header.fasta"
        malformed.write_text(">\nACGT\n", encoding="utf-8")

        with pytest.raises(ValueError, match="empty sequence identifier"):
            read_fasta(malformed)

    def test_header_without_sequence_before_next_header(self, tmp_path: Path):
        malformed = tmp_path / "missing_seq.fasta"
        malformed.write_text(">seq1\n>seq2\nACGT\n", encoding="utf-8")

        with pytest.raises(ValueError, match="Header 'seq1' has no sequence data"):
            read_fasta(malformed)

    def test_header_without_sequence_at_eof(self, tmp_path: Path):
        malformed = tmp_path / "eof_missing_seq.fasta"
        malformed.write_text(">seq1\nACGT\n>seq2\n", encoding="utf-8")

        with pytest.raises(ValueError, match="Header 'seq2' has no sequence data"):
            read_fasta(malformed)

    def test_invalid_characters_in_fasta_sequence(self, tmp_path: Path):
        malformed = tmp_path / "invalid_dna.fasta"
        malformed.write_text(">seq1\nACGTN\n", encoding="utf-8")

        with pytest.raises(ValueError, match=r"Invalid character\(s\) found in DNA sequence: 'N'"):
            read_fasta(malformed)

    def test_file_not_found_raises_error(self, tmp_path: Path):
        non_existent = tmp_path / "does_not_exist.fasta"
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            read_fasta(non_existent)

    def test_directory_path_raises_value_error(self, tmp_path: Path):
        with pytest.raises(ValueError, match="Expected a FASTA file path, but got a directory"):
            read_fasta(tmp_path)

    def test_invalid_path_type_raises_type_error(self):
        with pytest.raises(TypeError, match="file_path must be a string or Path"):
            read_fasta(12345)
