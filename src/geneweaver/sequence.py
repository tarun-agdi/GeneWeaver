"""
Sequence processing and validation utilities for GeneWeaver.
"""

from pathlib import Path
from typing import List, Set, Tuple, Union

VALID_DNA_BASES: Set[str] = {"A", "C", "G", "T"}
_DNA_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")

__all__ = [
    "VALID_DNA_BASES",
    "validate_dna_sequence",
    "reverse_complement",
    "read_fasta",
]


def validate_dna_sequence(sequence: str) -> str:
    """Validate and normalize a DNA or RNA sequence.

    Normalizes the input by converting lowercase characters to uppercase
    and replacing RNA uracil ('U'/'u') with DNA thymine ('T').
    Rejects any characters outside A, C, G, and T.

    Args:
        sequence: Nucleotide sequence string to validate.

    Returns:
        Normalized DNA sequence string containing only uppercase 'A', 'C', 'G', 'T'.

    Raises:
        TypeError: If sequence is not a string.
        ValueError: If sequence is empty or contains invalid characters.
    """
    if not isinstance(sequence, str):
        raise TypeError(f"Sequence must be a string, got {type(sequence).__name__}")

    if not sequence:
        raise ValueError("Sequence cannot be empty")

    normalized = sequence.upper().replace("U", "T")

    invalid_chars = sorted(set(normalized) - VALID_DNA_BASES)
    if invalid_chars:
        formatted = ", ".join(repr(c) for c in invalid_chars)
        raise ValueError(
            f"Invalid character(s) found in DNA sequence: {formatted}. "
            "Only A, C, G, T (and RNA U) are allowed."
        )

    return normalized


def reverse_complement(sequence: str) -> str:
    """Compute the reverse complement of a DNA or RNA sequence.

    Validates and normalizes the sequence using `validate_dna_sequence`,
    then returns the reverse complement sequence in 5' to 3' orientation
    (A <-> T, C <-> G).

    Args:
        sequence: Nucleotide sequence string (DNA or RNA).

    Returns:
        The normalized uppercase reverse-complement DNA sequence.

    Raises:
        TypeError: If sequence is not a string.
        ValueError: If sequence is empty or contains invalid characters.
    """
    normalized = validate_dna_sequence(sequence)
    return normalized.translate(_DNA_COMPLEMENT_TABLE)[::-1]


def read_fasta(file_path: Union[str, Path]) -> List[Tuple[str, str]]:
    """Read and validate sequences from a standard FASTA file.

    Parses single or multi-sequence FASTA files, concatenates multiline
    sequences, and validates each sequence using `validate_dna_sequence`.

    Args:
        file_path: Path to the FASTA file (str or pathlib.Path).

    Returns:
        A list of (sequence_id, validated_sequence) tuples preserving file order.

    Raises:
        TypeError: If file_path is not a str or Path.
        FileNotFoundError: If the file does not exist.
        ValueError: If the path is a directory, if the file is empty,
            or if the FASTA formatting or sequence data is malformed.
    """
    if not isinstance(file_path, (str, Path)):
        raise TypeError(f"file_path must be a string or Path, got {type(file_path).__name__}")

    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    if path.is_dir():
        raise ValueError(f"Expected a FASTA file path, but got a directory: {file_path}")

    records: List[Tuple[str, str]] = []
    current_header: Union[str, None] = None
    current_seq_lines: List[str] = []
    has_content = False
    line_num = 0

    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line_num += 1
            line = raw_line.strip()
            if not line:
                continue

            has_content = True

            if line.startswith(">"):
                if current_header is not None:
                    combined_seq = "".join(current_seq_lines)
                    if not combined_seq:
                        raise ValueError(
                            f"Malformed FASTA: Header '{current_header}' has no sequence data"
                        )
                    validated = validate_dna_sequence(combined_seq)
                    records.append((current_header, validated))
                    current_seq_lines = []

                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Malformed FASTA: Header at line {line_num} has empty sequence identifier"
                    )
                current_header = header
            else:
                if current_header is None:
                    raise ValueError(
                        f"Malformed FASTA: Expected '>' header line at line {line_num}, got sequence line: {line!r}"
                    )
                current_seq_lines.append(line)

    if not has_content:
        raise ValueError(f"FASTA file is empty: {file_path}")

    if current_header is not None:
        combined_seq = "".join(current_seq_lines)
        if not combined_seq:
            raise ValueError(
                f"Malformed FASTA: Header '{current_header}' has no sequence data"
            )
        validated = validate_dna_sequence(combined_seq)
        records.append((current_header, validated))

    return records
