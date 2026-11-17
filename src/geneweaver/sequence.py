"""
Sequence processing and validation utilities for GeneWeaver.
"""

from typing import Set

VALID_DNA_BASES: Set[str] = {"A", "C", "G", "T"}


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

