"""
PAM scanning and protospacer extraction utilities for GeneWeaver.
"""

from typing import List

from geneweaver.models import CandidateResult
from geneweaver.sequence import VALID_DNA_BASES, reverse_complement, validate_dna_sequence

__all__ = [
    "find_forward_pams",
    "scan_forward_pams",
    "extract_forward_protospacers",
    "extract_forward_candidates",
    "find_pam_candidates",
    "find_spcas9_candidates",
    "detect_pam_candidates",
    "scan_pam_candidates",
]


def _is_spcas9_pam(pam_window: str) -> bool:
    """Check if a 3-nucleotide window matches the SpCas9 NGG PAM pattern."""
    return (
        len(pam_window) == 3
        and pam_window[1:] == "GG"
        and pam_window[0] in VALID_DNA_BASES
    )


def find_forward_pams(sequence: str) -> List[dict]:
    """Scan the forward strand of a DNA sequence for SpCas9 NGG PAM sites.

    Identifies all occurrences of the NGG pattern on the forward strand,
    where N is any valid canonical DNA base (A, C, G, or T). Reuses
    the project's sequence validation to normalize case and RNA uracil,
    and reject invalid characters.

    Args:
        sequence: The DNA nucleotide sequence string to scan.

    Returns:
        List of dicts with keys:
            - 'pam': The 3-nucleotide PAM sequence string.
            - 'position': The 0-based starting index in the sequence.

    Raises:
        TypeError: If sequence is not a string.
        ValueError: If sequence is empty or contains invalid characters.
    """
    normalized_seq = validate_dna_sequence(sequence)
    results: List[dict] = []

    for i in range(len(normalized_seq) - 2):
        window = normalized_seq[i : i + 3]
        if _is_spcas9_pam(window):
            results.append({"pam": window, "position": i})

    return results


# Convenience alias
scan_forward_pams = find_forward_pams


def extract_forward_protospacers(
    sequence: str, protospacer_len: int = 20
) -> List[dict]:
    """Extract candidate 20-nt protospacers immediately upstream of forward NGG PAMs.

    Identifies forward-strand SpCas9 NGG PAMs and extracts the upstream protospacer.
    Skips any PAM candidate where fewer than `protospacer_len` nucleotides are available upstream.

    Args:
        sequence: The DNA nucleotide sequence string to scan.
        protospacer_len: Length of the protospacer to extract (default 20).

    Returns:
        List of dicts with keys:
            - 'protospacer': The extracted 20-nucleotide protospacer string.
            - 'pam': The 3-nucleotide PAM sequence string.
            - 'start': 0-based start index of the protospacer (inclusive).
            - 'end': 0-based end index of the protospacer (exclusive).
            - 'pam_position': 0-based starting index of the PAM in the sequence.
            - 'strand': Target strand ('+').

    Raises:
        TypeError: If sequence is not a string.
        ValueError: If sequence is empty, contains invalid characters, or if
            protospacer_len is not a positive integer.
    """
    if not isinstance(protospacer_len, int) or protospacer_len <= 0:
        raise ValueError(f"protospacer_len must be a positive integer, got {protospacer_len}")

    normalized_seq = validate_dna_sequence(sequence)
    pams = find_forward_pams(normalized_seq)
    results: List[dict] = []

    for pam_info in pams:
        p = pam_info["position"]
        if p < protospacer_len:
            continue

        start = p - protospacer_len
        end = p
        protospacer = normalized_seq[start:end]

        results.append(
            {
                "protospacer": protospacer,
                "pam": pam_info["pam"],
                "start": start,
                "end": end,
                "pam_position": p,
                "strand": "+",
            }
        )

    return results


# Convenience alias
extract_forward_candidates = extract_forward_protospacers


def find_pam_candidates(
    sequence: str,
    sequence_id: str = "seq",
    pam: str = "NGG",
    protospacer_len: int = 20,
    strands: str = "+-",
    sort_by_coordinate: bool = False,
) -> List[CandidateResult]:
    """Scan a DNA sequence for SpCas9 PAM (NGG) sites and extract candidate protospacers.

    Scans the forward and/or reverse-complement strand of the provided sequence.
    For each valid NGG PAM with a complete adjacent protospacer window, creates
    a CandidateResult record preserving 0-based half-open reference coordinates.

    Args:
        sequence: The DNA sequence to scan (can include RNA U or lowercase; validated).
        sequence_id: Identifier for the sequence (e.g. chromosome or FASTA header).
        pam: The PAM pattern (currently strictly 'NGG' for SpCas9).
        protospacer_len: Length of the adjacent protospacer to extract (default 20).
        strands: Strands to scan: '+', '-', or '+-' (default both).
        sort_by_coordinate: If True, sort results by reference start coordinate.
            If False, returns forward-strand candidates followed by reverse-strand candidates.

    Returns:
        List of CandidateResult instances for all identified complete target sites.

    Raises:
        TypeError: If arguments are of incorrect type.
        ValueError: If sequence or sequence_id is invalid, or if strands is unrecognized.
    """
    if not isinstance(sequence_id, str):
        raise TypeError(f"sequence_id must be a string, got {type(sequence_id).__name__}")
    if not sequence_id.strip():
        raise ValueError("sequence_id cannot be empty or whitespace only")

    if not isinstance(strands, str) or not any(s in strands for s in ("+", "-")):
        raise ValueError(f"strands must contain '+' and/or '-', got {strands!r}")

    if not isinstance(protospacer_len, int) or protospacer_len <= 0:
        raise ValueError(f"protospacer_len must be a positive integer, got {protospacer_len}")

    # Validate and normalize sequence (handles lowercase, U -> T, checks A/C/G/T)
    seq = validate_dna_sequence(sequence)
    seq_len = len(seq)
    pam_len = 3  # SpCas9 NGG
    min_window = protospacer_len + pam_len

    forward_candidates: List[CandidateResult] = []
    reverse_candidates: List[CandidateResult] = []

    # 1. Forward Strand Scanning
    if "+" in strands and seq_len >= min_window:
        for i in range(protospacer_len, seq_len - pam_len + 1):
            pam_window = seq[i : i + pam_len]
            if _is_spcas9_pam(pam_window):
                protospacer = seq[i - protospacer_len : i]
                start = i - protospacer_len
                end = i + pam_len
                forward_candidates.append(
                    CandidateResult(
                        sequence_id=sequence_id,
                        coordinates=(start, end),
                        strand="+",
                        protospacer=protospacer,
                        pam=pam_window,
                        mismatch_count=0,
                        mismatch_positions=(),
                        score=0.0,
                    )
                )

    # 2. Reverse-Complement Strand Scanning
    if "-" in strands and seq_len >= min_window:
        rc_seq = reverse_complement(seq)
        for j in range(protospacer_len, seq_len - pam_len + 1):
            pam_window = rc_seq[j : j + pam_len]
            if _is_spcas9_pam(pam_window):
                protospacer = rc_seq[j - protospacer_len : j]
                j_start = j - protospacer_len
                j_end = j + pam_len
                # Map reverse-strand coordinates back to 0-based half-open reference coordinates
                ref_start = seq_len - j_end
                ref_end = seq_len - j_start
                reverse_candidates.append(
                    CandidateResult(
                        sequence_id=sequence_id,
                        coordinates=(ref_start, ref_end),
                        strand="-",
                        protospacer=protospacer,
                        pam=pam_window,
                        mismatch_count=0,
                        mismatch_positions=(),
                        score=0.0,
                    )
                )
        # Order reverse candidates by reference start coordinate
        reverse_candidates.sort(key=lambda c: c.start)

    all_candidates = forward_candidates + reverse_candidates
    if sort_by_coordinate:
        all_candidates.sort(key=lambda c: (c.start, c.end, c.strand))

    return all_candidates


# Convenience aliases matching assignment terminology
find_spcas9_candidates = find_pam_candidates
detect_pam_candidates = find_pam_candidates
scan_pam_candidates = find_pam_candidates

