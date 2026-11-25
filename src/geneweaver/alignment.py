"""
CPU reference alignment and mismatch scoring engine for GeneWeaver.

Provides fixed-length, ungapped nucleotide comparison between a 20-nt CRISPR
guide sequence and candidate 20-nt protospacer sequences.

Note:
    The alignment score computed here (``matches - mismatches``) is a simple
    ungapped sequence-comparison metric for computational search and CPU/GPU
    benchmarking. It is NOT a biological cleavage-risk prediction.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Tuple, Union

from geneweaver.models import CandidateResult
from geneweaver.scanner import find_pam_candidates
from geneweaver.sequence import read_fasta, validate_dna_sequence

DEFAULT_GUIDE_LENGTH = 20

__all__ = [
    "DEFAULT_GUIDE_LENGTH",
    "validate_guide_sequence",
    "align_protospacer_cpu",
    "align_candidate_cpu",
    "align_candidates_cpu",
    "search_and_align_cpu",
    "align_fasta_cpu",
]


def validate_guide_sequence(
    guide: str, expected_length: int = DEFAULT_GUIDE_LENGTH
) -> str:
    """Validate and normalize a CRISPR guide sequence.

    Reuses ``validate_dna_sequence`` to normalize case and convert RNA 'U' to
    DNA 'T', then enforces the fixed guide length (default 20 nt).

    Args:
        guide: Guide RNA or DNA sequence string.
        expected_length: Required nucleotide length (default 20).

    Returns:
        Normalized uppercase DNA guide sequence of length ``expected_length``.

    Raises:
        TypeError: If ``guide`` is not a string.
        ValueError: If ``guide`` is empty, contains invalid bases, or does not
            match ``expected_length``.
    """
    if not isinstance(expected_length, int) or isinstance(expected_length, bool) or expected_length <= 0:
        raise ValueError(f"expected_length must be a positive integer, got {expected_length!r}")

    normalized = validate_dna_sequence(guide)
    if len(normalized) != expected_length:
        raise ValueError(
            f"Guide sequence must be exactly {expected_length} nucleotides, "
            f"got {len(normalized)}"
        )
    return normalized


def _validate_candidate_protospacer(
    protospacer: str, expected_length: int = DEFAULT_GUIDE_LENGTH
) -> str:
    """Validate and normalize a candidate protospacer sequence."""
    normalized = validate_dna_sequence(protospacer)
    if len(normalized) != expected_length:
        raise ValueError(
            f"Candidate protospacer must be exactly {expected_length} nucleotides, "
            f"got {len(normalized)}"
        )
    return normalized


def align_protospacer_cpu(
    guide: str,
    protospacer: str,
    expected_length: int = DEFAULT_GUIDE_LENGTH,
) -> Tuple[int, Tuple[int, ...], float]:
    """Perform fixed-length, ungapped CPU comparison between a guide and protospacer.

    Scoring system:
        - Match: +1
        - Mismatch: -1
        - ``alignment_score = matches - mismatches``

    This is a simple ungapped sequence-comparison metric, not a biological
    cleavage-risk prediction.

    Args:
        guide: User-provided 20-nt guide sequence (DNA or RNA).
        protospacer: Candidate 20-nt protospacer sequence (DNA or RNA).
        expected_length: Expected sequence length for both sequences (default 20).

    Returns:
        A tuple ``(mismatch_count, mismatch_positions, score)`` where:
            - ``mismatch_count`` is the total number of differing positions.
            - ``mismatch_positions`` is a tuple of 0-based differing indices in ascending order.
            - ``score`` is ``float(matches - mismatches)``.

    Raises:
        TypeError: If ``guide`` or ``protospacer`` is not a string.
        ValueError: If either sequence is invalid or has length != ``expected_length``.
    """
    norm_guide = validate_guide_sequence(guide, expected_length=expected_length)
    norm_protospacer = _validate_candidate_protospacer(
        protospacer, expected_length=expected_length
    )

    mismatch_positions = tuple(
        idx
        for idx, (g_base, p_base) in enumerate(zip(norm_guide, norm_protospacer))
        if g_base != p_base
    )
    mismatch_count = len(mismatch_positions)
    matches = expected_length - mismatch_count
    score = float(matches - mismatch_count)

    return mismatch_count, mismatch_positions, score


def align_candidate_cpu(
    guide: str,
    candidate: CandidateResult,
    expected_length: int = DEFAULT_GUIDE_LENGTH,
) -> CandidateResult:
    """Align a single CandidateResult against a guide sequence on CPU.

    Returns a new ``CandidateResult`` with updated ``mismatch_count``,
    ``mismatch_positions``, and ``score`` while preserving all original
    metadata without mutating the input object.

    Args:
        guide: User-provided 20-nt guide sequence.
        candidate: ``CandidateResult`` instance produced by the PAM scanner.
        expected_length: Expected protospacer/guide length (default 20).

    Returns:
        A new ``CandidateResult`` with populated alignment metrics.

    Raises:
        TypeError: If ``candidate`` is not a ``CandidateResult`` or ``guide`` is not a string.
        ValueError: If ``guide`` or ``candidate.protospacer`` fails validation.
    """
    if not isinstance(candidate, CandidateResult):
        raise TypeError(
            f"candidate must be a CandidateResult instance, got {type(candidate).__name__}"
        )

    mismatch_count, mismatch_positions, score = align_protospacer_cpu(
        guide, candidate.protospacer, expected_length=expected_length
    )

    return CandidateResult(
        sequence_id=candidate.sequence_id,
        coordinates=candidate.coordinates,
        strand=candidate.strand,
        protospacer=candidate.protospacer,
        pam=candidate.pam,
        mismatch_count=mismatch_count,
        mismatch_positions=mismatch_positions,
        score=score,
    )


def align_candidates_cpu(
    guide: str,
    candidates: Sequence[CandidateResult],
    expected_length: int = DEFAULT_GUIDE_LENGTH,
) -> List[CandidateResult]:
    """Align a sequence of CandidateResult objects against a guide on CPU.

    Validates the guide sequence up front (even if ``candidates`` is empty)
    and returns a new list of updated ``CandidateResult`` objects without
    mutating the inputs.

    Args:
        guide: User-provided 20-nt guide sequence.
        candidates: Sequence of ``CandidateResult`` objects to align.
        expected_length: Expected protospacer/guide length (default 20).

    Returns:
        List of new ``CandidateResult`` instances with updated mismatch metrics and scores.
    """
    norm_guide = validate_guide_sequence(guide, expected_length=expected_length)
    return [
        align_candidate_cpu(norm_guide, candidate, expected_length=expected_length)
        for candidate in candidates
    ]


def search_and_align_cpu(
    guide: str,
    sequence: str,
    sequence_id: str = "seq",
    strands: str = "+-",
    sort_by_coordinate: bool = False,
) -> List[CandidateResult]:
    """Run end-to-end CPU PAM scanning and guide alignment on a single DNA sequence.

    Args:
        guide: 20-nt guide RNA or DNA sequence.
        sequence: Reference DNA sequence to scan.
        sequence_id: Identifier for the reference sequence.
        strands: Strands to scan ('+', '-', or '+-').
        sort_by_coordinate: Whether to sort results by reference start coordinate.

    Returns:
        List of aligned ``CandidateResult`` objects.
    """
    norm_guide = validate_guide_sequence(guide, expected_length=DEFAULT_GUIDE_LENGTH)
    candidates = find_pam_candidates(
        sequence=sequence,
        sequence_id=sequence_id,
        protospacer_len=DEFAULT_GUIDE_LENGTH,
        strands=strands,
        sort_by_coordinate=sort_by_coordinate,
    )
    return align_candidates_cpu(norm_guide, candidates, expected_length=DEFAULT_GUIDE_LENGTH)


def align_fasta_cpu(
    guide: str,
    fasta_path: Union[str, Path],
    strands: str = "+-",
    sort_by_coordinate: bool = False,
) -> List[CandidateResult]:
    """Run end-to-end CPU PAM scanning and guide alignment on a FASTA file.

    Args:
        guide: 20-nt guide RNA or DNA sequence.
        fasta_path: Path to the reference FASTA file.
        strands: Strands to scan ('+', '-', or '+-').
        sort_by_coordinate: Whether to sort candidates within each sequence by coordinate.

    Returns:
        List of aligned ``CandidateResult`` objects across all FASTA records.
    """
    norm_guide = validate_guide_sequence(guide, expected_length=DEFAULT_GUIDE_LENGTH)
    records = read_fasta(fasta_path)
    all_results: List[CandidateResult] = []
    for seq_id, seq in records:
        all_results.extend(
            search_and_align_cpu(
                guide=norm_guide,
                sequence=seq,
                sequence_id=seq_id,
                strands=strands,
                sort_by_coordinate=sort_by_coordinate,
            )
        )
    return all_results
