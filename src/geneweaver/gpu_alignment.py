"""
GPU-accelerated (PyTorch) batched mismatch scoring prototype for GeneWeaver.

Implements fixed-length, ungapped nucleotide comparison between a 20-nt CRISPR
guide sequence and batches of candidate 20-nt protospacers using PyTorch tensor
operations (selected fallback over GASAL2-Py for portability).

Note:
    This module performs fixed-length ungapped mismatch scoring
    (``matches - mismatches``), NOT full gapped biological alignment or
    cleavage-risk prediction.

    When CUDA is unavailable (or when PyTorch is not installed locally and
    ``allow_cpu_fallback=True``), execution falls back to CPU tensor/NumPy
    matrix operations with the exact same input/output contract.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geneweaver.alignment import (
    DEFAULT_GUIDE_LENGTH,
    _validate_candidate_protospacer,
    validate_guide_sequence,
)
from geneweaver.models import CandidateResult
from geneweaver.scanner import find_pam_candidates

try:
    import torch

    _TORCH_AVAILABLE = True
except ImportError:  # pragma: no cover - depends on environment
    torch = None  # type: ignore[assignment]
    _TORCH_AVAILABLE = False

BASE_TO_INT: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}

__all__ = [
    "BASE_TO_INT",
    "is_torch_available",
    "is_cuda_available",
    "resolve_gpu_device",
    "align_candidates_gpu",
    "search_and_align_gpu",
]


def is_torch_available() -> bool:
    """Return True if PyTorch can be imported in the current environment."""
    return _TORCH_AVAILABLE


def is_cuda_available() -> bool:
    """Return True if PyTorch is installed and a CUDA GPU device is available."""
    return bool(_TORCH_AVAILABLE and torch.cuda.is_available())


def resolve_gpu_device(
    device: Optional[str] = None, allow_cpu_fallback: bool = True
) -> str:
    """Resolve the execution device for batched tensor scoring.

    Args:
        device: Optional explicit device ('cuda' or 'cpu'). If None, selects
            'cuda' when available, otherwise falls back to 'cpu' if allowed.
        allow_cpu_fallback: Whether to allow CPU fallback when CUDA is unavailable.

    Returns:
        Resolved device string ('cuda' or 'cpu').

    Raises:
        ValueError: If ``device`` is not None, 'cuda', or 'cpu'.
        RuntimeError: If CUDA is requested (or fallback is disabled) when CUDA
            or PyTorch is unavailable.
    """
    if device is not None and device not in ("cuda", "cpu"):
        raise ValueError(f"device must be 'cuda', 'cpu', or None, got {device!r}")

    if device == "cuda":
        if not is_cuda_available():
            raise RuntimeError(
                "CUDA device requested, but PyTorch CUDA is unavailable in this environment."
            )
        return "cuda"

    if device == "cpu":
        if not allow_cpu_fallback:
            raise RuntimeError("CPU fallback is disabled, but device='cpu' was requested.")
        return "cpu"

    # device is None: prefer CUDA when available
    if is_cuda_available():
        return "cuda"

    if allow_cpu_fallback:
        return "cpu"

    raise RuntimeError(
        "CUDA GPU is unavailable and allow_cpu_fallback=False."
    )


def _encode_sequence(seq: str) -> List[int]:
    """Encode a validated uppercase DNA sequence into integer base IDs."""
    return [BASE_TO_INT[base] for base in seq]


def _score_batch_torch(
    guide_encoded: List[int],
    batch_encoded: List[List[int]],
    expected_length: int,
    device: str,
) -> List[Tuple[int, Tuple[int, ...], float]]:
    """Score a batch of encoded protospacers using PyTorch tensor operations."""
    with torch.no_grad():
        guide_tensor = torch.tensor(guide_encoded, dtype=torch.uint8, device=device)
        cand_tensor = torch.tensor(batch_encoded, dtype=torch.uint8, device=device)

        # (batch_size, expected_length) boolean mask of differing positions
        mismatch_mask = cand_tensor != guide_tensor
        mismatch_counts = mismatch_mask.sum(dim=1)
        matches = expected_length - mismatch_counts
        scores = matches.to(torch.int64) - mismatch_counts.to(torch.int64)

        # Move compact results back to host once per batch
        mask_rows = mismatch_mask.cpu().tolist()
        counts_list = mismatch_counts.cpu().tolist()
        scores_list = scores.cpu().tolist()

    results: List[Tuple[int, Tuple[int, ...], float]] = []
    for mask_row, count, score in zip(mask_rows, counts_list, scores_list):
        positions = tuple(idx for idx, is_diff in enumerate(mask_row) if is_diff)
        results.append((int(count), positions, float(score)))
    return results


def _score_batch_numpy_fallback(
    guide_encoded: List[int],
    batch_encoded: List[List[int]],
    expected_length: int,
) -> List[Tuple[int, Tuple[int, ...], float]]:
    """Vectorized NumPy CPU fallback matching the PyTorch tensor contract."""
    guide_arr = np.asarray(guide_encoded, dtype=np.uint8)
    cand_arr = np.asarray(batch_encoded, dtype=np.uint8)

    mismatch_mask = cand_arr != guide_arr
    mismatch_counts = mismatch_mask.sum(axis=1)
    matches = expected_length - mismatch_counts
    scores = matches.astype(np.int64) - mismatch_counts.astype(np.int64)

    results: List[Tuple[int, Tuple[int, ...], float]] = []
    for mask_row, count, score in zip(
        mismatch_mask.tolist(), mismatch_counts.tolist(), scores.tolist()
    ):
        positions = tuple(idx for idx, is_diff in enumerate(mask_row) if is_diff)
        results.append((int(count), positions, float(score)))
    return results


def align_candidates_gpu(
    guide: str,
    candidates: Sequence[CandidateResult],
    expected_length: int = DEFAULT_GUIDE_LENGTH,
    batch_size: int = 1024,
    device: Optional[str] = None,
    allow_cpu_fallback: bool = True,
) -> List[CandidateResult]:
    """Align CandidateResult objects against a guide using batched tensor scoring.

    Uses PyTorch tensor operations on CUDA when available (or CPU fallback when
    permitted). Preserves candidate ordering and all original metadata without
    mutating input CandidateResult objects.

    Note:
        This performs fixed-length ungapped mismatch scoring, not gapped
        biological alignment or cleavage-risk prediction.

    Args:
        guide: User-provided 20-nt guide sequence (DNA or RNA).
        candidates: Iterable of ``CandidateResult`` objects to score.
        expected_length: Expected guide/protospacer length (default 20).
        batch_size: Number of candidates per tensor batch (default 1024).
        device: Explicit execution device ('cuda', 'cpu', or None for auto).
        allow_cpu_fallback: Whether to fall back to CPU if CUDA/PyTorch is unavailable.

    Returns:
        List of new ``CandidateResult`` objects with populated ``mismatch_count``,
        ``mismatch_positions``, and ``score``.

    Raises:
        TypeError: If ``guide`` or ``candidates`` has an invalid type.
        ValueError: If ``guide``, any candidate protospacer, ``batch_size``, or
            ``device`` is invalid.
        RuntimeError: If CUDA is required but unavailable.
    """
    norm_guide = validate_guide_sequence(guide, expected_length=expected_length)

    if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size <= 0:
        raise ValueError(f"batch_size must be a positive integer, got {batch_size!r}")

    if isinstance(candidates, (str, bytes)) or not isinstance(candidates, Iterable):
        raise TypeError(
            f"candidates must be an iterable of CandidateResult instances, got {type(candidates).__name__}"
        )

    candidate_list = list(candidates)
    resolved_device = resolve_gpu_device(
        device=device, allow_cpu_fallback=allow_cpu_fallback
    )

    if not candidate_list:
        return []

    encoded_candidates: List[List[int]] = []
    for candidate in candidate_list:
        if not isinstance(candidate, CandidateResult):
            raise TypeError(
                f"Each item in candidates must be a CandidateResult instance, got {type(candidate).__name__}"
            )
        norm_proto = _validate_candidate_protospacer(
            candidate.protospacer, expected_length=expected_length
        )
        encoded_candidates.append(_encode_sequence(norm_proto))

    guide_encoded = _encode_sequence(norm_guide)

    scored_metrics: List[Tuple[int, Tuple[int, ...], float]] = []
    for start_idx in range(0, len(encoded_candidates), batch_size):
        batch_encoded = encoded_candidates[start_idx : start_idx + batch_size]
        if is_torch_available():
            scored_metrics.extend(
                _score_batch_torch(
                    guide_encoded=guide_encoded,
                    batch_encoded=batch_encoded,
                    expected_length=expected_length,
                    device=resolved_device,
                )
            )
        else:
            scored_metrics.extend(
                _score_batch_numpy_fallback(
                    guide_encoded=guide_encoded,
                    batch_encoded=batch_encoded,
                    expected_length=expected_length,
                )
            )

    aligned_results: List[CandidateResult] = []
    for candidate, (mismatch_count, mismatch_positions, score) in zip(
        candidate_list, scored_metrics
    ):
        aligned_results.append(
            CandidateResult(
                sequence_id=candidate.sequence_id,
                coordinates=candidate.coordinates,
                strand=candidate.strand,
                protospacer=candidate.protospacer,
                pam=candidate.pam,
                mismatch_count=mismatch_count,
                mismatch_positions=mismatch_positions,
                score=score,
            )
        )

    return aligned_results


def search_and_align_gpu(
    guide: str,
    sequence: str,
    sequence_id: str = "seq",
    strands: str = "+-",
    sort_by_coordinate: bool = False,
    batch_size: int = 1024,
    device: Optional[str] = None,
    allow_cpu_fallback: bool = True,
) -> List[CandidateResult]:
    """Run end-to-end PAM scanning and batched GPU/fallback scoring on a sequence."""
    norm_guide = validate_guide_sequence(guide, expected_length=DEFAULT_GUIDE_LENGTH)
    candidates = find_pam_candidates(
        sequence=sequence,
        sequence_id=sequence_id,
        protospacer_len=DEFAULT_GUIDE_LENGTH,
        strands=strands,
        sort_by_coordinate=sort_by_coordinate,
    )
    return align_candidates_gpu(
        guide=norm_guide,
        candidates=candidates,
        expected_length=DEFAULT_GUIDE_LENGTH,
        batch_size=batch_size,
        device=device,
        allow_cpu_fallback=allow_cpu_fallback,
    )
