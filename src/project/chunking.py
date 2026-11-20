
from dataclasses import dataclass, replace

from .pam_scanner import Candidate, scan_pam


@dataclass
class SequenceChunk:
    sequence_id: str
    chunk_index: int
    start: int
    end: int
    sequence: str


def chunk_sequence(
    sequence_id: str,
    sequence: str,
    chunk_size: int = 1_000_000,
    overlap: int = 22
):
    """
    Split a DNA sequence into overlapping chunks.

    Coordinates use 0-based, end-exclusive indexing.
    """

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError("Invalid overlap")

    chunks = []

    start = 0
    chunk_index = 0
    sequence_length = len(sequence)

    while start < sequence_length:

        end = min(start + chunk_size, sequence_length)

        chunk = SequenceChunk(
            sequence_id=sequence_id,
            chunk_index=chunk_index,
            start=start,
            end=end,
            sequence=sequence[start:end]
        )

        chunks.append(chunk)

        if end == sequence_length:
            break

        start += chunk_size - overlap
        chunk_index += 1

    return chunks


def scan_chunked_sequence(
    sequence_id: str,
    sequence: str,
    chunk_size: int = 1_000_000,
    overlap: int = 22
) -> list[Candidate]:
    """
    Scan overlapping chunks and return globally
    positioned, deduplicated PAM candidates.
    """

    chunks = chunk_sequence(
        sequence_id,
        sequence,
        chunk_size,
        overlap
    )

    all_candidates = []
    seen = set()

    for chunk in chunks:

        candidates = scan_pam(
            chunk.sequence_id,
            chunk.sequence
        )

        for candidate in candidates:

            global_candidate = replace(
                candidate,
                start=candidate.start + chunk.start,
                end=candidate.end + chunk.start,
                pam_start=candidate.pam_start + chunk.start,
                pam_end=candidate.pam_end + chunk.start
            )

            key = (
                global_candidate.start,
                global_candidate.end,
                global_candidate.strand,
                global_candidate.pam_start,
                global_candidate.pam_end
            )

            if key not in seen:
                seen.add(key)
                all_candidates.append(global_candidate)

    return sorted(
        all_candidates,
        key=lambda item: (item.start, item.strand)
    )
