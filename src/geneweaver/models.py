"""
Data models for the GeneWeaver alignment engine.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Sequence, Tuple


@dataclass
class CandidateResult:
    """Represents a single CRISPR target candidate alignment result.

    Attributes:
        sequence_id: Identifier of the reference sequence (e.g. FASTA header).
        coordinates: 0-based, half-open reference coordinates (start, end).
        strand: Target strand orientation, either '+' (forward) or '-' (reverse).
        protospacer: The candidate protospacer sequence (typically 20-nt DNA).
        pam: The Protospacer Adjacent Motif sequence (e.g. 'NGG' / 'AGG').
        mismatch_count: Total number of mismatched positions against the guide.
        mismatch_positions: 0-based index positions of mismatches within the protospacer.
        score: Alignment score or ranking metric.
    """

    sequence_id: str
    coordinates: Tuple[int, int]
    strand: str
    protospacer: str
    pam: str
    mismatch_count: int
    mismatch_positions: Tuple[int, ...] = field(default_factory=tuple)
    score: float = 0.0

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        # 1. sequence_id validation
        if not isinstance(self.sequence_id, str):
            raise TypeError(f"sequence_id must be a string, got {type(self.sequence_id).__name__}")
        if not self.sequence_id.strip():
            raise ValueError("sequence_id cannot be empty or whitespace only")

        # 2. coordinates validation (0-based, half-open: start <= end, start >= 0)
        if not isinstance(self.coordinates, (tuple, list)):
            raise TypeError(f"coordinates must be a tuple or list of (start, end), got {type(self.coordinates).__name__}")
        if len(self.coordinates) != 2:
            raise ValueError(f"coordinates must contain exactly 2 integers (start, end), got {len(self.coordinates)}")

        start, end = self.coordinates
        if type(start) is not int or type(end) is not int:
            raise TypeError("coordinates start and end must both be integers")
        if start < 0:
            raise ValueError(f"coordinates start must be >= 0 (0-based), got {start}")
        if end < start:
            raise ValueError(f"coordinates end ({end}) must be >= start ({start}) for half-open interval")

        # Normalize coordinates to an immutable tuple of ints
        self.coordinates = (start, end)

        # 3. strand validation
        if not isinstance(self.strand, str):
            raise TypeError(f"strand must be a string, got {type(self.strand).__name__}")
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")

        # 4. protospacer validation
        if not isinstance(self.protospacer, str):
            raise TypeError(f"protospacer must be a string, got {type(self.protospacer).__name__}")
        if not self.protospacer:
            raise ValueError("protospacer cannot be empty")

        # 5. pam validation
        if not isinstance(self.pam, str):
            raise TypeError(f"pam must be a string, got {type(self.pam).__name__}")
        if not self.pam:
            raise ValueError("pam cannot be empty")

        # 6. mismatch_count validation
        if type(self.mismatch_count) is not int:
            raise TypeError(f"mismatch_count must be an integer, got {type(self.mismatch_count).__name__}")
        if self.mismatch_count < 0:
            raise ValueError(f"mismatch_count must be non-negative, got {self.mismatch_count}")
        if self.mismatch_count > len(self.protospacer):
            raise ValueError(
                f"mismatch_count ({self.mismatch_count}) cannot exceed protospacer length ({len(self.protospacer)})"
            )

        # 7. mismatch_positions validation
        if not isinstance(self.mismatch_positions, (tuple, list)):
            raise TypeError(
                f"mismatch_positions must be a tuple or list of integers, got {type(self.mismatch_positions).__name__}"
            )

        normalized_positions = []
        for pos in self.mismatch_positions:
            if type(pos) is not int:
                raise TypeError(f"Each mismatch position must be an integer, got {type(pos).__name__}")
            if pos < 0 or pos >= len(self.protospacer):
                raise ValueError(
                    f"mismatch position {pos} out of bounds for protospacer length {len(self.protospacer)} (0 to {len(self.protospacer) - 1})"
                )
            normalized_positions.append(pos)

        if len(normalized_positions) != self.mismatch_count:
            raise ValueError(
                f"Number of mismatch positions ({len(normalized_positions)}) does not match mismatch_count ({self.mismatch_count})"
            )

        # Normalize mismatch_positions to an immutable tuple
        self.mismatch_positions = tuple(normalized_positions)

        # 8. score validation
        if not isinstance(self.score, (int, float)) or isinstance(self.score, bool):
            raise TypeError(f"score must be a numeric value (int or float), got {type(self.score).__name__}")
        self.score = float(self.score)

    @property
    def start(self) -> int:
        """0-based start coordinate in the reference sequence (inclusive)."""
        return self.coordinates[0]

    @property
    def end(self) -> int:
        """0-based end coordinate in the reference sequence (exclusive)."""
        return self.coordinates[1]

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to a dictionary representation."""
        return {
            "sequence_id": self.sequence_id,
            "coordinates": self.coordinates,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "protospacer": self.protospacer,
            "pam": self.pam,
            "mismatch_count": self.mismatch_count,
            "mismatch_positions": list(self.mismatch_positions),
            "score": self.score,
        }

