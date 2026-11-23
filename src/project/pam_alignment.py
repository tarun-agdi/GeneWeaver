
from dataclasses import dataclass
from typing import Iterable


@dataclass
class OffTarget:
    sequence_id: str
    position: int
    strand: str
    protospacer: str
    pam: str
    mismatches: int
    score: float


def count_mismatches(guide: str, target: str) -> int:
    """
    Count mismatches between equal-length DNA sequences.
    """
    guide = guide.upper()
    target = target.upper()

    if len(guide) != len(target):
        raise ValueError("Guide and target must have equal lengths.")

    return sum(
        base1 != base2
        for base1, base2 in zip(guide, target)
    )


def align_pam_candidates(
    guide: str,
    candidates: Iterable,
    max_mismatches: int = 3
) -> list[OffTarget]:

    guide = guide.upper()

    if not guide or any(base not in "ATCG" for base in guide):
        raise ValueError("Guide must contain valid DNA bases.")

    if max_mismatches < 0:
        raise ValueError("max_mismatches cannot be negative.")

    results = []

    for candidate in candidates:

        # Candidate attributes produced by the PAM scanner
        target = candidate.protospacer.upper()

        mismatches = count_mismatches(guide, target)

        if mismatches <= max_mismatches:

            score = (
                1 - mismatches / len(guide)
            ) * 100

            results.append(
                OffTarget(
                    sequence_id=candidate.sequence_id,
                    position=candidate.position,
                    strand=candidate.strand,
                    protospacer=target,
                    pam=candidate.pam,
                    mismatches=mismatches,
                    score=round(score, 2)
                )
            )

    return sorted(
        results,
        key=lambda result: (
            result.mismatches,
            result.position
        )
    )
