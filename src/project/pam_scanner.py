
from dataclasses import dataclass


@dataclass
class Candidate:
    sequence_id: str
    start: int
    end: int
    strand: str
    protospacer: str
    pam: str
    pam_start: int
    pam_end: int


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans("ACGT", "TGCA")
    return sequence.translate(translation)[::-1]


def scan_pam(sequence_id: str, dna: str) -> list[Candidate]:

    dna = dna.upper()

    if not dna or any(base not in "ACGT" for base in dna):
        raise ValueError("DNA must contain only A, C, G and T")

    candidates = []
    n = len(dna)

    # Scan forward and reverse-complement strands
    for strand, seq in [
        ("+", dna),
        ("-", reverse_complement(dna))
    ]:

        # NGG PAM scanning
        for i in range(20, n - 2):

            pam = seq[i:i + 3]

            if len(pam) != 3 or pam[1:] != "GG":
                continue

            protospacer = seq[i - 20:i]

            # Convert coordinates to original reference
            if strand == "+":
                start = i - 20
                end = i
                pam_start = i
                pam_end = i + 3

            else:
                start = n - i
                end = n - i + 20
                pam_start = n - (i + 3)
                pam_end = n - i

            candidates.append(
                Candidate(
                    sequence_id=sequence_id,
                    start=start,
                    end=end,
                    strand=strand,
                    protospacer=protospacer,
                    pam=pam,
                    pam_start=pam_start,
                    pam_end=pam_end
                )
            )

    return candidates
