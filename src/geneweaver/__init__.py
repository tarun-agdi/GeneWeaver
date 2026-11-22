"""
GeneWeaver
==========

GPU-Accelerated CRISPR Alignment Engine.

GeneWeaver aims to provide an efficient and scalable framework
for identifying and aligning CRISPR guide RNA sequences against
reference DNA sequences.

Project scope:
    - CRISPR guide sequence analysis
    - PAM identification and strand-aware scanning
    - Sequence alignment and mismatch analysis
    - CPU-based alignment support
    - GPU-accelerated alignment support
    - Structured alignment results

The package is designed with modularity and extensibility
in mind, allowing CPU and GPU execution backends.
"""

from geneweaver.models import CandidateResult
from geneweaver.scanner import (
    detect_pam_candidates,
    extract_forward_candidates,
    extract_forward_protospacers,
    find_forward_pams,
    find_pam_candidates,
    find_spcas9_candidates,
    scan_forward_pams,
    scan_pam_candidates,
)
from geneweaver.sequence import (
    VALID_DNA_BASES,
    read_fasta,
    reverse_complement,
    validate_dna_sequence,
)

__version__ = "0.1.0"

__all__ = [
    "CandidateResult",
    "VALID_DNA_BASES",
    "validate_dna_sequence",
    "reverse_complement",
    "read_fasta",
    "find_forward_pams",
    "scan_forward_pams",
    "extract_forward_protospacers",
    "extract_forward_candidates",
    "find_pam_candidates",
    "find_spcas9_candidates",
    "detect_pam_candidates",
    "scan_pam_candidates",
]