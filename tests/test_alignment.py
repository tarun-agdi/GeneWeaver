"""
Unit and integration tests for CPU reference alignment & mismatch scoring (alignment.py).
"""

from pathlib import Path
import sys
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.alignment import (
    align_candidate_cpu,
    align_candidates_cpu,
    align_fasta_cpu,
    align_protospacer_cpu,
    search_and_align_cpu,
    validate_guide_sequence,
)
from geneweaver.models import CandidateResult


class TestAlignProtospacerCpu:
    """Tests for fixed-length ungapped CPU protospacer alignment."""

    GUIDE = "ACGTACGTACGTACGTACGT"  # 20 nt

    def test_perfect_match(self):
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, self.GUIDE
        )
        assert mismatch_count == 0
        assert mismatch_positions == ()
        assert score == 20.0

    def test_one_nucleotide_mismatch(self):
        # Prompt example: mismatch at index 4 (A -> T)
        candidate = "ACGTTCGTACGTACGTACGT"
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, candidate
        )
        assert mismatch_count == 1
        assert mismatch_positions == (4,)
        assert score == 18.0  # 19 matches - 1 mismatch

    def test_multiple_mismatches_ascending_positions(self):
        # Mismatches at indices 1 (C->A), 7 (T->A), 15 (T->C) -> 3 mismatches
        candidate = "AAGTACGAACGTACGCACGT"
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, candidate
        )
        assert mismatch_count == 3
        assert mismatch_positions == (1, 7, 15)
        assert list(mismatch_positions) == sorted(mismatch_positions)
        assert score == 14.0  # 17 matches - 3 mismatches

    def test_mismatch_at_first_position(self):
        candidate = "TCGTACGTACGTACGTACGT"  # index 0 differs (A -> T)
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, candidate
        )
        assert mismatch_count == 1
        assert mismatch_positions == (0,)
        assert score == 18.0

    def test_mismatch_at_last_position(self):
        candidate = "ACGTACGTACGTACGTACGA"  # index 19 differs (T -> A)
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, candidate
        )
        assert mismatch_count == 1
        assert mismatch_positions == (19,)
        assert score == 18.0

    def test_all_twenty_positions_mismatching(self):
        # Every base differs from ACGTACGTACGTACGTACGT -> TGCATGCATGCATGCATGCA
        candidate = "TGCATGCATGCATGCATGCA"
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            self.GUIDE, candidate
        )
        assert mismatch_count == 20
        assert mismatch_positions == tuple(range(20))
        assert score == -20.0  # 0 matches - 20 mismatches

    def test_rna_u_to_dna_t_and_lowercase_normalization(self):
        rna_guide = "acguacguacguacguacgu"
        dna_candidate = "ACGTACGTACGTACGTACGT"
        mismatch_count, mismatch_positions, score = align_protospacer_cpu(
            rna_guide, dna_candidate
        )
        assert mismatch_count == 0
        assert mismatch_positions == ()
        assert score == 20.0

        # Also RNA U in candidate
        mismatch_count2, mismatch_positions2, score2 = align_protospacer_cpu(
            self.GUIDE, "acguucguacguacguacgu"
        )
        assert mismatch_count2 == 1
        assert mismatch_positions2 == (4,)
        assert score2 == 18.0


class TestValidationErrors:
    """Tests for invalid guide and candidate inputs."""

    VALID_20 = "ACGTACGTACGTACGTACGT"

    def test_invalid_guide_sequences(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            validate_guide_sequence("")

        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            validate_guide_sequence("ACGTACGTACGTACGTACGN")

        with pytest.raises(TypeError, match="Sequence must be a string"):
            validate_guide_sequence(None)

    def test_guide_shorter_or_longer_than_20_nt(self):
        with pytest.raises(ValueError, match="exactly 20 nucleotides, got 19"):
            validate_guide_sequence("A" * 19)

        with pytest.raises(ValueError, match="exactly 20 nucleotides, got 21"):
            validate_guide_sequence("A" * 21)

    def test_invalid_candidate_sequences(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            align_protospacer_cpu(self.VALID_20, "")

        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            align_protospacer_cpu(self.VALID_20, "ACGTACGTACGTACGTACG-")

        with pytest.raises(TypeError, match="Sequence must be a string"):
            align_protospacer_cpu(self.VALID_20, 12345)

    def test_candidate_shorter_or_longer_than_20_nt(self):
        with pytest.raises(ValueError, match="Candidate protospacer must be exactly 20 nucleotides, got 18"):
            align_protospacer_cpu(self.VALID_20, "A" * 18)

        with pytest.raises(ValueError, match="Candidate protospacer must be exactly 20 nucleotides, got 23"):
            align_protospacer_cpu(self.VALID_20, "A" * 23)

    def test_non_candidate_result_rejected(self):
        with pytest.raises(TypeError, match="CandidateResult instance"):
            align_candidate_cpu(self.VALID_20, {"protospacer": self.VALID_20})  # type: ignore[arg-type]

    def test_invalid_expected_length_rejected(self):
        for bad_len in (0, -5, True, "20"):
            with pytest.raises(ValueError, match="expected_length must be a positive integer"):
                validate_guide_sequence(self.VALID_20, expected_length=bad_len)  # type: ignore[arg-type]

    def test_invalid_candidates_iterable_rejected(self):
        with pytest.raises(TypeError, match="candidates must be an iterable"):
            align_candidates_cpu(self.VALID_20, None)  # type: ignore[arg-type]

        with pytest.raises(TypeError, match="candidates must be an iterable"):
            align_candidates_cpu(self.VALID_20, self.VALID_20)  # type: ignore[arg-type]


class TestCandidateResultIntegration:
    """Tests for CandidateResult alignment, metadata preservation, and immutability."""

    GUIDE = "ACGTACGTACGTACGTACGT"

    def test_metadata_preservation_and_non_mutation(self):
        original = CandidateResult(
            sequence_id="chr7_locus",
            coordinates=(120, 143),
            strand="-",
            protospacer="ACGTTCGTACGTACGTACGA",  # mismatches at 4 and 19
            pam="CGG",
            mismatch_count=0,
            mismatch_positions=(),
            score=0.0,
        )

        updated = align_candidate_cpu(self.GUIDE, original)

        # Original object must remain untouched
        assert original.mismatch_count == 0
        assert original.mismatch_positions == ()
        assert original.score == 0.0

        # Updated object preserves all metadata and populates alignment metrics
        assert updated is not original
        assert updated.sequence_id == "chr7_locus"
        assert updated.coordinates == (120, 143)
        assert updated.start == 120
        assert updated.end == 143
        assert updated.strand == "-"
        assert updated.protospacer == "ACGTTCGTACGTACGTACGA"
        assert updated.pam == "CGG"
        assert updated.mismatch_count == 2
        assert updated.mismatch_positions == (4, 19)
        assert updated.score == 16.0  # 18 matches - 2 mismatches

        # Serialization still works cleanly
        d = updated.to_dict()
        assert d["mismatch_count"] == 2
        assert d["mismatch_positions"] == [4, 19]
        assert d["score"] == 16.0

    def test_multiple_candidates_different_mismatch_counts(self):
        c_perfect = CandidateResult(
            sequence_id="s1", coordinates=(0, 23), strand="+",
            protospacer=self.GUIDE, pam="AGG", mismatch_count=0,
        )
        c_one = CandidateResult(
            sequence_id="s1", coordinates=(30, 53), strand="+",
            protospacer="TCGTACGTACGTACGTACGT", pam="TGG", mismatch_count=0,
        )
        c_all = CandidateResult(
            sequence_id="s1", coordinates=(60, 83), strand="-",
            protospacer="TGCATGCATGCATGCATGCA", pam="GGG", mismatch_count=0,
        )

        aligned = align_candidates_cpu(self.GUIDE, [c_perfect, c_one, c_all])
        assert len(aligned) == 3

        assert aligned[0].mismatch_count == 0
        assert aligned[0].mismatch_positions == ()
        assert aligned[0].score == 20.0

        assert aligned[1].mismatch_count == 1
        assert aligned[1].mismatch_positions == (0,)
        assert aligned[1].score == 18.0

        assert aligned[2].mismatch_count == 20
        assert aligned[2].mismatch_positions == tuple(range(20))
        assert aligned[2].score == -20.0

    def test_empty_candidate_list_still_validates_guide(self):
        assert align_candidates_cpu(self.GUIDE, []) == []
        with pytest.raises(ValueError, match="exactly 20 nucleotides"):
            align_candidates_cpu("ACGT", [])

    def test_realigning_already_aligned_candidate_replaces_metrics(self):
        already_aligned = CandidateResult(
            sequence_id="s1",
            coordinates=(0, 23),
            strand="+",
            protospacer=self.GUIDE,
            pam="AGG",
            mismatch_count=2,
            mismatch_positions=(0, 1),
            score=16.0,
        )
        realigned = align_candidate_cpu(self.GUIDE, already_aligned)
        assert realigned.mismatch_count == 0
        assert realigned.mismatch_positions == ()
        assert realigned.score == 20.0
        # Original stays untouched
        assert already_aligned.mismatch_count == 2
        assert already_aligned.score == 16.0


class TestEndToEndCpuPipeline:
    """Tests integrating PAM scanning, FASTA reading, and CPU alignment."""

    def test_search_and_align_cpu_synthetic(self):
        guide = "ACGUACGUACGUACGUACGU"  # RNA form of ACGTACGTACGTACGTACGT
        # Forward perfect match at (0, 23) + spacer + reverse 1-mismatch target
        # Reverse protospacer "TCGTACGTACGTACGTACGT" -> rc is "ACGTACGTACGTACGTACGA"
        seq = (
            "ACGTACGTACGTACGTACGT" + "TGG"  # fwd perfect match
            + "ATAT"
            + "CCA" + "ACGTACGTACGTACGTACGA"  # rev 1-mismatch at index 0
        )
        results = search_and_align_cpu(guide=guide, sequence=seq, sequence_id="syn_1")
        assert len(results) == 2

        fwd = results[0]
        assert fwd.strand == "+"
        assert fwd.coordinates == (0, 23)
        assert fwd.mismatch_count == 0
        assert fwd.mismatch_positions == ()
        assert fwd.score == 20.0

        rev = results[1]
        assert rev.strand == "-"
        assert rev.coordinates == (27, 50)
        assert rev.mismatch_count == 1
        assert rev.mismatch_positions == (0,)
        assert rev.score == 18.0

    def test_align_fasta_cpu_with_sample_fasta(self):
        fasta_path = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"
        # Guide matching seq1_forward_target and seq2_reverse_target perfectly:
        guide = "GTCGATCGATCGATCGATCA"
        results = align_fasta_cpu(guide=guide, fasta_path=fasta_path)

        # sample.fasta has 4 total PAM candidates across seq1 (1), seq2 (1), seq3 (2), seq4 (0)
        assert len(results) == 4

        by_seq = {r.sequence_id.split()[0]: [] for r in results}
        for r in results:
            by_seq[r.sequence_id.split()[0]].append(r)

        # seq1_forward_target: perfect match (0 mismatches, score 20.0)
        seq1_res = by_seq["seq1_forward_target"][0]
        assert seq1_res.strand == "+"
        assert seq1_res.mismatch_count == 0
        assert seq1_res.score == 20.0

        # seq2_reverse_target: perfect match (0 mismatches, score 20.0)
        seq2_res = by_seq["seq2_reverse_target"][0]
        assert seq2_res.strand == "-"
        assert seq2_res.mismatch_count == 0
        assert seq2_res.score == 20.0

        # seq3_dual_target: both candidates have mismatches against this guide
        assert all(r.mismatch_count > 0 for r in by_seq["seq3_dual_target"])
