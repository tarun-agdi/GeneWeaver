"""
Unit tests for SpCas9 PAM scanning and candidate extraction (src/geneweaver/scanner.py).
"""

from pathlib import Path
import sys
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.models import CandidateResult
from geneweaver.scanner import find_pam_candidates, find_spcas9_candidates
from geneweaver.sequence import read_fasta, reverse_complement


class TestSpCas9PamDetection:
    """Tests for NGG PAM detection on forward and reverse strands."""

    def test_all_four_canonical_ngg_pams_forward(self):
        # 20-nt protospacer followed by AGG, CGG, GGG, TGG
        protospacer = "ACGTACGTACGTACGTACGT"
        for pam in ("AGG", "CGG", "GGG", "TGG"):
            seq = protospacer + pam
            results = find_pam_candidates(seq, sequence_id="test_chr")
            assert len(results) == 1
            assert results[0].strand == "+"
            assert results[0].protospacer == protospacer
            assert results[0].pam == pam
            assert results[0].coordinates == (0, 23)

    def test_non_ngg_patterns_rejected(self):
        protospacer = "ACGTACGTACGTACGTACGT"
        # None of these end in GG
        for non_pam in ("AGA", "AGC", "AGT", "TGA", "CGC", "AAG"):
            seq = protospacer + non_pam
            results = find_pam_candidates(seq, sequence_id="test_chr", strands="+")
            assert len(results) == 0

    def test_overlapping_pams_detected(self):
        # Protospacer + AGGG (contains AGG and GGG overlapping)
        protospacer = "ACGTACGTACGTACGTACGT"
        seq = protospacer + "AGGG"
        results = find_pam_candidates(seq, sequence_id="test_chr", strands="+")
        assert len(results) == 2
        # First PAM: AGG at (0, 23)
        assert results[0].pam == "AGG"
        assert results[0].coordinates == (0, 23)
        # Second PAM: GGG at (1, 24)
        assert results[1].pam == "GGG"
        assert results[1].coordinates == (1, 24)


class TestNoPamFound:
    """Tests for sequences lacking PAM sites."""

    def test_poly_at_sequence_no_pam(self):
        # AT repeats have no G or C, thus cannot form NGG or CCN on either strand
        seq = "AT" * 30
        results = find_pam_candidates(seq, sequence_id="test_no_pam")
        assert len(results) == 0

    def test_short_sequence_below_window_no_pam(self):
        # Shorter than 23 nt cannot host a complete 20-nt protospacer + 3-nt PAM
        results = find_pam_candidates("ACGTACGTACGTACGTACGG", sequence_id="short")
        assert len(results) == 0


class TestForwardStrandExtraction:
    """Tests for forward-strand extraction and coordinates."""

    def test_exact_23nt_window(self):
        protospacer = "A" * 20
        seq = protospacer + "TGG"
        results = find_pam_candidates(seq, sequence_id="fwd_exact", strands="+")
        assert len(results) == 1
        c = results[0]
        assert c.strand == "+"
        assert c.coordinates == (0, 23)
        assert c.start == 0
        assert c.end == 23
        assert c.protospacer == protospacer
        assert c.pam == "TGG"

    def test_embedded_forward_target(self):
        prefix = "TTTT"  # 4 nt
        protospacer = "ACGTACGTACGTACGTACGT"  # 20 nt
        pam = "CGG"  # 3 nt
        suffix = "AAAA"  # 4 nt
        seq = prefix + protospacer + pam + suffix  # 31 nt
        results = find_pam_candidates(seq, sequence_id="fwd_embedded", strands="+")
        assert len(results) == 1
        c = results[0]
        assert c.coordinates == (4, 27)
        assert c.start == 4
        assert c.end == 27
        assert c.protospacer == protospacer
        assert c.pam == pam


class TestReverseStrandExtractionAndCoordinates:
    """Tests for reverse-complement strand extraction and coordinate mapping."""

    def test_exact_23nt_reverse_target(self):
        # On reference strand: CCA (corresponds to TGG on reverse) + 20-nt protospacer complement
        # If reverse strand protospacer is "A"*20, reference sequence has "T"*20
        ref_seq = "CCA" + "T" * 20
        results = find_pam_candidates(ref_seq, sequence_id="rev_exact", strands="-")
        assert len(results) == 1
        c = results[0]
        assert c.strand == "-"
        assert c.coordinates == (0, 23)
        assert c.start == 0
        assert c.end == 23
        assert c.protospacer == "A" * 20
        assert c.pam == "TGG"
        # Mathematical invariance check: reverse_complement(ref_seq[start:end]) == protospacer + pam
        assert reverse_complement(ref_seq[c.start : c.end]) == c.protospacer + c.pam

    def test_embedded_reverse_target_coordinates(self):
        prefix = "GGGG"  # 4 nt
        pam_comp = "CCG"  # 3 nt (corresponds to reverse PAM "CGG")
        protospacer_comp = "TGATCGATCGATCGATCGAC"  # 20 nt
        suffix = "AAAA"  # 4 nt
        ref_seq = prefix + pam_comp + protospacer_comp + suffix  # 31 nt

        results = find_pam_candidates(ref_seq, sequence_id="rev_embedded", strands="-")
        assert len(results) == 1
        c = results[0]
        assert c.strand == "-"
        # Target locus on reference starts at prefix (4) and spans 23 nt (4 to 27)
        assert c.coordinates == (4, 27)
        assert c.start == 4
        assert c.end == 27
        assert c.pam == "CGG"
        assert c.protospacer == reverse_complement(protospacer_comp)
        # Verify slice invariance
        assert reverse_complement(ref_seq[c.start : c.end]) == c.protospacer + c.pam


class TestIncompleteWindows:
    """Tests that PAMs with incomplete protospacer windows are skipped."""

    def test_pam_at_start_of_forward_strand_skipped(self):
        # PAM "TGG" at index 0 followed by 20 nt; has no 20-nt protospacer before it
        seq = "TGG" + "A" * 20
        results = find_pam_candidates(seq, strands="+")
        assert len(results) == 0

    def test_pam_with_19nt_prefix_skipped(self):
        # Only 19 nt before PAM
        seq = "A" * 19 + "TGG"
        results = find_pam_candidates(seq, strands="+")
        assert len(results) == 0

    def test_reverse_pam_at_end_of_sequence_skipped(self):
        # CCN at 3' end of reference sequence corresponds to NGG at 5' end of reverse strand (incomplete)
        seq = "A" * 20 + "CCA"
        results = find_pam_candidates(seq, strands="-")
        assert len(results) == 0


class TestMultiplePamSitesAndStrands:
    """Tests for multi-site and dual-strand extraction."""

    def test_dual_strand_extraction(self):
        # Sequence with one forward site and one reverse site separated by a spacer
        fwd_protospacer = "ACGTACGTACGTACGTACGT"
        rev_protospacer_comp = "T" * 20
        seq = fwd_protospacer + "TGG" + "ATAT" + "CCA" + rev_protospacer_comp

        results = find_pam_candidates(seq, sequence_id="dual", strands="+-")
        assert len(results) == 2

        fwd = [r for r in results if r.strand == "+"]
        rev = [r for r in results if r.strand == "-"]
        assert len(fwd) == 1
        assert len(rev) == 1

        assert fwd[0].coordinates == (0, 23)
        assert fwd[0].pam == "TGG"
        assert fwd[0].protospacer == fwd_protospacer

        # rev site starts after fwd (23) + spacer (4) = 27
        assert rev[0].coordinates == (27, 50)
        assert rev[0].pam == "TGG"
        assert rev[0].protospacer == "A" * 20

    def test_strand_filtering_options(self):
        seq = "ACGTACGTACGTACGTACGTTGGATATCCATTTTTTTTTTTTTTTTTTTT"
        assert len(find_pam_candidates(seq, strands="+")) == 1
        assert len(find_pam_candidates(seq, strands="-")) == 1
        assert len(find_pam_candidates(seq, strands="+-")) == 2

    def test_sort_by_coordinate_flag(self):
        # Reverse site at (0, 23), forward site at (25, 48)
        seq = "CCA" + "T" * 20 + "AA" + "A" * 20 + "TGG"
        # Default: forward first, then reverse
        default_order = find_pam_candidates(seq, sort_by_coordinate=False)
        assert [c.strand for c in default_order] == ["+", "-"]

        # With sort_by_coordinate=True: (0, 23) '-' comes before (25, 48) '+'
        coord_order = find_pam_candidates(seq, sort_by_coordinate=True)
        assert [c.strand for c in coord_order] == ["-", "+"]
        assert coord_order[0].start < coord_order[1].start


class TestCandidateResultCompatibility:
    """Tests verifying compliance with CandidateResult contract."""

    def test_all_fields_populated_and_valid(self):
        seq = "ACGTACGTACGTACGTACGTTGG"
        results = find_pam_candidates(seq, sequence_id="chrX")
        assert len(results) == 1
        c = results[0]
        assert isinstance(c, CandidateResult)
        assert c.sequence_id == "chrX"
        assert c.coordinates == (0, 23)
        assert c.start == 0
        assert c.end == 23
        assert c.strand == "+"
        assert c.protospacer == "ACGTACGTACGTACGTACGT"
        assert c.pam == "TGG"
        assert c.mismatch_count == 0
        assert c.mismatch_positions == ()
        assert c.score == 0.0

        # to_dict verification
        d = c.to_dict()
        assert d["sequence_id"] == "chrX"
        assert d["coordinates"] == (0, 23)
        assert d["strand"] == "+"


class TestInputValidation:
    """Tests for defensive parameter checks."""

    def test_invalid_sequence_characters_raise_error(self):
        with pytest.raises(ValueError, match=r"Invalid character\(s\)"):
            find_pam_candidates("ACGTACGTACGTACGTACGTNAGG")

    def test_empty_sequence_raises_error(self):
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            find_pam_candidates("")

    def test_non_string_sequence_raises_error(self):
        with pytest.raises(TypeError, match="Sequence must be a string"):
            find_pam_candidates(12345)

    def test_empty_sequence_id_raises_error(self):
        with pytest.raises(ValueError, match="sequence_id cannot be empty"):
            find_pam_candidates("ACGTACGTACGTACGTACGTTGG", sequence_id="   ")

    def test_invalid_strands_specification_raises_error(self):
        with pytest.raises(ValueError, match="strands must contain"):
            find_pam_candidates("ACGTACGTACGTACGTACGTTGG", strands="forward")

    def test_invalid_protospacer_len_raises_error(self):
        with pytest.raises(ValueError, match="protospacer_len must be a positive integer"):
            find_pam_candidates("ACGTACGTACGTACGTACGTTGG", protospacer_len=0)


class TestSampleFastaIntegration:
    """Tests scanning candidates against data/sample.fasta."""

    def test_sample_fasta_detection(self):
        fasta_path = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"
        assert fasta_path.exists(), "data/sample.fasta should exist"

        records = read_fasta(fasta_path)
        records_dict = {rec[0].split()[0]: rec[1] for rec in records}

        # 1. seq1_forward_target -> 1 forward candidate
        cands1 = find_pam_candidates(records_dict["seq1_forward_target"], sequence_id="seq1")
        assert len(cands1) == 1
        assert cands1[0].strand == "+"
        assert cands1[0].pam == "AGG"

        # 2. seq2_reverse_target -> 1 reverse candidate
        cands2 = find_pam_candidates(records_dict["seq2_reverse_target"], sequence_id="seq2")
        assert len(cands2) == 1
        assert cands2[0].strand == "-"
        assert cands2[0].pam == "TGG"

        # 3. seq3_dual_target -> 1 forward + 1 reverse candidate
        cands3 = find_pam_candidates(records_dict["seq3_dual_target"], sequence_id="seq3")
        assert len(cands3) == 2
        strands = {c.strand for c in cands3}
        assert strands == {"+", "-"}

        # 4. seq4_no_pam -> 0 candidates
        cands4 = find_pam_candidates(records_dict["seq4_no_pam"], sequence_id="seq4")
        assert len(cands4) == 0

    def test_convenience_alias_equality(self):
        seq = "ACGTACGTACGTACGTACGTTGG"
        res1 = find_pam_candidates(seq)
        res2 = find_spcas9_candidates(seq)
        assert res1 == res2

