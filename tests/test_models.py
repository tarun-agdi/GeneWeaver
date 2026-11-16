"""
Unit tests for GeneWeaver data models (CandidateResult).
"""

import sys
from pathlib import Path
import pytest

# Ensure src/ is on sys.path
src_path = Path(__file__).resolve().parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from geneweaver.models import CandidateResult


class TestCandidateResultValid:
    """Tests for valid CandidateResult construction and behavior."""

    def test_valid_construction_perfect_match(self):
        result = CandidateResult(
            sequence_id="chr1",
            coordinates=(100, 123),
            strand="+",
            protospacer="ACGTACGTACGTACGTACGT",
            pam="TGG",
            mismatch_count=0,
            mismatch_positions=(),
            score=100.0,
        )
        assert result.sequence_id == "chr1"
        assert result.coordinates == (100, 123)
        assert result.start == 100
        assert result.end == 123
        assert result.strand == "+"
        assert result.protospacer == "ACGTACGTACGTACGTACGT"
        assert result.pam == "TGG"
        assert result.mismatch_count == 0
        assert result.mismatch_positions == ()
        assert result.score == 100.0

    def test_valid_construction_with_mismatches(self):
        result = CandidateResult(
            sequence_id="NC_000913.3",
            coordinates=(500, 523),
            strand="-",
            protospacer="ACGTACGTACGTACGTACGA",
            pam="CGG",
            mismatch_count=2,
            mismatch_positions=(0, 19),
            score=78.5,
        )
        assert result.strand == "-"
        assert result.mismatch_count == 2
        assert result.mismatch_positions == (0, 19)
        assert result.score == 78.5

    def test_list_normalization(self):
        """Verify lists for coordinates and mismatch_positions are converted to immutable tuples."""
        result = CandidateResult(
            sequence_id="seq_1",
            coordinates=[10, 33],
            strand="+",
            protospacer="AAAAAAAAAAAAAAAAAAAA",
            pam="AGG",
            mismatch_count=1,
            mismatch_positions=[5],
            score=90,
        )
        assert isinstance(result.coordinates, tuple)
        assert result.coordinates == (10, 33)
        assert isinstance(result.mismatch_positions, tuple)
        assert result.mismatch_positions == (5,)
        assert isinstance(result.score, float)
        assert result.score == 90.0

    def test_zero_length_interval_valid_half_open(self):
        """0-based half-open coordinates allow start == end for an empty slice boundary."""
        result = CandidateResult(
            sequence_id="seq_zero",
            coordinates=(50, 50),
            strand="+",
            protospacer="A",
            pam="NGG",
            mismatch_count=0,
            mismatch_positions=(),
            score=0.0,
        )
        assert result.start == 50
        assert result.end == 50

    def test_to_dict_serialization(self):
        result = CandidateResult(
            sequence_id="test_chr",
            coordinates=(0, 23),
            strand="+",
            protospacer="GGGGGGGGGGGGGGGGGGGG",
            pam="GGG",
            mismatch_count=0,
            mismatch_positions=(),
            score=1.0,
        )
        data = result.to_dict()
        assert data["sequence_id"] == "test_chr"
        assert data["coordinates"] == (0, 23)
        assert data["start"] == 0
        assert data["end"] == 23
        assert data["strand"] == "+"
        assert data["protospacer"] == "GGGGGGGGGGGGGGGGGGGG"
        assert data["pam"] == "GGG"
        assert data["mismatch_count"] == 0
        assert data["mismatch_positions"] == []
        assert data["score"] == 1.0


class TestCandidateResultInvalid:
    """Tests for invalid inputs and proper error handling."""

    def _base_kwargs(self):
        return {
            "sequence_id": "chr1",
            "coordinates": (10, 33),
            "strand": "+",
            "protospacer": "ACGTACGTACGTACGTACGT",
            "pam": "TGG",
            "mismatch_count": 0,
            "mismatch_positions": (),
            "score": 100.0,
        }

    # Sequence ID tests
    def test_invalid_sequence_id_type(self):
        kwargs = self._base_kwargs()
        kwargs["sequence_id"] = 123
        with pytest.raises(TypeError, match="sequence_id must be a string"):
            CandidateResult(**kwargs)

    def test_empty_sequence_id(self):
        kwargs = self._base_kwargs()
        kwargs["sequence_id"] = "   "
        with pytest.raises(ValueError, match="sequence_id cannot be empty"):
            CandidateResult(**kwargs)

    # Coordinates tests
    def test_coordinates_not_sequence(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = 100
        with pytest.raises(TypeError, match="coordinates must be a tuple or list"):
            CandidateResult(**kwargs)

    def test_coordinates_wrong_length(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = (10, 20, 30)
        with pytest.raises(ValueError, match="coordinates must contain exactly 2 integers"):
            CandidateResult(**kwargs)

    def test_coordinates_non_integer(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = (10.5, 33)
        with pytest.raises(TypeError, match="coordinates start and end must both be integers"):
            CandidateResult(**kwargs)

    def test_coordinates_boolean_rejected(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = (True, 33)
        with pytest.raises(TypeError, match="coordinates start and end must both be integers"):
            CandidateResult(**kwargs)

    def test_coordinates_negative_start(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = (-1, 22)
        with pytest.raises(ValueError, match="start must be >= 0"):
            CandidateResult(**kwargs)

    def test_coordinates_end_less_than_start(self):
        kwargs = self._base_kwargs()
        kwargs["coordinates"] = (50, 40)
        with pytest.raises(ValueError, match="end.*must be >= start"):
            CandidateResult(**kwargs)

    # Strand tests
    def test_invalid_strand_type(self):
        kwargs = self._base_kwargs()
        kwargs["strand"] = 1
        with pytest.raises(TypeError, match="strand must be a string"):
            CandidateResult(**kwargs)

    def test_invalid_strand_value(self):
        kwargs = self._base_kwargs()
        kwargs["strand"] = "forward"
        with pytest.raises(ValueError, match=r"strand must be '\+' or '-'"):
            CandidateResult(**kwargs)

    # Protospacer & PAM tests
    def test_invalid_protospacer_type(self):
        kwargs = self._base_kwargs()
        kwargs["protospacer"] = None
        with pytest.raises(TypeError, match="protospacer must be a string"):
            CandidateResult(**kwargs)

    def test_empty_protospacer(self):
        kwargs = self._base_kwargs()
        kwargs["protospacer"] = ""
        with pytest.raises(ValueError, match="protospacer cannot be empty"):
            CandidateResult(**kwargs)

    def test_invalid_pam_type(self):
        kwargs = self._base_kwargs()
        kwargs["pam"] = ["NGG"]
        with pytest.raises(TypeError, match="pam must be a string"):
            CandidateResult(**kwargs)

    def test_empty_pam(self):
        kwargs = self._base_kwargs()
        kwargs["pam"] = ""
        with pytest.raises(ValueError, match="pam cannot be empty"):
            CandidateResult(**kwargs)

    # Mismatch count & positions tests
    def test_mismatch_count_non_int(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = "zero"
        with pytest.raises(TypeError, match="mismatch_count must be an integer"):
            CandidateResult(**kwargs)

    def test_mismatch_count_negative(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = -1
        with pytest.raises(ValueError, match="mismatch_count must be non-negative"):
            CandidateResult(**kwargs)

    def test_mismatch_count_exceeds_length(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = 25  # Protospacer is 20-nt
        with pytest.raises(ValueError, match="cannot exceed protospacer length"):
            CandidateResult(**kwargs)

    def test_mismatch_positions_non_sequence(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_positions"] = 0
        with pytest.raises(TypeError, match="mismatch_positions must be a tuple or list"):
            CandidateResult(**kwargs)

    def test_mismatch_position_element_non_int(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = 1
        kwargs["mismatch_positions"] = ["1"]
        with pytest.raises(TypeError, match="Each mismatch position must be an integer"):
            CandidateResult(**kwargs)

    def test_mismatch_position_out_of_bounds_negative(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = 1
        kwargs["mismatch_positions"] = (-1,)
        with pytest.raises(ValueError, match="out of bounds"):
            CandidateResult(**kwargs)

    def test_mismatch_position_out_of_bounds_high(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = 1
        kwargs["mismatch_positions"] = (20,)  # Protospacer length is 20, max index 19
        with pytest.raises(ValueError, match="out of bounds"):
            CandidateResult(**kwargs)

    def test_mismatch_count_mismatch_positions_inconsistency(self):
        kwargs = self._base_kwargs()
        kwargs["mismatch_count"] = 2
        kwargs["mismatch_positions"] = (3,)
        with pytest.raises(ValueError, match="does not match mismatch_count"):
            CandidateResult(**kwargs)

    # Score tests
    def test_score_invalid_type_str(self):
        kwargs = self._base_kwargs()
        kwargs["score"] = "high"
        with pytest.raises(TypeError, match="score must be a numeric value"):
            CandidateResult(**kwargs)

    def test_score_boolean_rejected(self):
        kwargs = self._base_kwargs()
        kwargs["score"] = True
        with pytest.raises(TypeError, match="score must be a numeric value"):
            CandidateResult(**kwargs)

