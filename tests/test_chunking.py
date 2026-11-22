
import pytest

from src.project.chunking import (
    chunk_sequence,
    scan_chunked_sequence
)


def test_chunk_creation():

    dna = "A" * 100

    chunks = chunk_sequence(
        "test",
        dna,
        chunk_size=40,
        overlap=22
    )

    assert len(chunks) == 4
    assert chunks[0].start == 0
    assert chunks[1].start == 18


def test_chunk_coordinates():

    dna = "A" * 100

    chunks = chunk_sequence(
        "test",
        dna,
        chunk_size=40,
        overlap=22
    )

    for chunk in chunks:
        assert len(chunk.sequence) == chunk.end - chunk.start


def test_empty_sequence():

    chunks = chunk_sequence("test", "")

    assert chunks == []


def test_invalid_chunk_size():

    with pytest.raises(ValueError):
        chunk_sequence("test", "ACGT", chunk_size=0)


def test_invalid_overlap():

    with pytest.raises(ValueError):
        chunk_sequence(
            "test", "ACGT",
            chunk_size=20,
            overlap=20
        )


def test_pam_across_chunk_boundary():

    # Candidate begins at index 3 and spans 23 bases.
    dna = "CCC" + ("A" * 20) + "AGG" + ("C" * 20)

    candidates = scan_chunked_sequence(
        "test",
        dna,
        chunk_size=24,
        overlap=22
    )

    matching = [
        candidate for candidate in candidates
        if candidate.start == 3
        and candidate.strand == "+"
    ]

    assert len(matching) == 1
    assert matching[0].pam == "AGG"
