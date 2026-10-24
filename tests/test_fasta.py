"""Checks on FASTA parsing and writing."""

import io

import pytest

from revcomp.fasta import FastaError, Record, format_record, parse, read_sequences, write

SAMPLE = """>primer_F forward primer, 20 nt
ACGTACGTAC
GTACGTACGT

>primer_R
TTTTGGGGCC
"""


def test_parse_multiple_records():
    records = list(parse(SAMPLE))
    assert len(records) == 2
    assert records[0].name == "primer_F"
    assert records[0].description == "forward primer, 20 nt"
    assert records[0].sequence == "ACGTACGTACGTACGTACGT"
    assert records[1].name == "primer_R"
    assert records[1].description == ""
    assert records[1].sequence == "TTTTGGGGCC"


def test_parse_accepts_a_file_handle():
    records = list(parse(io.StringIO(SAMPLE)))
    assert [record.name for record in records] == ["primer_F", "primer_R"]


def test_parse_rejects_sequence_before_a_header():
    with pytest.raises(FastaError, match="before any"):
        list(parse("ACGT\n>late_header\nACGT\n"))


def test_parse_rejects_an_empty_header():
    with pytest.raises(FastaError, match="no identifier"):
        list(parse(">\nACGT\n"))


def test_read_sequences_accepts_a_bare_sequence():
    records = read_sequences("ACG TAC\nGT", default_name="oligo")
    assert len(records) == 1
    assert records[0].name == "oligo"
    assert records[0].sequence == "ACGTACGT"


def test_read_sequences_rejects_empty_input():
    with pytest.raises(FastaError, match="no input"):
        read_sequences("   \n")


def test_format_record_wraps_at_the_requested_width():
    record = Record(name="x", sequence="A" * 25, description="note")
    text = format_record(record, line_width=10)
    assert text.splitlines() == [">x note", "A" * 10, "A" * 10, "A" * 5]


def test_format_record_can_write_one_long_line():
    record = Record(name="x", sequence="A" * 25)
    assert format_record(record, line_width=0).splitlines() == [">x", "A" * 25]


def test_format_record_writes_an_empty_record_as_a_bare_header():
    assert format_record(Record(name="x", sequence="")) == ">x\n"


def test_write_returns_the_number_of_records():
    handle = io.StringIO()
    count = write(list(parse(SAMPLE)), handle, line_width=60)
    assert count == 2
    assert handle.getvalue().count(">") == 2


def test_round_trip_through_write_and_parse():
    original = list(parse(SAMPLE))
    handle = io.StringIO()
    write(original, handle)
    handle.seek(0)
    assert list(parse(handle)) == original
