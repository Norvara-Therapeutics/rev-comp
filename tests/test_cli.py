"""End-to-end checks on the command line interface."""

import io

import pytest

from revcomp import cli


def run(argv, stdin=""):
    """Invoke the CLI with captured streams; returns (code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    code = cli.main(argv, stdin=io.StringIO(stdin), stdout=out, stderr=err)
    return code, out.getvalue(), err.getvalue()


def sequences(fasta_text):
    """Pull the residues out of FASTA text, one string per record."""
    from revcomp.fasta import parse

    return [record.sequence for record in parse(fasta_text)]


def test_complement_from_stdin():
    code, out, err = run(["complement"], stdin="ATGGCCATTGTAATG\n")
    assert code == cli.EXIT_OK
    assert sequences(out) == ["CATTACAATGGCCAT"]
    assert out.startswith(">sequence reverse-complement")
    assert err == ""


def test_complement_from_a_positional_argument():
    code, out, _ = run(["complement", "ATGGCC"])
    assert code == cli.EXIT_OK
    assert sequences(out) == ["GGCCAT"]


def test_complement_modes():
    assert sequences(run(["complement", "--mode", "complement", "AAAC"])[1]) == ["TTTG"]
    assert sequences(run(["complement", "--mode", "reverse", "AAAC"])[1]) == ["CAAA"]


def test_complement_preserves_fasta_names_and_descriptions():
    fasta = ">primer_F forward, 6 nt\nATGGCC\n"
    code, out, _ = run(["complement"], stdin=fasta)
    assert code == cli.EXIT_OK
    assert out.splitlines()[0] == ">primer_F forward, 6 nt reverse-complement"


def test_transcribe_coding_and_template_strands():
    assert sequences(run(["transcribe", "ATGGCC"])[1]) == ["AUGGCC"]
    assert sequences(run(["transcribe", "--strand", "template", "GGCCAT"])[1]) == [
        "AUGGCC"
    ]
    assert sequences(run(["transcribe", "--back", "AUGGCC"])[1]) == ["ATGGCC"]


def test_translate_single_frame():
    code, out, err = run(["translate", "--frame", "1", "ATGGCCATTGTAATGGGC"])
    assert code == cli.EXIT_OK
    assert sequences(out) == ["MAIVMG"]
    assert "frame=+1 table=1" in out


def test_translate_all_six_frames():
    code, out, _ = run(["translate", "--frame", "all", "ATGGCCATTGTAATGGGCCGCTGA"])
    assert code == cli.EXIT_OK
    assert len(sequences(out)) == 6
    assert "frame=-3" in out


def test_rna_output_applies_in_every_complement_mode():
    assert sequences(run(["complement", "--rna", "ATGGCC"])[1]) == ["GGCCAU"]
    assert sequences(run(["complement", "--rna", "-m", "complement", "ATGGCC"])[1]) == [
        "UACCGG"
    ]
    assert sequences(run(["complement", "--rna", "-m", "reverse", "ATGGCC"])[1]) == [
        "CCGGUA"
    ]


def test_translate_accepts_negative_and_comma_separated_frames():
    code, out, _ = run(["translate", "--frame=-1,2", "ATGGCCATTGTAATGGGC"])
    assert code == cli.EXIT_OK
    assert len(sequences(out)) == 2
    assert "frame=-1" in out and "frame=+2" in out


def test_translated_records_get_unique_identifiers():
    """Six frames under one name would collide in any id-keyed tool."""
    code, out, _ = run(["translate", "--frame", "all", "ATGGCCATTGTAATGGGC"])
    names = [line.split()[0] for line in out.splitlines() if line.startswith(">")]
    assert names == [
        ">sequence_f1",
        ">sequence_f2",
        ">sequence_f3",
        ">sequence_r1",
        ">sequence_r2",
        ">sequence_r3",
    ]
    assert len(set(names)) == 6


def test_cds_mode_checks_and_trims_a_coding_sequence():
    code, out, err = run(["translate", "--cds", "ATGGCCCTGTGGTGA"])
    assert code == cli.EXIT_OK
    assert sequences(out) == ["MALW"]  # terminal stop omitted
    assert err == ""


def test_cds_mode_reports_what_is_wrong_and_names_the_record():
    cases = {
        ">a\nATGGCCCTGTGGTG\n": "not a multiple of three",
        ">a\nGGGGCCTAA\n": "not an initiation codon",
        ">a\nATGGCCCTG\n": "does not end in a stop codon",
        ">a\nATGTGAGCCTAA\n": "internal stop codon at codon 2 of 4",
    }
    for fasta, expected in cases.items():
        code, out, err = run(["translate", "--cds"], stdin=fasta)
        assert code == cli.EXIT_SEQUENCE
        assert expected in err
        assert err.startswith("rev-comp: error: a: ")
        assert out == ""


def test_cds_mode_conflicts():
    for argv in (
        ["translate", "--cds", "--to-stop", "ATGGCCTGA"],
        ["translate", "--cds", "--frame", "all", "ATGGCCTGA"],
    ):
        code, _, err = run(argv)
        assert code == cli.EXIT_ERROR
        assert "--cds" in err


def test_mixed_t_and_u_is_warned_about():
    code, out, err = run(["complement", "ATGCUUGA"])
    assert code == cli.EXIT_OK
    assert "contains both T and U" in err
    assert sequences(out) == ["TCAAGCAT"]


def test_translate_honours_the_genetic_code():
    assert sequences(run(["translate", "ATGTGAAGA"])[1]) == ["M*R"]
    assert sequences(run(["translate", "--table", "2", "ATGTGAAGA"])[1]) == ["MW*"]
    assert sequences(run(["translate", "-t", "bacterial", "ATGTGAAGA"])[1]) == ["M*R"]


def test_translate_warns_once_that_the_input_is_not_whole_codons():
    """Offset frames always leave bases over; the input length is the signal."""
    code, out, err = run(["translate", "ATGGCCA"])
    assert code == cli.EXIT_OK
    assert sequences(out) == ["MA"]
    assert err.count("not a multiple of 3") == 1

    code, out, err = run(["translate", "--frame", "all", "ATGGCCATTGTAATGGGC"])
    assert code == cli.EXIT_OK
    assert err == ""  # length 18 is whole codons, offset frames aside


def test_quiet_suppresses_warnings():
    code, _, err = run(["translate", "--quiet", "ATGGCCA"])
    assert code == cli.EXIT_OK
    assert err == ""


def test_invalid_characters_fail_with_a_readable_message():
    code, out, err = run(["complement", "ATGXCC"])
    assert code == cli.EXIT_SEQUENCE
    assert "'X' at position 4" in err
    assert "--permissive" in err
    assert out == ""


def test_permissive_mode_warns_and_continues():
    code, out, err = run(["complement", "--permissive", "ATGXCC"])
    assert code == cli.EXIT_OK
    assert sequences(out) == ["GGNCAT"]
    assert "replaced 1 character" in err


def test_malformed_fasta_is_reported():
    code, _, err = run(["complement"], stdin="ACGT\n>late\nACGT\n")
    assert code == cli.EXIT_SEQUENCE
    assert "before any" in err


def test_sequence_and_input_cannot_both_be_given(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">a\nACGT\n")
    code, _, err = run(["complement", "ACGT", "--input", str(path)])
    assert code == cli.EXIT_ERROR
    assert "not both" in err


def test_missing_input_file_is_reported():
    code, _, err = run(["complement", "--input", "/nonexistent/path.fasta"])
    assert code == cli.EXIT_ERROR
    assert "cannot read" in err


def test_file_input_and_output(tmp_path):
    source = tmp_path / "primer.fasta"
    source.write_text(">primer_F\nATGGCCATTGTAATG\n")
    destination = tmp_path / "rev_primer.fasta"
    code, out, _ = run(
        ["complement", "--input", str(source), "--output", str(destination)]
    )
    assert code == cli.EXIT_OK
    assert out == ""
    assert sequences(destination.read_text()) == ["CATTACAATGGCCAT"]


def test_line_width_option(tmp_path):
    code, out, _ = run(["complement", "--line-width", "10", "A" * 25])
    body = [line for line in out.splitlines() if not line.startswith(">")]
    assert [len(line) for line in body] == [10, 10, 5]


def test_multi_record_input_is_processed_record_by_record():
    fasta = ">a\nATG\n>b\nGGG\n"
    code, out, _ = run(["complement"], stdin=fasta)
    assert sequences(out) == ["CAT", "CCC"]


def test_tables_lists_every_code():
    code, out, _ = run(["tables"])
    assert code == cli.EXIT_OK
    assert "Standard" in out
    assert "Vertebrate Mitochondrial" in out
    assert len(out.splitlines()) == 27


def test_tables_shows_one_code_in_detail():
    code, out, _ = run(["tables", "--table", "2"])
    assert code == cli.EXIT_OK
    assert "NCBI table 2" in out
    assert "start codons" in out
    assert "TGA W" in out


def test_unknown_table_is_a_usage_error():
    with pytest.raises(SystemExit) as info:
        run(["translate", "--table", "99", "ATG"])
    assert info.value.code == 2


def test_bad_frame_is_a_usage_error():
    with pytest.raises(SystemExit) as info:
        run(["translate", "--frame", "9", "ATG"])
    assert info.value.code == 2


def test_stop_and_unknown_symbols_must_be_one_character():
    for argv in (
        ["translate", "--stop-symbol", "STOP", "ATGTGA"],
        ["translate", "--unknown", "???", "ATGNNN"],
    ):
        with pytest.raises(SystemExit) as info:
            run(argv)
        assert info.value.code == 2


def test_negative_line_width_is_rejected():
    with pytest.raises(SystemExit) as info:
        run(["complement", "--line-width", "-5", "ACGT"])
    assert info.value.code == 2


def test_back_transcription_cannot_be_combined_with_a_strand():
    code, out, err = run(["transcribe", "--back", "--strand", "template", "AUGGCC"])
    assert code == cli.EXIT_ERROR
    assert "--back cannot be combined with --strand" in err
    assert out == ""


def test_empty_record_is_written_as_a_bare_header():
    code, out, err = run(["complement"], stdin=">a\n>b\nACGT\n")
    assert code == cli.EXIT_OK
    assert out.splitlines() == [">a reverse-complement", ">b reverse-complement", "ACGT"]
    assert "result is empty" in err


def test_a_closed_pipe_is_not_an_error():
    class ClosedPipe(io.StringIO):
        def write(self, text):  # noqa: D102 - stands in for `| head`
            raise BrokenPipeError(32, "Broken pipe")

    err = io.StringIO()
    code = cli.main(
        ["complement", "ACGT"],
        stdin=io.StringIO(),
        stdout=ClosedPipe(),
        stderr=err,
    )
    assert code == cli.EXIT_OK
    assert err.getvalue() == ""


def test_no_command_prints_help():
    code, out, _ = run([])
    assert code == cli.EXIT_OK
    assert "complement" in out and "translate" in out
