"""Checks on complement, transcription, and translation."""

import random

import pytest

from revcomp.sequence import (
    CodingSequenceError,
    InvalidCharacter,
    back_transcribe,
    complement,
    reverse,
    reverse_complement,
    transcribe,
    translate,
    translate_frames,
    validate,
)


def test_reverse_complement_known_sequence():
    assert reverse_complement("ATGGCCATTGTAATG") == "CATTACAATGGCCAT"


def test_reverse_complement_is_an_involution():
    rng = random.Random(0)
    for _ in range(50):
        seq = "".join(rng.choice("ACGTRYN") for _ in range(rng.randint(1, 60)))
        assert reverse_complement(reverse_complement(seq)) == seq


def test_complement_does_not_reverse():
    assert complement("AAAC") == "TTTG"
    assert reverse_complement("AAAC") == "GTTT"


def test_complement_handles_ambiguity_codes():
    assert reverse_complement("RYSWKMBDHVN") == "NBDHVKMWSRY"
    assert complement("RYKM") == "YRMK"


def test_case_and_gaps_are_preserved():
    assert reverse_complement("atgGCC") == "GGCcat"
    assert complement("ACG-T.A") == "TGC-A.T"


def test_rna_output_mode():
    assert reverse_complement("ATGGCC", rna=True) == "GGCCAU"
    assert complement("ATGGCC", rna=True) == "UACCGG"


def test_reverse_only():
    assert reverse("ATGC") == "CGTA"
    with pytest.raises(InvalidCharacter):
        reverse("ATGX")


def test_sequence_errors_stay_hashable():
    """Exceptions must be usable in sets and as dict keys."""
    error = InvalidCharacter(positions=(1,), characters=("X",), alphabet_name="IUPAC DNA")
    assert hash(error)


def test_gaps_can_be_rejected():
    assert validate("AC-GT")[0] == "AC-GT"
    with pytest.raises(InvalidCharacter):
        validate("AC-GT", allow_gaps=False)


def test_whitespace_is_ignored():
    assert reverse_complement("ATG GCC\nATT") == reverse_complement("ATGGCCATT")


def test_t_and_u_are_accepted_in_the_same_alphabet():
    """There is one nucleotide alphabet, so neither base is rejected."""
    assert validate("ACGTU")[0] == "ACGTU"
    assert translate("AUGGCCATT") == "MAI"


def test_invalid_character_reports_position_and_character():
    with pytest.raises(InvalidCharacter) as info:
        reverse_complement("ATGXCC1")
    message = str(info.value)
    assert "'X' at position 4" in message
    assert "'1' at position 7" in message
    assert info.value.positions == (4, 7)


def test_invalid_character_message_truncates_long_lists():
    with pytest.raises(InvalidCharacter) as info:
        validate("ZZZZZZZZ")
    assert "and 3 more" in str(info.value)


def test_permissive_replaces_unknown_characters_with_n():
    cleaned, replaced = validate("ATGXcc", permissive=True)
    assert cleaned == "ATGNcc"
    assert replaced == (4,)


def test_transcription_of_the_coding_strand():
    assert transcribe("ATGGCCATTGTAATG") == "AUGGCCAUUGUAAUG"


def test_transcription_of_the_template_strand_is_reverse_complemented():
    """The transcript is complementary and antiparallel to the template."""
    coding = "ATGGCCATTGTAATG"
    template = reverse_complement(coding)
    assert transcribe(template, strand="template") == transcribe(coding)


def test_transcribe_rejects_an_unknown_strand():
    with pytest.raises(ValueError, match="coding"):
        transcribe("ATGC", strand="sense")


def test_back_transcription_round_trips():
    assert back_transcribe("AUGGCC") == "ATGGCC"
    assert back_transcribe(transcribe("ATGGCCATT")) == "ATGGCCATT"


def test_translation_of_a_known_coding_sequence():
    """The human preproinsulin signal peptide, translated in frame 1."""
    cds = (
        "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC"
        "CCAGCCGCAGCC"
    )
    assert translate(cds) == "MALWMRLLPLLALLALWGPDPAAA"


def test_translation_accepts_rna_and_lowercase():
    assert translate("AUGGCCAUU") == translate("atggccatt") == "MAI"


def test_stop_codons_and_to_stop():
    assert translate("ATGTAAATG") == "M*M"
    assert translate("ATGTAAATG", to_stop=True) == "M"
    assert translate("ATGTAAATG", stop_symbol="#") == "M#M"


def test_translation_uses_the_requested_table():
    """TGA is a stop and AGA arginine in table 1; table 2 swaps both."""
    assert translate("ATGTGAAGA") == "M*R"
    assert translate("ATGTGAAGA", table=2) == "MW*"
    assert translate("ATGTGAAGA", table="vertebrate-mitochondrial") == "MW*"


def test_ambiguous_codons_resolve_only_when_decidable():
    assert translate("GGN") == "G"  # all four glycine codons
    assert translate("TTR") == "L"  # TTA and TTG are both leucine
    assert translate("AAY") == "N"  # AAT and AAC are both asparagine
    assert translate("NNN") == "X"
    assert translate("TTN") == "X"  # TTT/TTC are F but TTA/TTG are L
    assert translate("TRA") == "*"  # TAA and TGA are both stops in table 1
    assert translate("TRG") == "X"  # TAG is a stop but TGG is tryptophan
    assert translate("TRA", table=2) == "X"  # TGA is tryptophan in table 2


def test_gapped_codons_translate_to_unknown():
    assert translate("ATG---ATT") == "M" + "X" + "I"


def test_trailing_partial_codon_is_dropped():
    assert translate("ATGGCCA") == "MA"


def test_cds_translation_forces_the_initiator_and_drops_the_stop():
    assert translate("ATGGCCCTGTGGTGA", cds=True) == "MALW"
    # CTG codes for leucine internally but initiates as methionine.
    assert translate("CTGGCCTAA", cds=True) == "MA"
    assert translate("CTGGCCTAA") == "LA*"


def test_cds_translation_respects_the_table_of_start_codons():
    assert translate("GTGGCCTAA", table=11, cds=True) == "MA"
    with pytest.raises(CodingSequenceError, match="initiation codon"):
        translate("GTGGCCTAA", table=1, cds=True)


@pytest.mark.parametrize(
    "seq,message",
    [
        ("ATGGCCCTGTGGTG", "multiple of three"),
        ("GGGGCCTAA", "initiation codon"),
        ("ATGGCCCTG", "does not end in a stop"),
        ("ATGTGAGCCTAA", "internal stop codon"),
        ("", "empty"),
    ],
)
def test_cds_translation_rejects_what_is_not_a_cds(seq, message):
    with pytest.raises(CodingSequenceError, match=message):
        translate(seq, cds=True)


def test_cds_and_to_stop_are_mutually_exclusive():
    with pytest.raises(ValueError, match="cannot both be set"):
        translate("ATGGCCTAA", cds=True, to_stop=True)


def test_frames_cover_both_strands():
    seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    results = translate_frames(seq, [1, 2, 3, -1, -2, -3])
    assert [item.frame for item in results] == [1, 2, 3, -1, -2, -3]
    assert results[0].protein == translate(seq)
    assert results[1].protein == translate(seq[1:])
    assert results[3].protein == translate(reverse_complement(seq))
    assert results[5].protein == translate(reverse_complement(seq)[2:])
    assert [item.strand for item in results] == ["+", "+", "+", "-", "-", "-"]


def test_frames_record_trailing_bases():
    seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"  # 39 bases
    results = translate_frames(seq, [1, 2, 3])
    assert [item.trailing_bases for item in results] == [0, 2, 1]


def test_translate_frames_rejects_a_bad_frame():
    with pytest.raises(ValueError, match="frame must be one of"):
        translate_frames("ATGC", [4])
