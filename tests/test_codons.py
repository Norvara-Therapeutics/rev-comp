"""Checks on the bundled NCBI translation tables."""

import pytest

from revcomp.codons import CODONS, available_tables, get_genetic_code, resolve_table_name


def test_every_table_assigns_all_64_codons():
    for code in available_tables():
        assert len(code.codons) == 64
        assert set(code.codons) == set(CODONS)


def test_standard_code_assignments():
    code = get_genetic_code(1)
    assert code.id == 1
    assert code.name == "Standard"
    assert code.codons["ATG"] == "M"
    assert code.codons["TTT"] == "F"
    assert code.codons["TGG"] == "W"
    assert code.stop_codons == {"TAA", "TAG", "TGA"}
    assert code.start_codons == {"TTG", "CTG", "ATG"}


def test_vertebrate_mitochondrial_differences():
    """Table 2 reassigns four codons relative to the standard code."""
    standard = get_genetic_code(1)
    mito = get_genetic_code(2)
    differences = {
        codon
        for codon in CODONS
        if standard.codons[codon] != mito.codons[codon]
    }
    assert differences == {"AGA", "AGG", "ATA", "TGA"}
    assert mito.codons["TGA"] == "W"
    assert mito.codons["ATA"] == "M"
    assert mito.codons["AGA"] == "*"
    assert mito.stop_codons == {"TAA", "TAG", "AGA", "AGG"}


def test_bacterial_table_has_extra_start_codons():
    code = get_genetic_code(11)
    assert code.codons == get_genetic_code(1).codons  # same assignments
    assert {"ATT", "ATC", "ATA", "GTG", "TTG"} <= code.start_codons
    assert not code.start_codons & code.stop_codons


def test_ciliate_table_reads_through_two_stops():
    code = get_genetic_code(6)
    assert code.codons["TAA"] == "Q"
    assert code.codons["TAG"] == "Q"
    assert code.stop_codons == {"TGA"}


@pytest.mark.parametrize(
    "text,expected",
    [
        ("1", 1),
        ("standard", 1),
        ("STANDARD", 1),
        ("bacterial", 11),
        ("vertebrate-mitochondrial", 2),
        ("vertebrate mitochondrial", 2),
        ("Euplotid Nuclear", 10),
        ("2", 2),
    ],
)
def test_resolve_table_name(text, expected):
    assert resolve_table_name(text) == expected


@pytest.mark.parametrize("text", ["0", "7", "99", "nonsense", ""])
def test_resolve_table_name_rejects_unknown(text):
    with pytest.raises(KeyError):
        resolve_table_name(text)


def test_get_genetic_code_rejects_unknown_id():
    with pytest.raises(KeyError):
        get_genetic_code(7)
