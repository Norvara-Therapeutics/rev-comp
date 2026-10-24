"""NCBI genetic code tables and IUPAC nucleotide alphabets.

The codon tables are transcribed verbatim from the NCBI translation-table
definition file (``ncbieaa``/``sncbieaa`` strings):

    https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt

Downloaded 2026-10-21; sha256 of the source file:
2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04

Codons are enumerated in NCBI order: base 1 cycles slowest and base 3 fastest
over ``T, C, A, G``, so the *n*-th character of ``ncbieaa`` is the amino acid
for the *n*-th codon of that enumeration. ``sncbieaa`` marks initiation codons
with the amino acid they are translated as when they initiate (``M`` in nearly
every table), stop codons with ``*``, and everything else with ``-``; only the
amino-acid marks denote starts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, FrozenSet, Mapping, Tuple

__all__ = [
    "BASES",
    "CODONS",
    "STOP",
    "GAP_CHARACTERS",
    "IUPAC_NUCLEOTIDES",
    "IUPAC_EXPANSION",
    "DNA_COMPLEMENT",
    "RNA_COMPLEMENT",
    "GeneticCode",
    "available_tables",
    "get_genetic_code",
    "resolve_table_name",
]

#: Unambiguous bases, in the order NCBI enumerates codons.
BASES: Tuple[str, ...] = ("T", "C", "A", "G")

#: All 64 codons in NCBI enumeration order.
CODONS: Tuple[str, ...] = tuple(
    b1 + b2 + b3 for b1 in BASES for b2 in BASES for b3 in BASES
)

#: Symbol used for a translated stop codon.
STOP = "*"

#: Characters treated as alignment gaps rather than sequence errors.
GAP_CHARACTERS: FrozenSet[str] = frozenset("-.~")

#: Nucleotide codes accepted in input: the four bases, U, and the IUPAC
#: ambiguity codes. T and U are both accepted whichever the input is meant to
#: be, so there is one alphabet rather than separate DNA and RNA ones.
IUPAC_NUCLEOTIDES: FrozenSet[str] = frozenset("ACGTURYSWKMBDHVN")

#: Expansion of each IUPAC ambiguity code into the bases it may represent.
IUPAC_EXPANSION: Mapping[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

def _with_case(pairs: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for src_base, dst_base in zip(pairs[::2], pairs[1::2]):
        mapping[src_base] = dst_base
        mapping[src_base.lower()] = dst_base.lower()
    for gap in GAP_CHARACTERS:
        mapping[gap] = gap
    return mapping

#: Watson-Crick complement over the full IUPAC alphabet, DNA output.
DNA_COMPLEMENT: Mapping[str, str] = _with_case(
    "ATTACGGCUARYYRSSWWKMMKBVDHHDVBNN"
)

#: Watson-Crick complement over the full IUPAC alphabet, RNA output.
RNA_COMPLEMENT: Mapping[str, str] = _with_case(
    "AUTACGGCUARYYRSSWWKMMKBVDHHDVBNN"
)


@dataclass(frozen=True)
class _Spec:
    id: int
    names: Tuple[str, ...]
    ncbieaa: str
    sncbieaa: str


_SPECS: Dict[int, _Spec] = {
    1: _Spec(
        id=1,
        names=('Standard', 'SGC0',),
        ncbieaa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='---M------**--*----M---------------M----------------------------',
    ),
    2: _Spec(
        id=2,
        names=('Vertebrate Mitochondrial', 'SGC1',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
        sncbieaa='----------**--------------------MMMM----------**---M------------',
    ),
    3: _Spec(
        id=3,
        names=('Yeast Mitochondrial', 'SGC2',),
        ncbieaa='FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**----------------------MM---------------M------------',
    ),
    4: _Spec(
        id=4,
        names=('Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma', 'SGC3',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--MM------**-------M------------MMMM---------------M------------',
    ),
    5: _Spec(
        id=5,
        names=('Invertebrate Mitochondrial', 'SGC4',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
        sncbieaa='---M------**--------------------MMMM---------------M------------',
    ),
    6: _Spec(
        id=6,
        names=('Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear', 'SGC5',),
        ncbieaa='FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--------------*--------------------M----------------------------',
    ),
    9: _Spec(
        id=9,
        names=('Echinoderm Mitochondrial; Flatworm Mitochondrial', 'SGC8',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        sncbieaa='----------**-----------------------M---------------M------------',
    ),
    10: _Spec(
        id=10,
        names=('Euplotid Nuclear', 'SGC9',),
        ncbieaa='FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**-----------------------M----------------------------',
    ),
    11: _Spec(
        id=11,
        names=('Bacterial, Archaeal and Plant Plastid',),
        ncbieaa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='---M------**--*----M------------MMMM---------------M------------',
    ),
    12: _Spec(
        id=12,
        names=('Alternative Yeast Nuclear',),
        ncbieaa='FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**--*----M---------------M----------------------------',
    ),
    13: _Spec(
        id=13,
        names=('Ascidian Mitochondrial',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
        sncbieaa='---M------**----------------------MM---------------M------------',
    ),
    14: _Spec(
        id=14,
        names=('Alternative Flatworm Mitochondrial',),
        ncbieaa='FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        sncbieaa='-----------*-----------------------M----------------------------',
    ),
    15: _Spec(
        id=15,
        names=('Blepharisma Macronuclear',),
        ncbieaa='FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------*---*--------------------M----------------------------',
    ),
    16: _Spec(
        id=16,
        names=('Chlorophycean Mitochondrial',),
        ncbieaa='FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------*---*--------------------M----------------------------',
    ),
    21: _Spec(
        id=21,
        names=('Trematode Mitochondrial',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        sncbieaa='----------**-----------------------M---------------M------------',
    ),
    22: _Spec(
        id=22,
        names=('Scenedesmus obliquus Mitochondrial',),
        ncbieaa='FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='------*---*---*--------------------M----------------------------',
    ),
    23: _Spec(
        id=23,
        names=('Thraustochytrium Mitochondrial',),
        ncbieaa='FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--*-------**--*-----------------M--M---------------M------------',
    ),
    24: _Spec(
        id=24,
        names=('Rhabdopleuridae Mitochondrial',),
        ncbieaa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
        sncbieaa='---M------**-------M---------------M---------------M------------',
    ),
    25: _Spec(
        id=25,
        names=('Candidate Division SR1 and Gracilibacteria',),
        ncbieaa='FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='---M------**-----------------------M---------------M------------',
    ),
    26: _Spec(
        id=26,
        names=('Pachysolen tannophilus Nuclear',),
        ncbieaa='FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**--*----M---------------M----------------------------',
    ),
    27: _Spec(
        id=27,
        names=('Karyorelict Nuclear',),
        ncbieaa='FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--------------*--------------------M----------------------------',
    ),
    28: _Spec(
        id=28,
        names=('Condylostoma Nuclear',),
        ncbieaa='FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**--*--------------------M----------------------------',
    ),
    29: _Spec(
        id=29,
        names=('Mesodinium Nuclear',),
        ncbieaa='FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--------------*--------------------M----------------------------',
    ),
    30: _Spec(
        id=30,
        names=('Peritrich Nuclear',),
        ncbieaa='FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='--------------*--------------------M----------------------------',
    ),
    31: _Spec(
        id=31,
        names=('Blastocrithidia Nuclear',),
        ncbieaa='FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='----------**-----------------------M----------------------------',
    ),
    32: _Spec(
        id=32,
        names=('Balanophoraceae Plastid',),
        ncbieaa='FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        sncbieaa='---M------*---*----M------------MMMM---------------M------------',
    ),
    33: _Spec(
        id=33,
        names=('Cephalodiscidae Mitochondrial',),
        ncbieaa='FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
        sncbieaa='---M-------*-------M---------------M---------------M------------',
    ),
}


@dataclass(frozen=True)
class GeneticCode:
    """One NCBI translation table.

    Attributes:
        id: NCBI translation table number (``1`` is the standard code).
        names: Names NCBI gives the table, most common first.
        codons: Mapping of each of the 64 unambiguous codons (DNA alphabet,
            uppercase) to a one-letter amino acid, with ``*`` for stop.
        start_codons: Codons annotated as initiation codons in this table.
        stop_codons: Codons translated as stop in this table.
    """

    id: int
    names: Tuple[str, ...]
    codons: Mapping[str, str]
    start_codons: FrozenSet[str]
    stop_codons: FrozenSet[str]

    @property
    def name(self) -> str:
        """The table's primary NCBI name."""
        return self.names[0]

    def __str__(self) -> str:
        return f"{self.id} ({self.name})"


def _build(spec: _Spec) -> GeneticCode:
    codons = dict(zip(CODONS, spec.ncbieaa))
    # In sncbieaa, '-' marks a non-initiator and '*' re-marks a stop codon;
    # only an amino-acid letter marks a genuine initiation codon.
    starts = frozenset(
        codon
        for codon, flag in zip(CODONS, spec.sncbieaa)
        if flag not in ("-", STOP)
    )
    stops = frozenset(codon for codon, aa in codons.items() if aa == STOP)
    return GeneticCode(
        id=spec.id,
        names=spec.names,
        codons=codons,
        start_codons=starts,
        stop_codons=stops,
    )


_CACHE: Dict[int, GeneticCode] = {}

#: Convenience aliases accepted wherever a table is named on the command line.
_ALIASES: Mapping[str, int] = {
    "standard": 1,
    "nuclear": 1,
    "universal": 1,
    "vertebrate-mitochondrial": 2,
    "vertebrate-mito": 2,
    "yeast-mitochondrial": 3,
    "mold-mitochondrial": 4,
    "invertebrate-mitochondrial": 5,
    "invertebrate-mito": 5,
    "ciliate-nuclear": 6,
    "echinoderm-mitochondrial": 9,
    "euplotid-nuclear": 10,
    "bacterial": 11,
    "archaeal": 11,
    "plant-plastid": 11,
    "alternative-yeast-nuclear": 12,
    "ascidian-mitochondrial": 13,
    "blepharisma-macronuclear": 15,
    "chlorophycean-mitochondrial": 16,
    "trematode-mitochondrial": 21,
}


def available_tables() -> Tuple[GeneticCode, ...]:
    """Every NCBI translation table bundled with this package, by id."""
    return tuple(get_genetic_code(table_id) for table_id in sorted(_SPECS))


def resolve_table_name(value: str) -> int:
    """Resolve a table id or alias given as text to an NCBI table number.

    Raises:
        KeyError: if ``value`` is neither a known id nor a known alias.
    """
    text = value.strip().lower().replace("_", "-").replace(" ", "-")
    if text.isdigit():
        table_id = int(text)
        if table_id not in _SPECS:
            raise KeyError(value)
        return table_id
    if text in _ALIASES:
        return _ALIASES[text]
    for spec in _SPECS.values():
        for name in spec.names:
            if name.lower().replace(" ", "-") == text:
                return spec.id
    raise KeyError(value)


def get_genetic_code(table: int | str = 1) -> GeneticCode:
    """Return a :class:`GeneticCode` by NCBI id, name, or alias.

    Raises:
        KeyError: if the table is unknown.
    """
    table_id = resolve_table_name(table) if isinstance(table, str) else int(table)
    if table_id not in _SPECS:
        raise KeyError(table_id)
    if table_id not in _CACHE:
        _CACHE[table_id] = _build(_SPECS[table_id])
    return _CACHE[table_id]
