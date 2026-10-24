"""Sequence operations: complement, transcription, and translation.

All functions are pure string transformations over the IUPAC nucleotide
alphabet and depend only on the standard library.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Iterable, List, Tuple

from revcomp.codons import (
    DNA_COMPLEMENT,
    GAP_CHARACTERS,
    IUPAC_EXPANSION,
    IUPAC_NUCLEOTIDES,
    RNA_COMPLEMENT,
    STOP,
    GeneticCode,
    get_genetic_code,
)

__all__ = [
    "SequenceError",
    "InvalidCharacter",
    "CodingSequenceError",
    "FRAMES",
    "validate",
    "complement",
    "reverse",
    "reverse_complement",
    "transcribe",
    "back_transcribe",
    "translate",
    "translate_frames",
    "FrameTranslation",
]

#: The six reading frames, forward strand first.
FRAMES: Tuple[int, ...] = (1, 2, 3, -1, -2, -3)

#: How many offending positions an error message lists before truncating.
_MAX_REPORTED = 5


class SequenceError(ValueError):
    """Base class for malformed input detected by this package."""


@dataclass(eq=False)
class InvalidCharacter(SequenceError):
    """One or more characters are outside the permitted alphabet.

    Attributes:
        positions: 1-based offsets of the offending characters, in order.
        characters: The offending characters, aligned with ``positions``.
        alphabet_name: Human-readable name of the expected alphabet.
        record: Name of the record the characters came from, if known.
    """

    positions: Tuple[int, ...]
    characters: Tuple[str, ...]
    alphabet_name: str
    record: str | None = None

    def __str__(self) -> str:
        shown = list(zip(self.positions, self.characters))[:_MAX_REPORTED]
        detail = ", ".join(f"{char!r} at position {pos}" for pos, char in shown)
        if len(self.positions) > _MAX_REPORTED:
            detail += f", and {len(self.positions) - _MAX_REPORTED} more"
        where = f" in {self.record!r}" if self.record else ""
        return (
            f"{len(self.positions)} character(s) outside the "
            f"{self.alphabet_name} alphabet{where}: {detail}"
        )


class CodingSequenceError(SequenceError):
    """Raised when ``cds`` translation finds the input is not a valid CDS."""


def validate(
    seq: str,
    *,
    allow_gaps: bool = True,
    permissive: bool = False,
    record: str | None = None,
) -> Tuple[str, Tuple[int, ...]]:
    """Check ``seq`` against the IUPAC nucleotide alphabet.

    Whitespace (including newlines) is always stripped. Characters outside the
    alphabet raise :class:`InvalidCharacter` unless ``permissive`` is set, in
    which case they are replaced by ``N`` and their positions are returned.

    Args:
        seq: Raw sequence text, possibly containing whitespace.
        allow_gaps: Treat ``-``, ``.`` and ``~`` as gaps rather than errors.
        permissive: Replace offending characters with ``N`` instead of raising.
        record: Record name used in error messages.

    Returns:
        The cleaned sequence and the 1-based residue positions that were
        replaced, counted after whitespace removal
        (always empty when ``permissive`` is false, since anything to replace
        would have raised).

    Raises:
        InvalidCharacter: If an unknown character is present and ``permissive``
            is false.
    """
    cleaned = "".join(seq.split())
    allowed = set(IUPAC_NUCLEOTIDES)
    if allow_gaps:
        allowed = allowed | set(GAP_CHARACTERS)
    positions: List[int] = []
    characters: List[str] = []
    out: List[str] = []
    for index, char in enumerate(cleaned, start=1):
        if char.upper() in allowed:
            out.append(char)
            continue
        positions.append(index)
        characters.append(char)
        out.append("n" if char.islower() else "N")
    if positions and not permissive:
        raise InvalidCharacter(
            positions=tuple(positions),
            characters=tuple(characters),
            alphabet_name="IUPAC nucleotide",
            record=record,
        )
    return "".join(out), tuple(positions)


def complement(seq: str, *, rna: bool = False) -> str:
    """Complement ``seq`` base by base, preserving order, case, and gaps.

    Args:
        seq: Nucleotide sequence over the IUPAC alphabet.
        rna: Emit ``U`` in place of ``T``.

    Raises:
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    cleaned, _ = validate(seq)
    table = RNA_COMPLEMENT if rna else DNA_COMPLEMENT
    try:
        return "".join(table[char] for char in cleaned)
    except KeyError as exc:  # pragma: no cover - validate covers the alphabet
        raise SequenceError(f"cannot complement {exc.args[0]!r}") from exc


def reverse(seq: str) -> str:
    """Reverse ``seq`` without complementing it.

    Raises:
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    cleaned, _ = validate(seq)
    return cleaned[::-1]


def reverse_complement(seq: str, *, rna: bool = False) -> str:
    """Return the reverse complement of ``seq``: the opposite strand, 5'->3'.

    Args:
        seq: Nucleotide sequence over the IUPAC alphabet.
        rna: Emit ``U`` in place of ``T``.

    Raises:
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    return complement(seq, rna=rna)[::-1]


def transcribe(seq: str, *, strand: str = "coding") -> str:
    """Transcribe DNA to RNA.

    Args:
        seq: DNA sequence over the IUPAC alphabet.
        strand: Which strand ``seq`` represents. ``"coding"`` (sense) is
            copied base for base with ``T`` replaced by ``U``; ``"template"``
            (antisense) is reverse complemented first, since the transcript is
            complementary and antiparallel to it.

    Raises:
        ValueError: If ``strand`` is not ``"coding"`` or ``"template"``.
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    if strand not in ("coding", "template"):
        raise ValueError(f"strand must be 'coding' or 'template', not {strand!r}")
    if strand == "template":
        return reverse_complement(seq, rna=True)
    cleaned, _ = validate(seq)
    return cleaned.replace("T", "U").replace("t", "u")


def back_transcribe(seq: str) -> str:
    """Convert an RNA sequence back to DNA by replacing ``U`` with ``T``.

    Raises:
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    cleaned, _ = validate(seq)
    return cleaned.replace("U", "T").replace("u", "t")


def _codon_amino_acid(codon: str, code: GeneticCode, unknown: str) -> str:
    """Translate one codon, resolving IUPAC ambiguity where it is decidable."""
    codon = codon.upper().replace("U", "T")
    direct = code.codons.get(codon)
    if direct is not None:
        return direct
    try:
        options = [IUPAC_EXPANSION[char] for char in codon]
    except KeyError:
        return unknown
    amino_acids = {
        code.codons["".join(combo)] for combo in itertools.product(*options)
    }
    if len(amino_acids) == 1:
        return amino_acids.pop()
    return unknown


def _translate_cds(seq: str, code: GeneticCode, unknown: str) -> str:
    """Translate a complete coding sequence, checking that it is one.

    The initiator is translated as methionine whatever it codes for elsewhere,
    and the terminal stop is not written.

    Raises:
        CodingSequenceError: If the sequence is not a whole number of codons,
            does not begin with an initiation codon of ``code``, contains an
            internal stop, or does not end in a stop.
    """
    if not seq:
        raise CodingSequenceError("the coding sequence is empty")
    if len(seq) % 3:
        raise CodingSequenceError(
            f"length {len(seq)} is not a multiple of three"
        )
    first = seq[:3].upper().replace("U", "T")
    if first not in code.start_codons:
        starts = ", ".join(sorted(code.start_codons))
        raise CodingSequenceError(
            f"{first} is not an initiation codon in table {code}; "
            f"this table starts at {starts}"
        )
    residues = ["M"]
    for start in range(3, len(seq), 3):
        amino_acid = _codon_amino_acid(seq[start : start + 3], code, unknown)
        if amino_acid == STOP:
            if start + 3 == len(seq):
                return "".join(residues)
            raise CodingSequenceError(
                f"internal stop codon at codon {start // 3 + 1} of "
                f"{len(seq) // 3}"
            )
        residues.append(amino_acid)
    raise CodingSequenceError("the coding sequence does not end in a stop codon")


@dataclass(frozen=True)
class FrameTranslation:
    """A translation of one reading frame.

    Attributes:
        frame: Frame number, ``1``-``3`` forward and ``-1``-``-3`` reverse.
        protein: The translated amino acid sequence.
        trailing_bases: Bases left over at the 3' end that did not fill a
            codon (0, 1, or 2).
    """

    frame: int
    protein: str
    trailing_bases: int

    @property
    def strand(self) -> str:
        """``"+"`` for a forward frame, ``"-"`` for a reverse one."""
        return "+" if self.frame > 0 else "-"


def translate(
    seq: str,
    *,
    table: int | str = 1,
    cds: bool = False,
    to_stop: bool = False,
    stop_symbol: str = STOP,
    unknown: str = "X",
) -> str:
    """Translate a nucleotide sequence in frame 1.

    Codons containing ambiguity codes are translated when every base they
    could represent gives the same amino acid (``GGN`` is glycine) and to
    ``unknown`` otherwise. A trailing partial codon is dropped.

    Args:
        seq: Nucleotide sequence; ``U`` and ``T`` are interchangeable.
        table: NCBI translation table id, name, or alias.
        cds: Treat ``seq`` as a complete coding sequence: require a whole
            number of codons that begins with an initiation codon of ``table``
            and ends with a stop, translate the initiator as ``M``, and omit
            the terminal stop.
        to_stop: Stop at the first stop codon and omit the stop symbol.
        stop_symbol: Symbol emitted for a stop codon.
        unknown: Symbol emitted for an undecidable codon.

    Raises:
        KeyError: If ``table`` is unknown.
        ValueError: If both ``cds`` and ``to_stop`` are set.
        InvalidCharacter: If ``seq`` contains an unknown character.
        CodingSequenceError: If ``cds`` is set and ``seq`` is not a valid CDS.
    """
    if cds and to_stop:
        raise ValueError("cds and to_stop cannot both be set")
    code = get_genetic_code(table)
    cleaned, _ = validate(seq)
    if cds:
        return _translate_cds(cleaned, code, unknown)
    residues: List[str] = []
    for start in range(0, len(cleaned) - len(cleaned) % 3, 3):
        amino_acid = _codon_amino_acid(cleaned[start : start + 3], code, unknown)
        if amino_acid == STOP:
            if to_stop:
                break
            residues.append(stop_symbol)
            continue
        residues.append(amino_acid)
    return "".join(residues)


def translate_frames(
    seq: str,
    frames: Iterable[int] = (1,),
    *,
    table: int | str = 1,
    cds: bool = False,
    to_stop: bool = False,
    stop_symbol: str = STOP,
    unknown: str = "X",
) -> List[FrameTranslation]:
    """Translate ``seq`` in each requested reading frame.

    Frames ``1``, ``2``, ``3`` start at the first, second, and third base of
    the given strand. Frames ``-1``, ``-2``, ``-3`` do the same on the reverse
    complement, so frame ``-1`` begins at the 3' end of the input.

    Args:
        seq: Nucleotide sequence over the IUPAC alphabet.
        frames: Frame numbers drawn from :data:`FRAMES`.
        table: NCBI translation table id, name, or alias.
        cds: Require each frame to be a complete coding sequence; see
            :func:`translate`.
        to_stop: Stop each frame at its first stop codon.
        stop_symbol: Symbol emitted for a stop codon.
        unknown: Symbol emitted for an undecidable codon.

    Returns:
        One :class:`FrameTranslation` per requested frame, in the order given.

    Raises:
        ValueError: If a frame number is not one of :data:`FRAMES`.
        InvalidCharacter: If ``seq`` contains an unknown character.
    """
    wanted = list(frames)
    bad = [frame for frame in wanted if frame not in FRAMES]
    if bad:
        raise ValueError(
            "frame must be one of "
            + ", ".join(str(frame) for frame in FRAMES)
            + f"; got {bad[0]}"
        )
    cleaned, _ = validate(seq)
    forward = cleaned
    backward = reverse_complement(cleaned) if any(f < 0 for f in wanted) else ""
    results: List[FrameTranslation] = []
    for frame in wanted:
        strand_seq = forward if frame > 0 else backward
        offset = abs(frame) - 1
        window = strand_seq[offset:]
        results.append(
            FrameTranslation(
                frame=frame,
                protein=translate(
                    window,
                    table=table,
                    cds=cds,
                    to_stop=to_stop,
                    stop_symbol=stop_symbol,
                    unknown=unknown,
                ),
                trailing_bases=len(window) % 3,
            )
        )
    return results
