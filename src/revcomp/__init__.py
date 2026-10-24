"""rev-comp: reverse complement, transcribe, and translate nucleotide sequences.

The public API mirrors the command line:

    >>> from revcomp import reverse_complement, transcribe, translate
    >>> reverse_complement("ATGGCCATTGTAATG")
    'CATTACAATGGCCAT'
    >>> transcribe("ATGGCC")
    'AUGGCC'
    >>> translate("ATGGCCATTGTAATGGGC")
    'MAIVMG'
"""

from __future__ import annotations

__version__ = "1.0.0rc1"

from revcomp.codons import GeneticCode, available_tables, get_genetic_code
from revcomp.fasta import FastaError, Record, parse, read_sequences, write
from revcomp.sequence import (
    FRAMES,
    CodingSequenceError,
    FrameTranslation,
    InvalidCharacter,
    SequenceError,
    back_transcribe,
    complement,
    reverse,
    reverse_complement,
    transcribe,
    translate,
    translate_frames,
    validate,
)

__all__ = [
    "__version__",
    "FRAMES",
    "CodingSequenceError",
    "FastaError",
    "FrameTranslation",
    "GeneticCode",
    "InvalidCharacter",
    "Record",
    "SequenceError",
    "available_tables",
    "back_transcribe",
    "complement",
    "get_genetic_code",
    "parse",
    "read_sequences",
    "reverse",
    "reverse_complement",
    "transcribe",
    "translate",
    "translate_frames",
    "validate",
    "write",
]
