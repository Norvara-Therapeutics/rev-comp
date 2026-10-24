"""Minimal FASTA reading and writing.

Permissive about blank lines and trailing whitespace; strict about sequence
text before the first header, the usual sign of a truncated paste.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import IO, Iterable, Iterator, List

__all__ = ["FastaError", "Record", "parse", "read_sequences", "format_record", "write"]

#: Default residues per line when writing FASTA.
DEFAULT_LINE_WIDTH = 60


class FastaError(ValueError):
    """Raised when FASTA input cannot be parsed."""


@dataclass
class Record:
    """One FASTA record.

    Attributes:
        name: The identifier, i.e. the header up to the first whitespace.
        description: Whatever followed the identifier on the header line.
        sequence: The residues, with all whitespace removed.
    """

    name: str
    sequence: str
    description: str = ""

    @property
    def header(self) -> str:
        """The header line without its leading ``>``."""
        if self.description:
            return f"{self.name} {self.description}"
        return self.name

    def with_suffix(self, suffix: str, sequence: str) -> "Record":
        """Return a copy carrying ``sequence`` and ``suffix`` in its header."""
        description = f"{self.description} {suffix}".strip() if suffix else self.description
        return Record(name=self.name, sequence=sequence, description=description)


def parse(source: IO[str] | str | Iterable[str]) -> Iterator[Record]:
    """Yield :class:`Record` objects from FASTA text, a file, or lines.

    Args:
        source: An open text handle, a string of FASTA text, or any iterable
            of lines.

    Yields:
        One record per header line, in file order.

    Raises:
        FastaError: If residues appear before the first header, or a header
            has an empty identifier.
    """
    lines = source.splitlines() if isinstance(source, str) else source
    name: str | None = None
    description = ""
    chunks: List[str] = []
    for lineno, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield Record(name=name, sequence="".join(chunks), description=description)
            header = line[1:].strip()
            if not header:
                raise FastaError(f"line {lineno}: FASTA header has no identifier")
            name, _, description = header.partition(" ")
            description = description.strip()
            chunks = []
        else:
            if name is None:
                excerpt = line if len(line) <= 20 else line[:20] + "..."
                raise FastaError(
                    f"line {lineno}: sequence data before any '>' header "
                    f"({excerpt!r})"
                )
            chunks.append(line)
    if name is not None:
        yield Record(name=name, sequence="".join(chunks), description=description)


def read_sequences(text: str, *, default_name: str = "sequence") -> List[Record]:
    """Read FASTA records, or a single bare sequence, from ``text``.

    Input that does not begin with ``>`` is treated as one unnamed sequence so
    that a sequence can be pasted or piped in directly.

    Args:
        text: FASTA text or raw residues.
        default_name: Record name given to a bare sequence.

    Raises:
        FastaError: If FASTA-looking input is malformed, or input is empty.
    """
    stripped = text.strip()
    if not stripped:
        raise FastaError("no input sequence")
    # A '>' anywhere at the start of a line means this was meant to be FASTA,
    # even if the first record's header is missing; parsing it gives a far
    # clearer diagnostic than treating '>' as an illegal residue.
    if any(line.startswith(">") for line in stripped.splitlines()):
        records = list(parse(stripped))
        if not records:
            raise FastaError("no FASTA records found")
        return records
    return [Record(name=default_name, sequence="".join(stripped.split()))]


def format_record(record: Record, *, line_width: int = DEFAULT_LINE_WIDTH) -> str:
    """Render one record as FASTA text ending in a newline.

    Args:
        record: The record to render.
        line_width: Residues per line; ``0`` writes one long line. An empty
            record is written as a bare header.
    """
    if not record.sequence:
        return ">" + record.header + "\n"
    if line_width > 0:
        body = [
            record.sequence[start : start + line_width]
            for start in range(0, len(record.sequence), line_width)
        ]
    else:
        body = [record.sequence]
    return ">" + record.header + "\n" + "\n".join(body) + "\n"


def write(
    records: Iterable[Record],
    handle: IO[str],
    *,
    line_width: int = DEFAULT_LINE_WIDTH,
) -> int:
    """Write ``records`` to ``handle`` as FASTA and return the count written."""
    written = 0
    for record in records:
        handle.write(format_record(record, line_width=line_width))
        written += 1
    return written
