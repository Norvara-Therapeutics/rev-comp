"""Command line interface for ``rev-comp``."""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
from typing import IO, List, Sequence

from revcomp import __version__
from revcomp.codons import available_tables, get_genetic_code, resolve_table_name
from revcomp.fasta import DEFAULT_LINE_WIDTH, FastaError, Record, read_sequences, write
from revcomp.sequence import (
    FRAMES,
    CodingSequenceError,
    InvalidCharacter,
    SequenceError,
    back_transcribe,
    complement,
    reverse,
    reverse_complement,
    transcribe,
    translate_frames,
    validate,
)

PROGRAM = "rev-comp"

#: Process exit codes.
EXIT_OK = 0
EXIT_ERROR = 1  # I/O or usage problem reported by this program
EXIT_SEQUENCE = 3  # input was read but is not valid sequence

_FRAME_GROUPS = {
    "all": list(FRAMES),
    "forward": [1, 2, 3],
    "reverse": [-1, -2, -3],
}


def _parse_frames(value: str) -> List[int]:
    """Parse the ``--frame`` argument into an ordered list of frame numbers."""
    tokens = [token.strip() for token in value.replace(" ", ",").split(",")]
    frames: List[int] = []
    for token in tokens:
        if not token:
            continue
        key = token.lower()
        if key in _FRAME_GROUPS:
            frames.extend(_FRAME_GROUPS[key])
            continue
        try:
            frame = int(token)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"{token!r} is not a frame; use 1, 2, 3, -1, -2, -3, "
                "'forward', 'reverse', or 'all'"
            ) from None
        if frame not in FRAMES:
            raise argparse.ArgumentTypeError(
                f"frame {frame} is out of range; use 1, 2, 3, -1, -2, or -3 "
                "(write negatives as --frame=-1 so they are not read as options)"
            )
        frames.append(frame)
    if not frames:
        raise argparse.ArgumentTypeError("no reading frame given")
    seen: List[int] = []
    for frame in frames:
        if frame not in seen:
            seen.append(frame)
    return seen


def _single_char(value: str) -> str:
    """Validate a residue symbol written for stops or undecidable codons."""
    if len(value) != 1:
        raise argparse.ArgumentTypeError(
            f"expected a single character, not {value!r}"
        )
    return value


def _line_width(value: str) -> int:
    """Validate the residues-per-line argument."""
    try:
        width = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not a number") from None
    if width < 0:
        raise argparse.ArgumentTypeError(
            "line width cannot be negative; use 0 for one line per record"
        )
    return width


def _table(value: str) -> int:
    """Resolve a ``--table`` argument to an NCBI translation table id."""
    try:
        return resolve_table_name(value)
    except KeyError:
        raise argparse.ArgumentTypeError(
            f"unknown genetic code {value!r}; run '{PROGRAM} tables' to list them"
        ) from None


def _common_io(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "sequence",
        nargs="?",
        help="sequence given directly on the command line; omit to read --input",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="PATH",
        help="FASTA file or raw sequence file to read (default: standard input)",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="file to write (default: standard output)",
    )
    parser.add_argument(
        "-w",
        "--line-width",
        type=_line_width,
        default=DEFAULT_LINE_WIDTH,
        metavar="N",
        help="residues per output line; 0 for one line per record "
        f"(default: {DEFAULT_LINE_WIDTH})",
    )
    parser.add_argument(
        "--name",
        default="sequence",
        metavar="NAME",
        help="record name to use when the input is a bare sequence "
        "(default: sequence)",
    )
    parser.add_argument(
        "--permissive",
        action="store_true",
        help="replace characters outside the IUPAC alphabet with N and warn, "
        "instead of failing",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="suppress warnings on standard error",
    )


def build_parser() -> argparse.ArgumentParser:
    """Construct the full argument parser."""
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description=(
            "Reverse complement, transcribe, and translate nucleotide "
            "sequences. Reads FASTA or a bare sequence, writes FASTA."
        ),
        epilog=(
            "Examples:\n"
            f"  {PROGRAM} complement --input primer.fasta --output rev_primer.fasta\n"
            f"  {PROGRAM} transcribe --input gene.fasta --output transcript.fasta\n"
            f"  {PROGRAM} translate --input gene.fasta --frame 1 --table standard\n"
            f"  {PROGRAM} translate --frame all --table 2 ATGGCCATTGTAATG\n"
            f"  echo ATGCGT | {PROGRAM} complement\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"{PROGRAM} {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    complement_parser = subparsers.add_parser(
        "complement",
        help="reverse complement a sequence (or just complement, or just reverse)",
        description=(
            "Reverse complement a nucleotide sequence: the opposite strand read "
            "5'->3'. Ambiguity codes are complemented to their IUPAC "
            "counterparts (R<->Y, K<->M, B<->V, D<->H) and case is preserved."
        ),
    )
    _common_io(complement_parser)
    complement_parser.add_argument(
        "-m",
        "--mode",
        choices=("revcomp", "complement", "reverse"),
        default="revcomp",
        help="revcomp: reverse and complement (default); complement: complement "
        "in place; reverse: reverse without complementing",
    )
    complement_parser.add_argument(
        "--rna",
        action="store_true",
        help="write the output with U in place of T",
    )
    complement_parser.set_defaults(handler=_run_complement)

    transcribe_parser = subparsers.add_parser(
        "transcribe",
        help="convert DNA to RNA (or back)",
        description=(
            "Transcribe DNA to RNA. By default the input is taken to be the "
            "coding (sense) strand, which is copied base for base with T "
            "replaced by U. Pass --strand template if the input is the "
            "template (antisense) strand, in which case it is reverse "
            "complemented first."
        ),
    )
    _common_io(transcribe_parser)
    transcribe_parser.add_argument(
        "-s",
        "--strand",
        choices=("coding", "template"),
        default=None,
        help="which strand the input represents (default: coding)",
    )
    transcribe_parser.add_argument(
        "--back",
        action="store_true",
        help="reverse the operation: convert RNA to DNA by replacing U with T",
    )
    transcribe_parser.set_defaults(handler=_run_transcribe)

    translate_parser = subparsers.add_parser(
        "translate",
        help="translate nucleotides to protein in one or more reading frames",
        description=(
            "Translate a nucleotide sequence. Codons containing ambiguity "
            "codes are translated when every base they could stand for gives "
            "the same amino acid (GGN is glycine) and to X otherwise. A "
            "trailing partial codon is dropped with a warning."
        ),
    )
    _common_io(translate_parser)
    translate_parser.add_argument(
        "-f",
        "--frame",
        type=_parse_frames,
        default=[1],
        metavar="FRAMES",
        help="reading frames: 1, 2, 3, -1, -2, -3, a comma-separated list, or "
        "'forward', 'reverse', 'all' for all six (default: 1). Write negative "
        "frames as --frame=-1",
    )
    translate_parser.add_argument(
        "-t",
        "--table",
        type=_table,
        default=1,
        metavar="CODE",
        help="NCBI genetic code as a number, name, or alias, e.g. 1, "
        "'standard', 2, 'vertebrate-mitochondrial', 11, 'bacterial' "
        "(default: 1)",
    )
    translate_parser.add_argument(
        "--cds",
        action="store_true",
        help="treat the input as a complete coding sequence: require a whole "
        "number of codons from an initiation codon of --table to a stop, "
        "translate the initiator as M, and omit the terminal stop",
    )
    translate_parser.add_argument(
        "--to-stop",
        action="store_true",
        help="end each frame at its first stop codon, omitting the stop symbol",
    )
    translate_parser.add_argument(
        "--stop-symbol",
        type=_single_char,
        default="*",
        metavar="CHAR",
        help="symbol written for a stop codon (default: *)",
    )
    translate_parser.add_argument(
        "--unknown",
        type=_single_char,
        default="X",
        metavar="CHAR",
        help="symbol written for a codon whose amino acid is undecidable "
        "(default: X)",
    )
    translate_parser.set_defaults(handler=_run_translate)

    tables_parser = subparsers.add_parser(
        "tables",
        help="list the available genetic codes",
        description="List the NCBI translation tables accepted by --table.",
    )
    tables_parser.add_argument(
        "-t",
        "--table",
        type=_table,
        metavar="CODE",
        help="show the full codon assignment for one table",
    )
    tables_parser.set_defaults(handler=_run_tables)

    return parser


class UsageError(Exception):
    """Raised for argument combinations argparse cannot express."""


def _read_input(args: argparse.Namespace, stdin: IO[str]) -> List[Record]:
    if args.sequence is not None and args.input:
        raise UsageError("give a sequence or --input, not both")
    if args.sequence is not None:
        text = args.sequence
    elif args.input and args.input != "-":
        with open(args.input, "r", encoding="utf-8") as handle:
            text = handle.read()
    else:
        text = stdin.read()
    return read_sequences(text, default_name=args.name)


def _clean(record: Record, args: argparse.Namespace, warn) -> Record:
    """Validate one record's residues, honouring --permissive."""
    sequence, replaced = validate(
        record.sequence,
        permissive=args.permissive,
        record=record.name,
    )
    if replaced:
        warn(
            f"{record.name}: replaced {len(replaced)} character(s) outside the "
            f"IUPAC alphabet with N (first at position {replaced[0]})"
        )
    upper = sequence.upper()
    if "T" in upper and "U" in upper:
        warn(f"{record.name}: contains both T and U; treating U as T")
    return Record(name=record.name, sequence=sequence, description=record.description)


def _run_complement(args: argparse.Namespace, records: List[Record], warn) -> List[Record]:
    out: List[Record] = []
    for record in records:
        if args.mode == "revcomp":
            sequence = reverse_complement(record.sequence, rna=args.rna)
            suffix = "reverse-complement"
        elif args.mode == "complement":
            sequence = complement(record.sequence, rna=args.rna)
            suffix = "complement"
        else:
            sequence = reverse(record.sequence)
            if args.rna:
                sequence = sequence.replace("T", "U").replace("t", "u")
            suffix = "reversed"
        out.append(record.with_suffix(suffix, sequence))
    return out


def _run_transcribe(args: argparse.Namespace, records: List[Record], warn) -> List[Record]:
    if args.back and args.strand is not None:
        raise UsageError("--back cannot be combined with --strand")
    args.strand = args.strand or "coding"
    out: List[Record] = []
    for record in records:
        if args.back:
            sequence = back_transcribe(record.sequence)
            suffix = "cDNA"
        else:
            sequence = transcribe(record.sequence, strand=args.strand)
            suffix = "mRNA" if args.strand == "coding" else "mRNA from template strand"
        out.append(record.with_suffix(suffix, sequence))
    return out


def _frame_tag(frame: int) -> str:
    """Short, unique identifier suffix for a reading frame: f1-f3, r1-r3."""
    return ("f" if frame > 0 else "r") + str(abs(frame))


def _run_translate(args: argparse.Namespace, records: List[Record], warn) -> List[Record]:
    if args.cds and args.to_stop:
        raise UsageError("--cds cannot be combined with --to-stop")
    if args.cds and len(args.frame) != 1:
        raise UsageError("--cds applies to a single reading frame")
    code = get_genetic_code(args.table)
    out: List[Record] = []
    for record in records:
        # Offset frames necessarily leave bases over, so the informative fact
        # is whether the input itself is a whole number of codons.
        leftover = len(record.sequence) % 3
        if leftover and not args.cds:
            warn(
                f"{record.name}: length {len(record.sequence)} is not a multiple "
                f"of 3; trailing bases are dropped from each frame"
            )
        try:
            translations = translate_frames(
                record.sequence,
                args.frame,
                table=args.table,
                cds=args.cds,
                to_stop=args.to_stop,
                stop_symbol=args.stop_symbol,
                unknown=args.unknown,
            )
        except CodingSequenceError as exc:
            raise CodingSequenceError(f"{record.name}: {exc}") from None
        for item in translations:
            out.append(
                Record(
                    name=f"{record.name}_{_frame_tag(item.frame)}",
                    sequence=item.protein,
                    description=(
                        f"{record.description} frame={item.frame:+d} "
                        f"table={code.id}"
                    ).strip(),
                )
            )
    return out


def _run_tables(args: argparse.Namespace, stdout: IO[str]) -> int:
    if args.table:
        code = get_genetic_code(args.table)
        stdout.write(f"NCBI table {code.id}: {'; '.join(code.names)}\n")
        stdout.write(f"start codons: {', '.join(sorted(code.start_codons))}\n")
        stdout.write(f"stop codons:  {', '.join(sorted(code.stop_codons))}\n\n")
        for first in "TCAG":
            for third in "TCAG":
                row = []
                for second in "TCAG":
                    codon = first + second + third
                    row.append(f"{codon} {code.codons[codon]}")
                stdout.write("   ".join(row) + "\n")
            stdout.write("\n")
        return EXIT_OK
    for code in available_tables():
        stdout.write(f"{code.id:>3}  {'; '.join(code.names)}\n")
    return EXIT_OK


def main(
    argv: Sequence[str] | None = None,
    *,
    stdin: IO[str] | None = None,
    stdout: IO[str] | None = None,
    stderr: IO[str] | None = None,
) -> int:
    """Run the command line interface and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    stdin = stdin if stdin is not None else sys.stdin
    stdout = stdout if stdout is not None else sys.stdout
    stderr = stderr if stderr is not None else sys.stderr

    if not getattr(args, "command", None):
        parser.print_help(stdout)
        return EXIT_OK

    quiet = getattr(args, "quiet", False)

    def warn(message: str) -> None:
        if not quiet:
            stderr.write(f"{PROGRAM}: warning: {message}\n")

    if args.command == "tables":
        return _run_tables(args, stdout)

    try:
        records = _read_input(args, stdin)
    except UsageError as exc:
        stderr.write(f"{PROGRAM}: error: {exc}\n")
        return EXIT_ERROR
    except FastaError as exc:
        stderr.write(f"{PROGRAM}: error: {exc}\n")
        return EXIT_SEQUENCE
    except OSError as exc:
        source = repr(args.input) if args.input else "standard input"
        stderr.write(f"{PROGRAM}: error: cannot read {source}: {exc.strerror}\n")
        return EXIT_ERROR

    try:
        cleaned = [_clean(record, args, warn) for record in records]
        results = args.handler(args, cleaned, warn)
    except UsageError as exc:
        stderr.write(f"{PROGRAM}: error: {exc}\n")
        return EXIT_ERROR
    except InvalidCharacter as exc:
        stderr.write(f"{PROGRAM}: error: {exc}\n")
        stderr.write(
            f"{PROGRAM}: pass --permissive to replace them with N and continue\n"
        )
        return EXIT_SEQUENCE
    except SequenceError as exc:
        stderr.write(f"{PROGRAM}: error: {exc}\n")
        return EXIT_SEQUENCE

    for record in results:
        if not record.sequence:
            warn(f"{record.name}: result is empty")

    try:
        if args.output and args.output != "-":
            with open(args.output, "w", encoding="utf-8") as handle:
                write(results, handle, line_width=args.line_width)
        else:
            write(results, stdout, line_width=args.line_width)
            stdout.flush()
    except BrokenPipeError:
        # Downstream closed the pipe, as `| head` does: a normal end, not an
        # error. Point the stream at the null device so that the interpreter's
        # final flush does not fail again on the way out.
        with contextlib.suppress(Exception):
            os.dup2(os.open(os.devnull, os.O_WRONLY), stdout.fileno())
        return EXIT_OK
    except OSError as exc:
        destination = repr(args.output) if args.output else "standard output"
        stderr.write(f"{PROGRAM}: error: cannot write {destination}: {exc.strerror}\n")
        return EXIT_ERROR
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
