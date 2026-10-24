# rev-comp

Reverse complement, transcription and multi-frame translation from the command
line. Pure Python with no third-party dependencies, reading FASTA or a bare
sequence and writing FASTA.

## Install

Install the prerelease from PyPI:

```bash
python -m pip install --pre norvara-rev-comp
```

The Python package and command keep the shorter project name:

```python
import revcomp
```

```bash
rev-comp --version
```

## Quickstart

```bash
rev-comp translate --cds --input examples/phix174_geneJ.fasta
```

```
>NC_001422.1:848-964_f1 Escherichia phage phiX174 gene J, complete cds frame=+1 table=1
MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF
```

```bash
rev-comp complement -i examples/phix174.fasta -o revcomp_phix174.fasta
rev-comp translate --frame all -t 11 -i examples/phix174_geneJ.fasta
rev-comp complement ATGGCCATTGTAATG        # or pipe a sequence on stdin
```

`examples/phix174.fasta`, shipped with the repository and the source
distribution, is the complete phiX174 genome (RefSeq NC_001422.1);
`examples/phix174_geneJ.fasta` is its annotated gene J, nucleotides 848..964.

## Options

```
-i, --input PATH         FASTA or raw sequence to read (default: stdin)
-o, --output PATH        file to write (default: stdout)
-w, --line-width N       residues per output line; 0 for one line (default: 60)
    --name NAME          record name for a bare input sequence
    --permissive         replace characters outside the IUPAC alphabet with N
-q, --quiet              suppress warnings on standard error

complement - reverse complement, or complement or reverse alone
-m, --mode {revcomp,complement,reverse}
    --rna                write U in place of T

transcribe - DNA to RNA from either strand, or RNA back to DNA
-s, --strand {coding,template}
    --back               convert RNA to DNA instead

translate - one frame, six frames, or a checked coding sequence
-f, --frame FRAMES       1, 2, 3, -1, -2, -3, a comma-separated list, or
                         'forward', 'reverse', 'all' (default: 1)
-t, --table CODE         genetic code by number, name or alias (default: 1)
    --cds                require a complete coding sequence
    --to-stop            end each frame at its first stop codon
    --stop-symbol CHAR   symbol written for a stop (default: *)
    --unknown CHAR       symbol written for an undecidable codon (default: X)
```

`rev-comp tables` lists the 27 bundled genetic codes; `rev-comp tables --table N`
prints one in full.

## Reading the output

- Names and descriptions are carried through with the operation appended, so
  `>primer_F forward primer` becomes `>primer_F forward primer
  reverse-complement`. Translated records are named for their frame, `name_f1`
  to `name_r3`, so six-frame output has unique identifiers.
- Frames `1`, `2`, `3` are offsets into the sequence as given; `-1`, `-2`, `-3`
  are offsets into its reverse complement, so frame `-1` begins at the 3' end.
- `transcribe` assumes the **coding** strand and substitutes U for T. For a
  **template** strand pass `--strand template`, which reverse complements first:
  the transcript is antiparallel to the template, so plain substitution would
  give the right bases in the wrong order.
- `--cds` requires a whole number of codons from an initiation codon of the
  chosen table to a terminal stop, translates the initiator as methionine, omits
  the stop, and reports a sequence that is not a valid CDS rather than
  translating it anyway.
- Ambiguity codes are complemented to their IUPAC counterparts (R↔Y, K↔M, B↔V,
  D↔H) and translated when decidable (`GGN` is glycine) and to `X` otherwise.
  Case is preserved, gaps (`-`, `.`, `~`) pass through, and T and U are
  interchangeable, though a sequence holding both is warned about.
- Characters outside the IUPAC alphabet are reported by position and record and
  nothing is written; `--permissive` substitutes N and warns instead. Exit
  codes: `0` success, `1` I/O or conflicting arguments, `2` usage error, `3`
  input read but not valid sequence.

## Genetic codes

All 27 translation tables are bundled, transcribed verbatim from the `ncbieaa`
and `sncbieaa` strings of
[gc.prt](https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt), whose sha256 is
recorded in `codons.py`. Each table carries its own stop codons and its own
initiation codons, so `--cds` means something different under each one:

`--table` takes a number, an NCBI name, or an alias such as `standard`,
`bacterial` or `vertebrate-mito`. Tables 1 and 11 assign the same amino acids
and differ only in which codons may initiate, so `--cds` accepts a GTG or TTG
start under table 11 and rejects it under table 1.

## Limitations

- Input is held in memory, which suits oligos, genes, plasmids and viral genomes
  rather than vertebrate chromosomes. Validation completes before anything is
  written, so an invalid record leaves no partial output.
- Sequences are treated as linear and no coordinate subsetting is offered, so
  slice before piping in. Alignment gaps are accepted but not removed, and
  translate to `X` where they fall.
- Translation is frame arithmetic, not gene finding: no splice sites, ribosome
  binding sites or coding potential are considered.

## Development

```bash
pip install -e ".[dev]"
pytest -q
python -m pyflakes src/revcomp tests
```

## License

Licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0).
See [LICENSE](LICENSE) for the complete terms.
