"""Checks against an annotated reference.

Gene J of Escherichia phage phiX174 (RefSeq NC_001422.1, nucleotides 848..964)
and its recorded protein product, alongside the complete genome shipped in
``examples/``. Together they exercise the genetic code, the frame arithmetic and
the strand mapping on real sequence rather than on constructed inputs.
"""

from pathlib import Path

from revcomp.fasta import parse
from revcomp.sequence import reverse_complement, translate, translate_frames

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

GENE_J_PROTEIN = "MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF"
GENE_J_START, GENE_J_END = 848, 964


def read(name):
    """Read the single record held in an example file."""
    with open(EXAMPLES / name, encoding="utf-8") as handle:
        (record,) = parse(handle)
    return record


def test_the_genome_example_is_the_complete_record():
    genome = read("phix174.fasta")
    assert genome.name == "NC_001422.1"
    assert len(genome.sequence) == 5386


def test_the_gene_example_is_the_annotated_slice_of_the_genome():
    genome = read("phix174.fasta")
    gene = read("phix174_geneJ.fasta")
    assert gene.sequence == genome.sequence[GENE_J_START - 1 : GENE_J_END]
    assert len(gene.sequence) == 117


def test_cds_translation_reproduces_the_deposited_protein():
    gene = read("phix174_geneJ.fasta")
    assert translate(gene.sequence, cds=True) == GENE_J_PROTEIN
    assert translate(gene.sequence) == GENE_J_PROTEIN + "*"


def test_the_gene_is_recovered_from_the_reverse_complement():
    """Frame -1 of the reverse complement is frame +1 of the gene."""
    gene = read("phix174_geneJ.fasta")
    (minus,) = translate_frames(reverse_complement(gene.sequence), [-1])
    assert minus.protein == GENE_J_PROTEIN + "*"


def test_reverse_complement_round_trips_the_whole_genome():
    genome = read("phix174.fasta")
    assert reverse_complement(reverse_complement(genome.sequence)) == genome.sequence
