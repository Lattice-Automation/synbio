"""Example of annotating a Bio.SeqRecord with common plasmid Bio.SeqFeatures."""

from Bio.SeqIO import parse
from synbio.features import annotate

record = next(parse("plasmid.fa", "fasta"))
record_with_features = annotate(record, identity=0.96)
