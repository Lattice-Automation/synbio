"""Example of Gibson cloning fragments together via PCR."""

import os

# 'synbio' uses the enzymes provided in Bio.Restriction
from Bio.SeqIO import parse
from synbio.assembly import gibson

insert = next(parse(os.path.join("..", "data", "gibson", "pdsred2.fa"), "fasta"))
backbone = next(parse(os.path.join("..", "data", "gibson", "pDusk.fa"), "fasta"))

# assembled the SeqRecords together, creating primers for each SeqRecord so that
# after PCR they will anneal to one another
fragments = [insert, backbone]
plasmid, primer_pairs = gibson(fragments)

# inspect the plasmids that will be created
for fragment, primers in zip(fragments, primer_pairs):
    print("PCR", fragment.id, "with FWD", primers.fwd, "and REV", primers.rev)
