"""Example of a Gibson Assembly with two fragments."""

import os

from Bio.SeqIO import parse

from synbio.designs import Plasmid
from synbio.protocols import Gibson

# read in the SeqRecords from files
insert = next(parse(os.path.join("..", "data", "gibson", "BBa_K1085023.fa"), "fasta"))
backbone = next(parse(os.path.join("..", "data", "gibson", "pSB1C3.fa"), "fasta"))

# concatenate the insert to the backbone
design = Plasmid([insert, backbone])

# create a protocol using Gibson as the sole composite step and run
# filter on KanR in a backbone
protocol = Gibson(design=design)

# export the composite plasmid
protocol.to_fasta("plasmid.fasta")

# export plate layouts
protocol.to_csv("plates.csv")

# export human readable protocol
protocol.to_txt("protocol.txt")
