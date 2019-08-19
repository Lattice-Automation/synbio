"""Example of cloning together SeqRecords after digestion with enzymes."""

import os

# 'synbio' uses the enzymes provided in Bio.Restriction
from Bio.Restriction import NotI, BamHI
from Bio.SeqIO import parse
from synbio.assembly import clone

# insert has "DsRed2" between NotI and BamHI
insert = next(parse(os.path.join("..", "data", "cloning", "pdsred2.gb"), "genbank"))

# backbone has "KanR" resistance between BamHI and NotI
backbone = next(parse(os.path.join("..", "data", "cloning", "pdusk.gb"), "genbank"))

# simulate a digestion and ligation between the backbone and the insert via NotI and BamHI
# keep only the plasmids that have both the DsRed2 and KanR features
plasmids = clone([insert, backbone], [NotI, BamHI], include=["DsRed2", "KanR"])

# inspect the plasmids that will be created
for plasmid in plasmids:
    print(plasmid.id, f"{len(plasmid)}bp", plasmid.description)
