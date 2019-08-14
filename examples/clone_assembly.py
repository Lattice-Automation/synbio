"""Example of cloning together SeqRecords directly."""

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
# keep only the assembled plasmids that have both the DsRed2 and KanR features
assemblies = clone([insert, backbone], [NotI, BamHI], include=["DsRed2", "KanR"])

# inspect the plasmids that will be created
for plasmids, records in assemblies:
    for p in plasmids:
        print(p.id, len(p), "bp was cloned from:", ", ".join(r.id for r in records))
