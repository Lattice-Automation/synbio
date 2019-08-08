"""Example of a CombinatorialBins MoClo assembly with steps and output.

All possible assemblies from combinations of all SeqRecords across
multiple bins are tested and those that will successfully linearize
are combined into composite plasmids.
"""

import os

from Bio.SeqIO import parse

from synbio import Protocol, CombinatorialBins
from synbio.composite import GoldenGate


def read_all_records():
    """Gather all SeqRecords from "goldengate" dir in examples."""

    records = []
    for file in os.listdir(os.path.join(".", "goldengate")):
        if file.endswith(".gb"):
            records.extend(parse(os.path.join(".", "goldengate", file), "genbank"))
    return records


# create a combinatorial library design from multiple "bins"
design = CombinatorialBins()
gg_records = read_all_records()
for type in ["promoter", "RBS", "CDS", "terminator"]:
    record_bin = [r for r in gg_records if any(f.type == type for f in r.features)]
    design.append(record_bin)  # add a new cominatorial bin

# create a protocol using GoldenGate as the sole composite step and run
protocol = Protocol(name="CombinatorialBins Golden Gate", design=design)

# filter on plasmids with a KanR feature in the backbone and at least 5 consituent SeqRecords
protocol.add(GoldenGate(resistance="KanR", min_count=5))
protocol.run()

# export all the output plasmids to a multi-FASTA
protocol.to_fasta("plasmids.fasta")

# export plate layouts
protocol.to_csv("plates.csv")

# export protocol
protocol.to_txt("protocol.txt")

# because this is a larger assembly, make robotic instructions
protocol.to_picklists("robotic_picklist.gwl", platform="tecan")
