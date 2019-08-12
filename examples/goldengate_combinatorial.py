"""Example of a Combinatorial Golden Gate assembly with steps and output."""

import os

from Bio.SeqIO import parse

from synbio import Combinatorial
from synbio.protocols import GoldenGate


def read_all_records():
    """Gather all SeqRecords from "goldengate" dir in examples."""

    gg_dir = os.path.join(".", "goldengate")

    records = []
    for file in os.listdir(gg_dir):
        if file.endswith(".gb"):
            records.extend(parse(os.path.join(gg_dir, file), "genbank"))
    return records


# create a combinatorial library design from multiple "bins"
design = Combinatorial(read_all_records())

# create a protocol using Golden Gate as the sole composite step and run
protocol = GoldenGate(
    name="Combinatorial Golden Gate", design=design, resistance="KanR", min_count=5
)
protocol.run()

# export all the output plasmids to a multi-FASTA
protocol.to_fasta("composite_parts.fasta")

# export plate layouts
protocol.to_csv("plate_layouts.csv")

# export human protocol
protocol.to_txt("protocol.txt")

# export a hamilton picklist
protocol.to_picklists("robotic_picklist.gwl", platform="tecan")
