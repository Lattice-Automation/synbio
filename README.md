# synbio

`synbio` is a library for designing and assembling DNA. Users can design plasmids or libraries and export multi-step build protocols. Input SeqRecords. Output assembly SeqRecords, protocols, plate maps, and robotic picklists.

Documentation is available at [https://lattice-automation.github.io/synbio](https://lattice-automation.github.io/synbio)

## Installation

```bash
pip install synbio
```

## Models

`synbio` only expects the user to define their `Design` and `Protocol`. Several protocols are pre-defined.

### Designs

All are in `synbio.designs`:

- `Combinatorial` - list of SeqRecords to combinatorially anneal into all valid assemblies
- `CombinatorialBins` - list of bins of SeqRecords for combinatorial assembly of records between bins
- `Plasmid` - single list of SeqRecords to combine into a plasmid
- `PlasmidLibrary` - list of list of SeqRecords to combine into plasmids

### Protocols

All are in `synbio.protocols`:

- `Gibson` - Gibson assembly based on [NEB's E5510](https://www.neb.com/protocols/2012/12/11/gibson-assembly-protocol-e5510)
- `GoldenGate` - Golden Gate assembly based on [NEB's E1601](https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix#Product%20Information)

## Example

In the example below, the user specifies a combinatorial library design. All SeqRecords are tested for circularization with other SeqRecords. New and valid plasmids are assembled.

Behind the scenes, `synbio` is filtering all combinations of SeqRecords from the design that will circularize into valid plasmids (via [circuits in a graph](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0544-x/figures/1)). After running the `protocol`, users can export plate maps (`to_csv()`), composite plasmids (`to_fasta()`, `to_genbank()`), and assembly instructions (`to_txt()`, `to_picklists()`).

```python
"""Example of a Combinatorial Golden Gate assembly with human and robot output protocols."""

import os
from Bio.SeqIO import parse
from synbio.designs import Combinatorial
from synbio.protocols import GoldenGate

def read_all_records():
    gg_dir = os.path.join(".", "data", "goldengate")
    records = []
    for file in os.listdir(gg_dir):
        if file.endswith(".gb"):
            records.extend(parse(os.path.join(gg_dir, file), "genbank"))
    return records

# create a combinatorial library design from all valid combinations
design = Combinatorial(read_all_records())

# create a protocol using Golden Gate as the sole composite step and run
protocol = GoldenGate(
    name="Combinatorial Golden Gate",
    design=design,
    include=["KanR"],  # only keep circularized plasmids with a KanR SeqFeature
    min_count=5,  # only keep circularized plasmids from >=5 SeqRecords
)
protocol.to_fasta("plasmids.fasta")  # export multi-FASTA
protocol.to_csv("plates.csv")  # export plate layouts
protocol.to_txt("protocol.txt")  # export human protocol
protocol.to_picklists("picklist", platform="hamilton")  # export a hamilton picklist
```

plasmids.fasta:

```txt
>J23100_AB+B0032m_BC+C0012m_CD+B0015_DE+DVK_AE
GGAGTTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGCTACTAGAGTCACACAGGAAAG
TACTAAATGATGGTGAATGTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGT
...
```

plates.csv:

```csv
Setup Wells with volumes (uL) shown:
Plate:1,1,2,3,4,5,6,7,8,9,10,11,12
A,B0015_DE(4),C0080_CD(18),R0010_AB(54),water(36)
B,B0015_DE(160),DVK_AE(160),cre_CD(18),water(156)
...
```

protocol.txt:

```txt
Combinatorial GoldenGate:
1. Setup PCR plate with (volumes) shown:
	1.1. Dilute plasmid DNA to 75 ng/µL in 'water'
	1.2. Create 'assembly-mix' from 1:1 T4 Ligase Buffer (10X) and NEB Golden Gate Assembly Mix
...
```

picklist.gwl:

```txt
A;Plate:2;;;15;;2.0;;;
D;Plate:3;;;80;;2.0;;;
W;;;;;;;;;
...
```

## Annotation

In addition to DNA assembly, `synbio` exposes a plasmid annotation function in `synbio.features`. An example below shows a SeqRecord being augmented additional SeqFeatures from a curated database of common plasmid features.

```python
from Bio.SeqIO import parse
from synbio.features import annotate

record = next(parse("plasmid.fa", "fasta"))
record_with_features = annotate(record, identity=0.96)
```

## Alternatives

This is a non-exhaustive list. Contact me for a comparison of these libraries/platforms and `synbio`.

- [Aquarium](https://www.aquarium.bio/) is an extensive library/application for LIMS, protocol definition/execution, and workflow design. A lab operating system.
- [Autoprotocol](https://github.com/autoprotocol/autoprotocol-python) is a specification standard for experiments in the life sciences.
- [BioBricks](https://github.com/liaupm/BioBlocks) is a general focus, web-based editor for describing experiments in Biology.
- [Biocoder](https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-4-13) is a C++ library with extensive protocol step definition capabilities.
- [Plateo](https://github.com/Edinburgh-Genome-Foundry/Plateo) is a python library for planning, running and checking laboratory experiments. Great for parsing and exporting plates and picklists form multiple formats.
- [pydna](https://github.com/BjornFJohansson/pydna) is a python DNA assembly simulation library with a human-readable description of clone and assembly strategies.
