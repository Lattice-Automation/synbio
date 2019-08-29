"""Gibson Assembly Composite Step."""

from typing import Dict, List, Tuple

from Bio.SeqRecord import SeqRecord

from ..assembly import gibson
from ..containers import Container, Well
from ..designs import Design
from ..instructions import Temperature
from ..mix import Mix
from ..primers import Primers
from ..protocol import Protocol
from ..reagents import Reagent
from ..steps import Setup, Pipette, ThermoCycle, Incubate, HeatShock


PCR_MIX = Mix(
    {
        SeqRecord: 1.0,
        Primers: 1.0,
        Reagent("dNTPs"): 0.5,
        Reagent("10X Standard Taq Buffer"): 2.5,
        Reagent("Taq Polymerase"): 0.125,
    },
    fill_with=Reagent("water"),
    fill_to=25,
)

GIBSON_MIX = Mix(
    {Reagent("Gibson master mix 2X"): 10.0, SeqRecord: 1.0},
    fill_with=Reagent("water"),
    fill_to=20.0,
)


class Gibson(Protocol):
    """Gibson Assembly.

    Create primer pairs that will anneal each SeqRecord to its two neighbors
    so that, after PCR, the fragments will circularize via a Gibson Assembly.

    Based on the Gibson Assembly protocol outlined on the NEB site for kit (e5510):
    https://www.neb.com/protocols/2012/12/11/gibson-assembly-protocol-e5510

    Keyword Args:
        hifi: whether to use NEB's HiFi assembly method
        gibson_mix: the assembly mix to use when mixing the Gibson wells. Based on NEB's (e5510)
    """

    def __init__(
        self,
        design: Design = Design(),
        name: str = "",
        hifi: bool = False,
        pcr_mix: Mix = PCR_MIX,
        gibson_mix: Mix = GIBSON_MIX,
        separate_reagents: bool = False,
    ):
        super().__init__(name=name, design=design, separate_reagents=separate_reagents)

        self.hifi = hifi
        self.pcr_mix = pcr_mix
        self.gibson_mix = gibson_mix

        # for mapping input gibson wells to their output contents
        self.gibson_product: Dict[Container, Container] = {}

    def run(self) -> "Protocol":
        """Run each step of the protcol. Build up the output records and instructions."""

        pcr_wells, gibson_wells = self._create_mixed_wells()

        # extend for 1 minute per kb: https://www.neb.com/protocols/0001/01/01/taq-dna-polymerase-with-standard-taq-buffer-m0273
        max_length_dna = 0
        for gibson_well in gibson_wells:
            for content in [c for c in gibson_well if isinstance(c, SeqRecord)]:
                max_length_dna = max(max_length_dna, len(content))
        extension_time = (float(max_length_dna) / 1000) * 60

        for step in [
            Setup(target=pcr_wells),
            Pipette(name="Mix DNA, assembly-mix and water", target=pcr_wells),
            ThermoCycle(
                name="PCR fragments to add junctions between fragments",
                temps=[
                    Temperature(temp=95, time=30),
                    Temperature(temp=60, time=60),
                    Temperature(temp=68, time=extension_time),  # hold at 4 degrees
                ],
                cycles=30,
            ),
            Pipette(
                name="Mix PCR'ed fragments together with Gibson master mix and water",
                target=gibson_wells,
            ),
            Incubate(
                Temperature(temp=50, time=900), mutate=self.mutate
            ),  # set the SeqRecords
        ] + HeatShock:
            step(self)

        return self

    def _create_mixed_wells(self) -> Tuple[List[Container], List[Container]]:
        """Create two sets of wells. One for PCR'ing fragments, another for mixing them.

        Returns:
            (List[Container], List[Container]) -- two lists. First is of PCR
                wells for PCR'ing the fragments in preparation for Gibson,
                the second is for mixing the PCR'ed fragments together for
                circularization
        """

        pcr_wells: List[Container] = []
        gibson_wells: List[Container] = []
        for records in self.design:
            plasmid, primer_pairs = gibson(records, hifi=self.hifi)

            for i, primers in enumerate(primer_pairs):
                # create a well that mixes the primers with the input fragment
                pcr_contents, pcr_volumes = self.pcr_mix([records[i], primers])
                pcr_well = Well(contents=pcr_contents, volumes=pcr_volumes)
                pcr_wells.append(pcr_well)

            gibson_contents, gibson_volumes = self.gibson_mix(records)
            gibson_well = Well(contents=gibson_contents, volumes=gibson_volumes)
            gibson_wells.append(gibson_well)

            # map input well to it's output contents after Gibson Assembly
            self.gibson_product[gibson_well] = Well(
                plasmid, volumes=[sum(gibson_volumes)]
            )

        if not pcr_wells:
            raise RuntimeError(f"Failed to create any Gibson assemblies.")

        return sorted(pcr_wells), sorted(gibson_wells)

    def mutate(self, well: Container) -> Container:
        """Given the contents of a well, return single SeqRecord after Gibson Assembly."""

        if well in self.gibson_product:
            return self.gibson_product[well]

        raise RuntimeError(f"Gibson assembly of an unknown container {well}")
