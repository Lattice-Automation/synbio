"""Gibson Assembly Composite Step."""

from typing import Dict, List, Tuple

from Bio.SeqRecord import SeqRecord

from ..assembly import gibson
from ..containers import Container, Well
from ..instructions import Temperature
from ..mix import Mix
from ..primers import Primers
from ..protocol import Protocol
from ..reagents import Reagent
from ..species import Species
from ..steps import Setup, Pipette, Add, ThermoCycle, Incubate, Move


class Gibson(Protocol):
    """Gibson Assembly.

    Based on the Gibson Assembly protocol outlined on the NEB site for kit (e5510):
    https://www.neb.com/protocols/2012/12/11/gibson-assembly-protocol-e5510

    Arguments:
        hifi {bool} -- whether to use NEB's HiFi assembly method

    Keyword Arguments:
        mix {Mix} -- the assembly mix to use when mixing the Gibson wells.
            Based on NEB's (e5510)
    """

    def __init__(
        self,
        *args,
        pcr_mix: Mix = Mix(
            {
                SeqRecord: 1.0,
                Primers: 1.0,
                Reagent("dNTPs"): 0.5,
                Reagent("10X Standard Taq Buffer"): 2.5,
                Reagent("Taq Polymerase"): 0.125,
            },
            fill_with=Reagent("water"),
            fill_to=25,
        ),
        gibson_mix: Mix = Mix(
            {Reagent("Gibson master mix 2X"): 10.0, SeqRecord: 1.0},
            fill_with=Reagent("water"),
            fill_to=20.0,
        ),
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.pcr_mix = pcr_mix
        self.gibson_mix = gibson_mix

        # for mapping input gibson wells to their output contents
        self.gibson_product: Dict[Container, Container] = {}

    def run(self) -> "Protocol":
        """Run each step of the protcol. Build up the output records and instructions."""

        pcr_wells, gibson_wells = self._create_mixed_wells()

        for step in [
            Setup(
                name="Setup PCR plate with fragments and primers at the (volumes) shown",
                target=pcr_wells,
            ),
            ThermoCycle(
                name="PCR fragments with added homology",
                temps=[
                    Temperature(temp=95, time=30),  # TODO: calculate these
                    Temperature(temp=60, time=60),
                    Temperature(temp=68, time=720),  # hold at 4 degrees
                ],
                cycles=30,
            ),
            Pipette(
                name="Mix PCR'ed fragments together for Gibson, along with Gibson master mix",
                target=gibson_wells,
            ),
            Incubate(
                Temperature(temp=50, time=900), mutate=self.mutate  # set the SeqRecords
            ),
            Move(name="Move 3 uL from each mixture well to new plate(s)", volume=3.0),
            Add(
                name="Add 10 uL of competent E. coli to each well",
                add=Species("competent_e_coli"),
                volume=10.0,
            ),
            ThermoCycle(name="Heat shock", temps=[Temperature(temp=42, time=30)]),
            Add(
                name="Add 150 uL of SOC media to each well",
                add=Reagent("soc_media"),
                volume=150.0,
            ),
            Incubate(name="Incubate", temp=Temperature(temp=37, time=3600)),
        ]:
            step.execute(self)

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
            plasmid, primer_pairs = gibson(records)

            for i, primers in enumerate(primer_pairs):
                # create a well that mixes the primers with the input fragment
                pcr_contents, pcr_volumes = self.pcr_mix([records[i], primers])
                pcr_well = Well(contents=pcr_contents, volumes=pcr_volumes)
                pcr_wells.append(pcr_well)

            gibson_contents, gibson_volumes = self.gibson_mix(records)
            gibson_well = Well(contents=gibson_contents, volumes=gibson_volumes)
            gibson_wells.append(gibson_well)

            # map input well to it's output contents after Gibson Assembly
            self.gibson_product[gibson_well] = Well(plasmid)

        if not pcr_wells:
            raise RuntimeError(f"Failed to create any Gibson assemblies")

        return sorted(pcr_wells), sorted(gibson_wells)

    def mutate(self, well: Container) -> Container:
        """Given the contents of a well, return single SeqRecord after digest/ligation."""

        if well in self.gibson_product:
            return self.gibson_product[well]

        raise RuntimeError(f"Gibson assembly of an unknown container {well}")
