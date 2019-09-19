"""GoldenGate assembly design process and steps."""

from typing import List

from Bio.Restriction import BsaI, BpiI
from Bio.Restriction.Restriction import RestrictionType
from Bio.SeqRecord import SeqRecord

from .clone import Clone
from ..assembly import goldengate
from ..containers import Container, Well
from ..designs import Design
from ..instructions import Temperature
from ..mix import Mix
from ..reagents import Reagent
from ..steps import Setup, Pipette, ThermoCycle, HeatShock


GOLDEN_GATE_MIX = Mix(
    {Reagent("master mix"): 4.0, SeqRecord: 2.0},
    fill_with=Reagent("water"),
    fill_to=20.0,
)


class GoldenGate(Clone):
    """GoldenGate assembly.

    Takes the design of a protocol and finds combinations of
    SeqRecords that will circularize into valid new plasmids.

    This Protocol is derived from NEB's guide:
    https://www.neb.com/protocols/2018/10/02/golden-gate-assembly-protocol-for-using-neb-golden-gate-assembly-mix-e1601

    Full responsibilities include

    * subselecting the input designs that will form valid plasmids after digestion with BsaI and BpiI
    * adding NEB Golden Gate Assembly Mix to the valid designs as Contents for a Container
    * computing the final plasmid (SeqRecord) after ligation of digested fragments
    * memoize the sorted enzyme + ids -> SeqRecord from step #3, pass as a
    * "mutate" method for the step after ThermoCycle()
    * add steps to carry out the rest of the assembly (heat shock, incubate, etc)

    Keyword Args:
        name: Name of this protocol
        enzymes: List of digest enzymes to catalyze SeqRecords with
        include: Include only plasmids with a feature matching something
            in the include list use in backbone selection (default: {None})
        min_count: The minimum number of SeqRecords in an assembly for it to
            be considered valid. smaller assemblies are ignored
        separate_reagents: Whether to separate reagent plate from other wells
    """

    def __init__(
        self,
        design: Design = Design(),
        name: str = "",
        enzymes: List[RestrictionType] = [BsaI, BpiI],
        include: List[str] = None,
        min_count: int = -1,
        separate_reagents: bool = False,
    ):
        super().__init__(name=name, design=design, separate_reagents=separate_reagents)

        self.min_count = min_count
        self.enzymes = enzymes
        self.mix = GOLDEN_GATE_MIX

    def run(self):
        """Filter designs to those that will form valid and new Golden Gate plasmids.

        Run each step and accumulate Golden Gate Assembly steps and layouts.
        """

        # get all the unique contents and set them up in their own wells
        mixed_wells = self._create_mixed_wells()

        for step in [
            Setup(
                target=mixed_wells,
                instructions=[
                    "Dilute plasmid DNA to 75 ng/µL in 'water'",
                    "Create 'assembly-mix' from 1:1 T4 Ligase Buffer (10X) and NEB Golden Gate Assembly Mix",
                ],
            ),
            Pipette(name="Mix DNA with assembly-mix and water", target=mixed_wells),
            ThermoCycle(
                [
                    Temperature(temp=37, time=3600),
                    Temperature(temp=4, time=-1),  # hold at 4 degrees
                ],
                mutate=self.mutate,  # set the SeqRecords
            ),
        ] + HeatShock:
            step(self)

    def _create_mixed_wells(self) -> List[Container]:
        """Return the valid circularizable assemblies.

        Also build up the dictionary for `self.wells_to_construct`, a map from sorted
        Fragment IDs to the SeqRecord that they will form after digestion and ligation.

        Returns:
            List[Container] -- list of wells to mix fragments for GoldenGate
        """

        mixed_wells: List[Container] = []
        for plasmids, fragments in goldengate(
            self.design,
            self.enzymes,
            include=self.include,
            min_count=self.min_count,
            linear=self.design.linear,
        ):
            # add reaction mix and water
            well_contents, well_volumes = self.mix(fragments)

            # create a well that mixes the assembly mix, plasmids, and reagents
            well = Well(contents=well_contents, volumes=well_volumes)

            # used in self.mutate
            self.wells_to_construct[well] = Well(plasmids, [sum(well_volumes)])
            mixed_wells.append(well)

        if not mixed_wells:
            raise RuntimeError(f"Failed to create any GoldenGate assemblies")

        return sorted(mixed_wells)
