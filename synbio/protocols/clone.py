"""Clone fragments together with digestion and ligation."""

from statistics import mean
from typing import Dict, List

from Bio.Restriction.Restriction import RestrictionType
from Bio.SeqRecord import SeqRecord

from ..assembly import clone_many_combinatorial
from ..containers import Container, Well
from ..designs import Design
from ..instructions import Temperature
from ..mix import Mix
from ..protocol import Protocol
from ..reagents import Reagent
from ..steps import Setup, Pipette, ThermoCycle, HeatShock


CLONING_MIX = Mix(
    {Reagent("10X NEBuffer"): 5.0, SeqRecord: 1.0, RestrictionType: 1.0},
    fill_with=Reagent("water"),
    fill_to=50.0,
)


class Clone(Protocol):
    """Clone SeqRecords together using BioPython enzymes.

    Digest the SeqRecords with all the Enzymes provided, find valid circularized
    assemblies, and create a protocol for preparing and ligating the fragments.

    This protocol is based on NEB's clone guide:
    https://www.neb.com/tools-and-resources/usage-guidelines/clone-guide

    Keyword Args:
        enzymes: the list of enzymes to digest the SeqRecords with
        include: include only plasmids with a feature matching something
            in the include list use in backbone selection (default: {None})
        mix: the assembly mix to use when mixing the assemblies with enzymes
        min_count: the minimum number of SeqRecords in an assembly for it to
            be considered valid. smaller assemblies are ignored
    """

    def __init__(
        self,
        name: str = "",
        design: Design = Design(),
        enzymes: List[RestrictionType] = None,
        mix: Mix = CLONING_MIX,
        include: List[str] = None,
        min_count: int = -1,
        separate_reagents: bool = False,
    ):
        super().__init__(name=name, design=design, separate_reagents=separate_reagents)

        self.enzymes = enzymes or []
        self.include = include
        self.mix = mix
        self.min_count = min_count
        self.wells_to_construct: Dict[Container, Container] = {}

    def run(self):
        """Filter designs to those that will form valid and new plasmids.

        Run each Clone step on the protocol. See:
        https://www.neb.com/protocols/2018/10/02/golden-gate-assembly-protocol-for-using-neb-golden-gate-assembly-mix-e1601
        """

        # get all the unique contents and set them up in their own wells
        mixed_wells = self._create_mixed_wells()

        # get the mean incubation temperature from the enzymes
        incubate_temp = mean([e.opt_temp for e in self.enzymes])

        for step in [
            Setup(target=mixed_wells),
            Pipette(
                name="Mix DNA with the enzymes, NEBuffer, and water", target=mixed_wells
            ),
            ThermoCycle(
                [
                    Temperature(temp=incubate_temp, time=3600),
                    Temperature(temp=65, time=300),
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
            a list of wells to mix fragments for Clone
        """

        if not self.enzymes:
            raise ValueError("Clone protocol lacks list of BioPython Enzymes")

        mixed_wells: List[Container] = []
        for plasmids, fragments in clone_many_combinatorial(
            self.design,
            enzymes=self.enzymes,
            include=self.include,
            min_count=self.min_count,
            linear=self.design.linear,
        ):
            # add reaction mix and water
            well_contents, well_volumes = self.mix(fragments + self.enzymes)

            # create a well that mixes the assembly mix, plasmids, and reagents
            well = Well(contents=well_contents, volumes=well_volumes)

            # used in self.mutate
            self.wells_to_construct[well] = Well(plasmids, [sum(well_volumes)])
            mixed_wells.append(well)

        if not mixed_wells:
            raise RuntimeError(f"Failed to create any Clone assemblies")

        return sorted(mixed_wells)

    def mutate(self, well: Container) -> Container:
        """Given the contents of a well, return single SeqRecord after digest/ligation."""

        if well in self.wells_to_construct:
            return self.wells_to_construct[well]

        raise KeyError(f"{well} not recognized as a Clone assembly")
