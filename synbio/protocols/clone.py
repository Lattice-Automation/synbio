"""Clone fragments together with digestion and ligation."""

from statistics import mean
from typing import Dict, List

from Bio.Restriction.Restriction import RestrictionType
from Bio.SeqRecord import SeqRecord

from ..assembly import subclone
from ..containers import Container, Well
from ..instructions import Temperature
from ..mix import Mix
from ..protocol import Protocol
from ..reagents import Reagent
from ..species import Species
from ..steps import Setup, Pipette, Add, ThermoCycle, Incubate, Move


class Clone(Protocol):
    """Clone SeqRecords together using the enzymes provided.

    Digest the SeqRecords with all the Enzymes provided, find valid circularized
    assemblies, and create a protocol for preparing and ligating the fragments.

    This protocol is based on NEB's cloning guide:
    https://www.neb.com/tools-and-resources/usage-guidelines/cloning-guide

    Keyword Arguments:
        resistance {str} -- resistance to use in backbone selection (default: {""}),
            remove all circularizable assemblies missing resistance to this backbone
        mix {Mix} -- the assembly mix to use when mixing the assemblies with enzymes
        min_count {int} -- the minimum number of SeqRecords in an assembly for it to
            be considered valid. smaller assemblies are ignored
    """

    def __init__(
        self,
        *args,
        enzymes: List[RestrictionType] = None,
        mix: Mix = Mix(
            {Reagent("10X NEBuffer"): 5.0, SeqRecord: 1.0, RestrictionType: 1.0},
            fill_with=Reagent("water"),
            fill_to=50.0,
        ),
        resistance: str = "",
        min_count: int = -1,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.enzymes = enzymes or []
        self.resistance = resistance
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
            Setup(name="Setup PCR plate with (volumes) shown", target=mixed_wells),
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

    def _create_mixed_wells(self) -> List[Container]:
        """Return the valid circularizable assemblies.

        Also build up the dictionary for `self.wells_to_construct`, a map from sorted
        Fragment IDs to the SeqRecord that they will form after digestion and ligation.

        Returns:
            List[Container] -- list of wells to mix fragments for Clone
        """

        mixed_wells: List[Container] = []
        for plasmid, fragments in subclone(
            self.design,
            enzymes=[],
            resistance=self.resistance,
            min_count=self.min_count,
        ):
            # add reaction mix and water
            well_contents, well_volumes = self.mix(fragments + self.enzymes)

            # create a well that mixes the assembly mix, plasmids, and reagents
            well = Well(contents=well_contents, volumes=well_volumes)

            # used in self.mutate
            self.wells_to_construct[well] = Well(plasmid, [sum(well_volumes)])
            mixed_wells.append(well)

        if not mixed_wells:
            raise RuntimeError(f"Failed to create any Clone assemblies")

        return sorted(mixed_wells)

    def mutate(self, well: Container) -> Container:
        """Given the contents of a well, return single SeqRecord after digest/ligation."""

        if well in self.wells_to_construct:
            return self.wells_to_construct[well]

        raise KeyError(f"{well} not recognized as a Clone assembly")
