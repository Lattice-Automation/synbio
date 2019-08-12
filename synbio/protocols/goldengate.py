"""GoldenGate assembly design process and steps."""

from typing import Dict, List

from Bio.SeqRecord import SeqRecord

from ..assembly import goldengate
from ..containers import Container, content_id, Well
from ..instructions import Temperature
from ..mix import Mix
from ..protocol import Protocol
from ..reagents import Reagent
from ..species import Species
from ..steps import Setup, Pipette, Add, ThermoCycle, Incubate, Move


class GoldenGate(Protocol):
    """GoldenGate assembly.

    Takes the design of a protocol and finds combinations of
    SeqRecords that will circularize into valid new plasmids.

    Full responsibilities include:
        1. subselecting the input designs that will form valid plasmids
            after digestion with BsaI and BpiI
        2. adding NEB Golden Gate Assembly Mix to the valid designs as Contents for a Container
        3. computing the final plasmid (SeqRecord) after ligation of digested fragments
        4. memoize the sorted enzyme + ids -> SeqRecord from step #3, pass as a
            "mutate" method for the step after ThermoCycle()
        5. add steps to carry out the rest of the assembly (heat shock, incubate, etc)

    Keyword Arguments:
        resistance {str} -- resistance to use in backbone selection (default: {"KanR"}),
            remove all circularizable assemblies missing resistance to this backbone
        mix {Mix} -- the assembly mix to use when mixing the GoldenGate assemblies
        min_count {int} -- the minimum number of SeqRecords in an assembly for it to
            be considered valid. smaller assemblies are ignored
    """

    def __init__(
        self,
        *args,
        resistance: str = "",
        mix: Mix = Mix(
            {Reagent("master mix"): 4.0, SeqRecord: 2.0},
            fill_with=Reagent("water"),
            fill_to=20.0,
        ),
        min_count: int = -1,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        # search for this. must be in some SeqRecords features
        self.resistance = resistance
        self.mix = mix
        self.min_count = min_count
        self.id_to_well: Dict[str, Container] = {}

    def run(self):
        """Filter designs to those that will form valid and new GoldenGate devices.

        Run each GoldenGate step on the protocol. See:
        https://www.neb.com/protocols/2018/10/02/golden-gate-assembly-protocol-for-using-neb-golden-gate-assembly-mix-e1601
        """

        # get all the unique contents and set them up in their own wells
        mixed_wells = self._create_mixed_wells()

        for step in [
            Setup(
                name="Setup PCR plate with (volumes) shown",
                target=mixed_wells,
                instructions=[
                    "Dilute plasmid DNA to 75 ng/µL in 'water'",
                    "Create 'assembly-mix' from 1:1 T4 Ligase Buffer (10X) and NEB Golden Gate Assembly Mix",
                ],
            ),
            Pipette(
                name="Mix plasmids DNA (2 uL) 'assembly-mix' (4 uL) and 'water' (14 uL)",
                target=mixed_wells,
            ),
            ThermoCycle(
                [
                    Temperature(temp=37, time=3600),
                    Temperature(temp=60, time=300),
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

        Also build up the dictionary for `self.id_to_well`, a map from sorted
        Fragment IDs to the SeqRecord that they will form after digestion and ligation.

        Returns:
            List[Container] -- list of wells to mix fragments for GoldenGate
        """

        mixed_wells: List[Container] = []
        for assembly in goldengate(self.design, resistance=self.resistance):
            # add reaction mix and water
            well_contents, well_volumes = self.mix(assembly)

            # create a well that mixes the assembly mix, plasmids, and reagents
            well = Well(contents=well_contents, volumes=well_volumes)

            # used in self.mutate
            self.id_to_well[str(well.id)] = well
            mixed_wells.append(well)

        if not mixed_wells:
            raise RuntimeError(f"Failed to create any GoldenGate assemblies")

        return sorted(mixed_wells)

    def mutate(self, well: Container) -> Container:
        """Given the contents of a well, return single SeqRecord after digest/ligation."""

        well_id = str(well.id)
        if well_id not in self.id_to_well:
            raise KeyError(f"{well_id} not recognized as a GoldenGate assembly")

        # create hew new composite SeqRecord after digestion/ligation
        records = [c for c in well if isinstance(c, SeqRecord)]
        record_ids = [content_id(r) for r in records]
        record_combined = records[0].upper()
        for other_record in records[1:]:
            record_combined += other_record.upper()
        record_combined.id = "|".join(record_ids)

        # return a new well with the new combined SeqRecord
        return well.create([record_combined])
