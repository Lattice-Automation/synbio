"""Containers hold SeqRecords, Primers, Enzymes, etc."""

import math
import string
from typing import Dict, Iterable, List, Optional, Union, Tuple
import uuid

from Bio.Restriction.Restriction import RestrictionType
from Bio.SeqRecord import SeqRecord

from .instructions import Instruction
from .primers import Primers
from .reagents import Reagent
from .species import Species

Content = Union[SeqRecord, RestrictionType, Primers, Reagent, Species]
"""The content of a container can be a sequence, enzyme, or primers."""


def content_id(content: Content) -> str:
    """Return a "unique" ID for each set of content.

    Args:
        content: The contents of some container

    Returns:
        The unique ID associated with the content

    Raises:
        TypeError: If unrecognized content type
    """

    if isinstance(content, SeqRecord):
        if content.id != "<unknown id>":
            return content.id
        return str(content.seq)
    if isinstance(content, RestrictionType):
        return str(content)  # get enzyme cut seq
    if isinstance(content, Primers):
        return "primers:" + str(content.fwd) + ";" + str(content.rev)
    if isinstance(content, Reagent):
        return content.name
    if isinstance(content, Species):
        return content.name
    raise TypeError(content)


class Container:
    """A container with contents.

    TODO: make contents a set for uniqueness

    Keyword Args:
        contents: the contents of the container (default: {None})
        volumes: volumes of each content (default: {None})
    """

    volume_dead = -1
    """Volume that should be left unused at bottom of container."""

    volume_max = -1
    """Max volume within each container."""

    rows = 1
    """Rows in the container (or its parents, as with Wells and their Plate)"""

    cols = 1
    """Cols in the container (or its parents, as with Wells and their Plate)"""

    def __init__(
        self,
        contents: Union[Content, List[Content]] = None,
        volumes: List[float] = None,
    ):
        self.id = uuid.uuid1()

        if not contents:
            self.contents: List[Content] = []
        elif not isinstance(contents, list):
            self.contents = [contents]
        else:
            self.contents = contents

        self.volumes = volumes if volumes else [-1] * len(self.contents)
        self.withdrawn = 0.0  # volume with withdraws during pipette sim

    def add(self, contents: Union[Content, List[Content]]):
        """Add more content to this container

        Args:
            contents: Single content or list
        """

        if isinstance(contents, list):
            self.contents.extend(contents)
        else:
            self.contents.append(contents)

    @classmethod
    def create(cls, contents: List[Content], **kwargs):
        """Copy this container type (cls) into a new instance with contents."""

        return cls(contents, **kwargs)

    def volume(self) -> float:
        """Return total volume of contents."""

        return sum(self.volumes)

    def withdraw(self, volume: float):
        """Simulate a volume withdraw from this container."""

        self.withdrawn += volume

        if self.withdrawn > self.volume():
            excess = self.withdrawn - self.volume()
            raise RuntimeError(f"withdrew {excess} uL too much from {self}")

    def empty(self, withdraw: float = 0) -> bool:
        """Whether there's enough volume left for a withdraw of volume size."""

        return self.withdrawn + withdraw > self.volume()

    def __hash__(self):
        """Return the unique id of this container."""

        return hash(self.id)

    def __len__(self):
        """Container's length is length of its content."""

        return len(self.contents)

    def __str__(self):
        return f"{type(self).__name__}:" + ",".join(content_id(c) for c in self)

    def __contains__(self, content: Content) -> bool:
        """Return whether the content is in this well."""

        return content_id(content) in {content_id(c) for c in self.contents}

    def __iter__(self):
        """Iterate over the contents of the container."""

        return iter(self.contents)

    def __getitem__(self, key) -> Content:
        """Return the content at the location passed."""

        return self.contents[key]

    def __lt__(self, other: "Container") -> bool:
        """Return whether this container should come before the other."""

        def min_rank_content(container: Container) -> int:
            min_rank = 1000
            for content in container:
                if isinstance(content, SeqRecord):
                    min_rank = min(0, min_rank)
                elif isinstance(content, Reagent):
                    min_rank = min(1, min_rank)
                elif isinstance(content, Species):
                    min_rank = min(2, min_rank)
                elif isinstance(content, RestrictionType):
                    min_rank = min(3, min_rank)
            return min_rank

        def container_id(container: Container) -> str:
            container_records = [c for c in container if isinstance(c, SeqRecord)]
            if container_records:
                return "".join(content_id(c) for c in container_records)
            return "".join(sorted(content_id(c) for c in container))

        if min_rank_content(self) < min_rank_content(other):
            return True

        return container_id(self) < container_id(other)


class Well(Container):
    """A single well in a plate."""

    volume_max = 200
    volume_dead = 15
    rows = 8
    cols = 12


class Tube(Container):
    """A single tube for culturing or larger liquids.

    Based on Eppendorf 1.5 mL
    https://labware.opentrons.com/opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap
    """

    volume_max = 2000
    volume_dead = 30
    rows = 4
    cols = 6


class Reservoir(Container):
    """A single lab reservoir large volume liquids.

    Based on Agilent 1 Well Reservoir
    https://labware.opentrons.com/agilent_1_reservoir_290ml
    """

    volume_max = 290000
    volume_dead = 1000
    rows = 1
    cols = 1


class Fridge(Container):
    """Ambiguous; a fridge in a lab."""

    volume_max = 1_000_000_000  # like suitcase from harry potter

    def __init__(self, contents: Union[Content, List[Content]] = None):
        """Set the contents of the Well."""

        super().__init__()

        self.id: uuid.UUID = uuid.uuid1()

        if not isinstance(contents, list):
            self.contents = [contents if contents else []]
        else:
            self.contents = contents

    def empty(self, withdraw: float = 0) -> bool:
        """This isn't a LIMS system... Fridge is always full."""

        return False

    def withdraw(self, volume: float):
        """Do nothing, Fridge is infinite."""


class Layout:
    """Layout on a bench or a robotic platform

    Allocates positions to each grouped container (plate/reservoir, etc)

    Attributes:
        containers: The list of containers to put in the plates
        src_containers: List of source containers (default: {None})
        existing_plates: The number of already existing plates
        log_volume: Whether to log volume to csv output (ex water(20)) (default: {False})
        separate_reagents: Whether to separate reagent plate from other wells
    """

    def __init__(
        self,
        containers: Iterable[Container],
        src_containers: Iterable[Container] = None,
        existing_plates: int = 0,
        log_volume: bool = False,
        separate_reagents: bool = False,
    ):
        self.containers = sorted(containers)
        self.existing_plates = existing_plates
        self.log_volume = log_volume
        self.separate_reagents = separate_reagents

        def get_containers(containers: Optional[Iterable[Container]], ctype: type):
            # keep order consistent, see Container.__lt__ for sort method
            return sorted([c for c in (containers or []) if isinstance(c, ctype)])

        # get reservoirs
        src_reservoirs = get_containers(src_containers, Reservoir)
        self.reservoirs = get_containers(self.containers, Reservoir)

        # get the tubes
        src_tubes = get_containers(src_containers, Tube)
        self.tubes = get_containers(self.containers, Tube)

        # only the wells go in plates
        src_wells = get_containers(src_containers, Well)
        self.wells = get_containers(self.containers, Well)

        # pre-computed maps from container to plate name, well index, well name
        self.container_to_plate_name: Dict[Container, str] = {}
        self.container_to_well_index: Dict[Container, int] = {}
        self.container_to_well_name: Dict[Container, str] = {}

        # if we're also adding source plates, we have to shift well index upwards
        if not self.wells:
            return

        well = self.wells[0]
        well_count = well.rows * well.cols
        dest_shift = 0
        if src_wells:
            dest_shift = (
                self._plate_count(self.reservoirs, self.tubes, self.wells) * well_count
            )
            src_wells = sorted(src_wells)
            self._set_well_meta(src_wells, 0)
        self._set_well_meta(self.wells, dest_shift)

    @classmethod
    def from_instruction(
        cls,
        instruction: Instruction,
        src_containers: bool = False,
        existing_plates: int = 0,
        log_volume: bool = False,
        separate_reagents: bool = False,
    ) -> "Layout":
        """Create a Layout from an instruction with transfers.

        Args:
            instruction: The instruction with transfers

        Keyword Args:
            src_containers: Whether to also map out source containers
            existing_plates: The number of already existing plates
            log_volume: Whether to log each wells volume during to_csv
            separate_reagents: Whether to separate reagent plate from other wells

        Returns:
            A new Layout for plates, reservoirs, etc
        """

        if not instruction.transfers:
            raise ValueError(f"instruction lacks transfers: {instruction}")

        dest_wells = {t.dest for t in instruction.transfers}
        if src_containers:
            return cls(
                dest_wells,
                src_containers={t.src for t in instruction.transfers},
                existing_plates=existing_plates,
                log_volume=log_volume,
                separate_reagents=separate_reagents,
            )

        return cls(
            dest_wells,
            existing_plates=existing_plates,
            log_volume=log_volume,
            separate_reagents=separate_reagents,
        )

    def to_csv(self) -> str:
        """Return self in CSV representation for a CSV/Excel file.

        If there are more wells than the max within a plate, multiple
        plates are returned within a row with the plate name to the top-left of each plate.
        Each cells holds the ids of all well containers.

        Returns:
            A CSV representation of row of containers, usually wells in plates
        """

        if not self.wells:
            return ""

        # initalize the cells with plates for the wells
        cells = self._wells_to_cells()

        if self.reservoirs:
            for i, reservoir in enumerate(self.reservoirs):
                res_name = "Reservoir:" + str(i + 1)
                res_volume = sum(reservoir.volumes)
                res_contents = "|".join(content_id(c) for c in reservoir)
                cells[0] += ["", f"{res_name}", f"{res_contents}({res_volume})"]

        plate_output = ""
        for row in cells:
            plate_output += ",".join(row) + "\n"
        return plate_output + "\n"

    def _set_well_meta(self, containers: List[Well], shift: int):
        """Save meta about a container: plate, well name and index."""

        if not containers:
            return

        well = containers[0]
        rows = well.rows
        cols = well.cols

        col_names = {n: str(n + 1) for n in range(cols + 1)}
        row_names = string.ascii_uppercase[: rows + 1]
        wells_per_plate = rows * cols

        def set_well_meta(wells: List[Well], well_shift: int):
            for i, container in enumerate(wells):
                i += well_shift
                plate_index = math.floor(i / wells_per_plate)
                plate_name = "Plate:" + str(plate_index + 1 + self.existing_plates)

                well_index = i % wells_per_plate + 1  # 1-based well index: A2 == 2
                row_name = row_names[i % rows]
                col_name = col_names[math.floor(plate_index / rows)]
                well_name = row_name + col_name

                self.container_to_plate_name[container] = plate_name
                self.container_to_well_index[container] = well_index
                self.container_to_well_name[container] = well_name

        if self.separate_reagents:
            wells_reagents, wells_other = self._separate_reagents(containers)
            set_well_meta(wells_other, shift)
            shift += (
                math.ceil(float(len(wells_other)) / wells_per_plate) * wells_per_plate
            )
            set_well_meta(wells_reagents, shift)
        else:
            set_well_meta(containers, shift)

    def _wells_to_cells(self) -> List[List[str]]:
        """Convert a list of wells to a list of list of strings for each well

        Returns:
            A list of list of strings, each a cell in CSV worksheet
        """

        well = self.wells[0]
        well_count = well.rows * well.cols
        rows = string.ascii_uppercase[: well.rows]

        cells: List[List[str]] = [[] for _ in range(well.rows + 1)]

        def add_wells(wells: List[Well]):
            for i, container in enumerate(wells):
                if i % well_count == 0:
                    if cells[0]:  # there are already other plates
                        for j in range(well.rows + 1):
                            cells[j].append("")

                    # add a plate name header
                    cells[0].append(self.container_to_plate_name[container])

                    # add column headers
                    for j in range(well.cols):
                        cells[0].append(str(j + 1))

                    # add row headers
                    for j in range(well.rows):
                        cells[j + 1].append(rows[j])

                contents = "|".join(  # add contents to a well
                    [
                        content_id(c)
                        if not self.log_volume
                        else content_id(c) + f"({round(container.volumes[k], 1)})"
                        for k, c in enumerate(container)
                    ]
                )

                cells[i % well.rows + 1].append(contents)

            # append a bunch of commas for the empty wells
            cells_remaining = well_count - len(wells) % well_count
            row = len(wells) % well.rows  # the last row
            for _ in range(cells_remaining):
                row += 1
                cells[row % well.rows + 1].append("")

        if self.separate_reagents:
            wells_reagents, wells_others = self._separate_reagents(self.wells)
            add_wells(wells_others)
            add_wells(wells_reagents)
        else:
            add_wells(self.wells)

        if not self.wells:
            cells = [[]]

        return cells

    def __len__(self) -> int:
        """Return the number of plate spaces this object requires on a robotic deck."""

        return self._plate_count(self.reservoirs, self.tubes, self.wells)

    def _plate_count(
        self, reservoirs: List[Reservoir], tubes: List[Tube], wells: List[Well]
    ) -> int:
        """Return the total plate space for the reservoirs, tubes and wells.

        Each reservoir takes up one plate space.
        Tubes and wells take up a variable number of space, depending on
        the number or rows/columns of their respective rack/plate

        Args:
            reservoirs: Large liquid reservoirs
            tubes: Tubes for a rack
            wells: Wells for a plate

        Returns:
            the number of "plate" spaces this layout requires
        """

        count = 0
        count += len(reservoirs)

        if tubes:
            tubes_per_rack = tubes[0].rows * tubes[0].cols
            count += math.ceil(float(len(tubes)) / tubes_per_rack)

        if wells:
            wells_per_plate = wells[0].rows * wells[0].cols
            if self.separate_reagents:
                wells_reagents, wells_other = self._separate_reagents(wells)
                count += math.ceil(float(len(wells_reagents)) / wells_per_plate)
                count += math.ceil(float(len(wells_other)) / wells_per_plate)
            else:
                count += math.ceil(float(len(wells)) / wells_per_plate)

        return count

    def _separate_reagents(
        self, containers: List[Container]
    ) -> Tuple[List[Container], List[Container]]:
        """Separate containers with only reagents from those with a mix of other things.
        
        Args:
            containers: The containers to separate by those with just reagents vs others
        
        Returns:
            A tuple with 1. containers with just reagents 2. containers with other contents
        """

        containers_reagents: List[Container] = []
        containers_other: List[Container] = []

        for container in containers:
            if all(isinstance(c, Reagent) for c in container):
                containers_reagents.append(container)
            else:
                containers_other.append(container)
        return containers_reagents, containers_other
