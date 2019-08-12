"""Containers hold SeqRecords, Primers, Enzymes, etc."""

import math
import string
from typing import Dict, Iterable, List, Union, Optional
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

    Arguments:
        content {Content} -- contents of some container

    Returns:
        str -- the unique ID associated with the content

    Raises:
        TypeError: if unrecognized content type
    """

    if isinstance(content, SeqRecord):
        if content.id != "<unknown id>":
            return content.id
        return str(content.seq)
    elif isinstance(content, RestrictionType):
        return str(content)  # get enzyme cut seq
    elif isinstance(content, Primers):
        return "primers:" + str(content.fwd) + ";" + str(content.rev)
    elif isinstance(content, Reagent):
        return content.name
    elif isinstance(content, Species):
        return content.name
    raise TypeError(content)


class Container:
    """A container with contents.

    TODO: make contents a set

    Keyword Arguments:
        contents {List[Content]} -- the contents of the container (default: {None})
        volumes {List[float]} -- volumes of each content (default: {None})
    """

    volume_max = -1
    rows = 1
    cols = 1

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
        if len(self.volumes) != len(self.contents):
            raise ValueError

        self.withdrawn = 0.0  # volume with withdraws during pipette sim

    def add(self, contents: Union[Content, List[Content]]):
        """Add more content to this container

        Arguments:
            contents {Union[Content, List[Content]]} -- single content or list
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

        def content_rank(content: Content) -> str:
            """Unique id based on content type."""

            if isinstance(content, SeqRecord):
                return "0" + content_id(content)
            if isinstance(content, Reagent):
                return "1" + content_id(content)
            return "2" + content_id(content)

        def container_id(container: Container) -> str:
            return "".join(sorted(content_rank(c) for c in container))

        return container_id(self) < container_id(other)


class Well(Container):
    """A single well in a plate."""

    volume_max = 200
    rows = 8
    cols = 12


class Tube(Container):
    """A single tube for culturing or larger liquids.

    Based on Eppendorf 1.5 mL
    https://labware.opentrons.com/opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap
    """

    volume_max = 2000
    rows = 4
    cols = 6


class Reservoir(Container):
    """A single lab reservoir large volume liquids.

    Based on Agilent 1 Well Reservoir
    https://labware.opentrons.com/agilent_1_reservoir_290ml
    """

    volume_max = 290000
    rows = 1
    cols = 1


class Fridge(Container):
    """Ambiguous; something to retrieve from the fridge."""

    default_id = uuid.uuid1()  # default fridge is the same
    volume_max = 1_000_000_000  # like suitcase from harry potter

    def __init__(self, contents: Union[Content, List[Content]] = None):
        """Set the contents of the Well."""

        self.id: uuid.UUID = self.default_id

        if not isinstance(contents, list):
            self.contents = [contents if contents else []]
        else:
            self.contents = contents

    def empty(self, withdraw: float = 0) -> bool:
        """This isn't a LIMS system... Fridge is always full."""

        return False

    def withdraw(self, volume: float):
        """Do nothing, Fridge is infinite."""

        pass


class Layout:
    """Layout on a bench or a robotic platform

    Allocates positions to each grouped container (plate/reservoir, etc)

    Arguments:
        containers {Iterable[Container]} -- the list of containers to put in the plates

    Keyword Arguments:
        src_containers {Iterable[Container]} -- list of source containers (default: {None})
        existing_plates {int} -- the number of already existing plates
        log_volume {bool} -- whether to log volume to csv output (ex water(20)) (default: {False})
    """

    def __init__(
        self,
        containers: Iterable[Container],
        src_containers: Iterable[Container] = None,
        existing_plates: int = 0,
        log_volume: bool = False,
    ):
        # keep order consistent, see Container.__lt__ for sort method
        self.containers = sorted(containers)
        self.existing_plates = existing_plates
        self.log_volume = log_volume

        def get_containers(containers: Optional[Iterable[Container]], ctype: type):
            return [c for c in (containers or []) if isinstance(c, ctype)]

        # get reservoirs
        src_reservoir = get_containers(src_containers, Reservoir)
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
            src_well_plate = src_wells[0].rows * src_wells[0].cols
            dest_shift += well_count * math.ceil(len(src_wells) / src_well_plate)
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
    ) -> "Layout":
        """Create a Layout from an instruction with transfers.

        Arguments:
            instruction {Instruction} -- the instruction with transfers
            src_containers {bool} -- whether to also map out source containers
            existing_plates {int} -- the number of already existing plates
            log_volume {bool} -- whether to log each wells volume during to_csv

        Returns:
            Layout -- new Layout for plates, reservoirs, etc
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
            )

        return cls(dest_wells, existing_plates=existing_plates, log_volume=log_volume)

    def to_csv(self) -> str:
        """Return self in CSV representation for a CSV/Excel file.

        If there are more wells than the max within a plate, multiple
        plates are returned within a row with the plate name to the top-left of each plate.
        Each cells holds the ids of all well containers.

        Returns:
            str -- CSV representation of plate row
        """

        if not self.wells:
            return ""

        well = self.wells[0]
        wells = well.rows * well.cols
        rows = string.ascii_uppercase[: well.rows]

        csv_wells: List[List[str]] = [[] for _ in range(well.rows + 1)]
        for i, container in enumerate(self.wells):
            if i % wells == 0:
                if i != 0:
                    for j in range(well.rows + 1):
                        csv_wells[j].append("")

                # add a plate name header
                csv_wells[0].append(
                    "Plate:" + str(math.floor(i / wells) + 1 + self.existing_plates)
                )

                for j in range(well.cols):  # add column headers
                    csv_wells[0].append(str(j + 1))

                for j in range(well.rows):  # add row headers
                    csv_wells[j + 1].append(rows[j])

            contents = " ".join(  # add contents to a well
                [
                    content_id(c)
                    if not self.log_volume
                    else content_id(c) + f"({round(container.volumes[k], 1)})"
                    for k, c in enumerate(container)
                ]
            )
            csv_wells[i % well.rows + 1].append(contents)

        plate_output = ""
        for row in csv_wells:
            plate_output += ",".join(row) + "\n"
        return plate_output + "\n"

    def _set_well_meta(self, containers: List[Well], shift: int):
        """Cache meta about a container: plate, well name and index."""

        if not containers:
            return

        well = containers[0]
        rows = well.rows
        cols = well.cols

        col_names = {n: str(n + 1) for n in range(cols + 1)}
        row_names = string.ascii_uppercase[: rows + 1]
        well_count = rows * cols

        for i, container in enumerate(containers):
            i += shift
            plate_index = math.floor(i / well_count)
            plate_name = "Plate:" + str(plate_index + 1 + self.existing_plates)

            well_index = i % well_count + 1  # 1-based well index: A2 == 2
            row_name = row_names[i % rows]
            col_name = col_names[math.floor(plate_index / rows)]
            well_name = row_name + col_name

            self.container_to_plate_name[container] = plate_name
            self.container_to_well_index[container] = well_index
            self.container_to_well_name[container] = well_name

    def _plate_count(
        self, reservoirs: List[Reservoir], tubes: List[Tube], wells: List[Well]
    ) -> int:
        """Return the total plate space for the reservoirs, tubes and wells.

        Each reservoir takes up one plate space.
        Tubes and wells take up a variable number of space, depending on
        the number or rows/columns of their respective rack/plate
        
        Arguments:
            reservoirs {List[Reservoir]} -- large liquid reservoirs
            tubes {List[Tube]} -- list of tubes for a rack
            wells {List[Well]} -- list of wells for a plate
        
        Returns:
            int -- the number of "plate" spaces this layout requires
        """

        count = 0
        count += len(reservoirs)

        if tubes:
            tubes_per_rack = tubes[0].rows * tubes[0].cols
            count += math.ceil(float(len(tubes)) / tubes_per_rack)

        if wells:
            wells_per_plate = wells[0].rows * wells[0].cols
            count += math.ceil(float(len(wells)) / wells_per_plate)

        return count

    def __len__(self) -> int:
        """Return the number of plate spaces this object requires on a robotic deck."""

        return self._plate_count(self.reservoirs, self.tubes, self.wells)
