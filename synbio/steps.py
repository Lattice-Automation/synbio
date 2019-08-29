"""Steps: lab processes to assemble a design."""

from collections import defaultdict
from typing import Callable, List, Optional, Dict, Sequence

from .containers import Container, content_id, Content, Fridge
from .instructions import Temperature, Transfer, Instruction
from .protocol import Protocol, Step
from .reagents import Reagent
from .species import Species


class Setup(Step):
    """Create a list of transfers to make a list of setup containers

    All transfers for a Setup step come from the magical and bottomless 'Fridge'

    Attributes:
        target: target container list/ordering (default: {None})
        dest: the type of target container (default: {first target container})
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step,
    """

    def __init__(
        self,
        target: Sequence[Container],
        dest: Container = None,
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        if not target:
            raise ValueError

        self.target = target
        self.dest = dest
        self.name = name
        self.instructions = instructions if instructions else []

    def __call__(self, protocol: Protocol):
        """Create setup containers with enough contents to fill the target containers."""

        # accumulate the list of transfer volumes necessary
        id_to_content: Dict[str, Content] = {}
        id_to_volumes: Dict[str, List[float]] = defaultdict(list)
        for container in self.target:
            for i, content in enumerate(container):
                cid = content_id(content)
                id_to_content[cid] = content
                id_to_volumes[cid].append(container.volumes[i])

        # create the setup containers and transfers based on volume/count
        container = self.dest or self.target[0]
        volume_max = container.volume_max
        setup: List[Container] = []
        for cid, content in id_to_content.items():
            volumes = id_to_volumes[cid]

            new_container = container.create(
                content, volumes=[volumes.pop() + container.volume_dead]
            )
            while volumes:
                volume = volumes.pop()

                if new_container.volume() + volume > volume_max:
                    setup.append(new_container)
                    new_container = container.create(content, volumes=[volume])

                new_container.volumes[0] += volume
            setup.append(new_container)

        # make a transfer to fill each setup container
        transfers = [Transfer(src=Fridge(c), dest=c, volume=c.volume()) for c in setup]
        dest_name = type(setup[0]).__name__
        instructions = [
            f"Setup each {dest_name} with volumes (uL) specified"
        ] + self.instructions

        instruction = Instruction(
            name=self.name, transfers=transfers, instructions=instructions
        )

        protocol.add_instruction(instruction)
        protocol.containers = sorted(setup)

        return protocol


class Pipette(Step):
    """Pipette the input containers to match a 'target' list of containers

    Create transfers to map the input containers to a new list of output containers (target)

    Attributes:
        target: target container list/ordering (default: {None})
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step
    """

    def __init__(
        self,
        target: Sequence[Container],
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        self.target = target
        self.name = name
        self.instructions = instructions if instructions else []

    def __call__(self, protocol: Protocol):
        """Transfer from protcol's current containers to match the target containers."""

        # create a map from current content to its source containers
        content_to_srcs: Dict[str, List[Container]] = defaultdict(list)
        for container in protocol.containers:
            for content in container:
                cid = content_id(content)
                content_to_srcs[cid].append(container)

        # for each target container, make a list of transfers to pipette
        # from the input container to output containers
        # create new containers at will here when others of source contents run out
        transfers: List[Transfer] = []
        for target_container in self.target:
            for i, content in enumerate(target_container):
                cid = content_id(content)

                if cid not in content_to_srcs:
                    content_to_srcs[cid] = [Fridge(content)]  # it's coming from fridge

                volume = target_container.volumes[i]  # volume of transfer

                # use the first container that isn't empty
                src = next(w for w in content_to_srcs[cid] if not w.empty(volume))
                src.withdraw(volume)

                transfer = Transfer(src=src, dest=target_container, volume=volume)
                transfers.append(transfer)

        # update the protocol
        protocol.add_instruction(
            Instruction(
                name=self.name, transfers=transfers, instructions=self.instructions
            )
        )
        protocol.containers = sorted(list(self.target))

        return protocol


class Move(Step):
    """Move a fixed volume from each container to another container of same or new type

    Attributes:
        volume: the amount in millileters to move to the new container
        type: the type of new container to move the contents to
    """

    def __init__(self, volume: float, dest: Container = None, name: str = ""):
        super().__init__()

        self.name = name
        self.volume = volume
        self.dest = dest

    def __call__(self, protocol: Protocol):
        """Create the pipette instructions for the move."""

        new_containers: List[Container] = []
        transfers: List[Transfer] = []
        for container in protocol.containers:
            new_container = (
                self.dest.create(container.contents)
                if self.dest
                else container.create(container.contents)
            )
            new_containers.append(new_container)

            transfer = Transfer(src=container, dest=new_container, volume=self.volume)
            transfers.append(transfer)

        src_name = type(protocol.containers[0]).__name__
        dest_name = type(new_containers[0]).__name__

        instructions = [
            f"Move {self.volume} uL from each {src_name} to new {dest_name}s"
        ]

        protocol.add_instruction(
            Instruction(name=self.name, transfers=transfers, instructions=instructions)
        )
        protocol.containers = sorted(new_containers)

        return protocol


class Add(Step):
    """Add contents to the existing containers.

    Attributes:
        add: the src container or content to add
        volume: the volume of the new content to add
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step
    """

    def __init__(
        self,
        add: Content,
        volume: float,
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        assert add, "Must select Content to add to each container in Add Step"

        if isinstance(add, Container):
            self.add = add
        else:
            self.add = Fridge(add)
        self.volume = volume
        self.name = name
        self.instructions = instructions if instructions else []

    def __call__(self, protocol: Protocol):
        """Add new content from a src container to each other container."""

        transfers: List[Transfer] = []
        for container in protocol.containers:
            transfer = Transfer(src=self.add, dest=container, volume=self.volume)
            transfers.append(transfer)

        container = protocol.containers[0]
        instructions = self.instructions or [
            f"Add {self.volume} uL of {content_id(self.add[0])} to each {type(container).__name__}"
        ]

        instruction = Instruction(
            name=self.name, transfers=transfers, instructions=instructions
        )
        protocol.add_instruction(instruction)

        return protocol


class ThermoCycle(Step):
    """Thermo cycle for PCR, digestion/ligation, etc.

    Make a list of temperature instructions with the temperature
    to set it at, the length of time and the number of cycles.

    Example:
        An example of a ThermoCycle for PCR:
        >>> ThermoCycle(cycles=30, temps=[
        ...    Temperature(temp=97, time=5 * 60),  # denature
        ...    Temperature(temp=55, time=30),  # annealing
        ...    Temperature(temp=72, time=60),  # extension
        ... ])

    Attributes:
        temps: list of temperatures in a gradient
        name: the name of this step in the protocol
        cycles: the number of thermo cycles (default: {1})
        mutate:
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
        instructions: list of additional instructions to add
        extension: last extension after other thermocycle steps.
            Example is a final 5 minute extension at the end of PCR that's common
    """

    def __init__(
        self,
        temps: List[Temperature],
        name: str = "",
        cycles: int = 1,
        mutate: Optional[Callable[[Container], Container]] = None,
        extension: Temperature = None,
        instructions: List[str] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = temps
        self.cycles = cycles
        self.mutate = mutate
        self.extension = extension
        self.instructions = instructions if instructions else []

    def __call__(self, protocol: Protocol):
        """Create an instruction for temps and their durations.

        If a mutate function was passed, mutate the contents
        of the containers appropriately.

        Args:
            protocol: the protocol to add a thermoCycle step to
        """

        if self.mutate:
            new_containers: List[Container] = []
            for container in protocol.containers:
                new_containers.append(self.mutate(container))
            protocol.containers = sorted(new_containers)

        if self.cycles > 1:
            self.instructions += [f"For {self.cycles} cycles:"]

        instruction = Instruction(
            name=self.name, temps=self.temps, instructions=self.instructions
        )
        protocol.add_instruction(instruction)

        return protocol


class Incubate(Step):
    """Incubate contents for some time at a set temperature.

    Incubate the contents of containers and in a fridge or some other incubator.

    Attributes:
        name: the name of this step in the protocol
        temp: the target the container should be in (default: {None})
        mutate:
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
    """

    def __init__(
        self,
        temp: Temperature,
        name: str = "Incubate",
        mutate: Optional[Callable[[Container], Container]] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = [temp]
        self.mutate = mutate

    def __call__(self, protocol: Protocol):
        """Create an instruction for temperatures and their durations.

        Args:
            protocol: the protocol to add an incubtion step to
        """

        if self.mutate:
            new_containers: List[Container] = []
            for container in protocol.containers:
                new_containers.append(self.mutate(container))
            protocol.containers = sorted(new_containers)

        protocol.add_instruction(Instruction(name=self.name, temps=self.temps))

        return protocol


HeatShock: List[Step] = [
    Move(volume=3.0),
    Add(add=Species("E coli"), volume=10.0),
    ThermoCycle(name="Heat shock", temps=[Temperature(temp=42, time=30)]),
    Add(add=Reagent("SOC"), volume=150.0),
    Incubate(temp=Temperature(temp=37, time=3600)),
]
"""A composite HeatShock step for getting DNA into E coli.

>>> [
...    Move(volume=3.0),
...    Add(add=Species("E coli"), volume=10.0),
...    ThermoCycle(name="Heat shock", temps=[Temperature(temp=42, time=30)]),
...    Add(add=Reagent("SOC"), volume=150.0),
...    Incubate(temp=Temperature(temp=37, time=3600)),
... ]
"""
