"""Steps: lab processes to assemble a design."""

from collections import defaultdict
from typing import Callable, List, Optional, Dict, Union, Sequence

from .containers import Container, content_id, Content, Fridge
from .instructions import Temperature, Transfer, Instruction
from .protocol import Protocol, Step


class Setup(Step):
    """Create a list of transfers to make a list of setup containers

    All transfers for a Setup step come from the magical and bottomless 'Fridge'

    Keyword Arguments:
        target {List[Container]} -- target container list/ordering (default: {None})
        dest {Container} -- the type of target container (default: {first target container})
        name {str} -- the name of this step in the protocol
        instructions {List[str]} -- extra instructions to add to this step
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

    def execute(self, protocol: Protocol):
        """Create setup containers with enough contents to fill the target containers."""

        # accumulate the list of transfer volumes necessary
        id_to_content: Dict[str, Content] = {}
        id_to_volumes: Dict[str, List[float]] = defaultdict(list)
        for container in self.target:
            for i, content in enumerate(container):
                id_to_content[content_id(content)] = content
                id_to_volumes[content_id(content)].append(container.volumes[i])

        # create the setup containers and transfers based on volume/count
        container = self.dest or self.target[0]
        volume_max = container.volume_max
        setup_containers: List[Container] = []
        for cid, content in id_to_content.items():
            volumes = id_to_volumes[cid]

            new_container = container.create(content, volumes=[volumes.pop()])
            while volumes:
                volume = volumes.pop()

                if new_container.volume() + volume > volume_max:
                    setup_containers.append(new_container)
                    new_container = container.create(content, volumes=[volume])

                new_container.volumes[0] += volume
            setup_containers.append(new_container)

        transfers = [
            Transfer(src=Fridge(c), dest=c, volume=c.volume()) for c in setup_containers
        ]

        protocol.add_instruction(
            Instruction(
                name=self.name, transfers=transfers, instructions=self.instructions
            )
        )
        protocol.containers = setup_containers


class Pipette(Step):
    """Pipette the input containers to match a 'target' list of containers

    Create transfers to map the input containers to a new list of output containers (target)

    Keyword Arguments:
        target {List[Container]} -- target container list/ordering (default: {None})
        name {str} -- the name of this step in the protocol
        instructions {List[str]} -- extra instructions to add to this step
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

    def execute(self, protocol: Protocol):
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

                self.transfers.append(transfer)

        # update the protocol
        protocol.add_instruction(
            Instruction(
                name=self.name, transfers=self.transfers, instructions=self.instructions
            )
        )
        protocol.containers = list(self.target)


class Move(Step):
    """Move a fixed volume from each container to another container of same or new type

    Arguments:
        volume {float} -- the amount in millileters to move to the new container
        type {type} -- the type of new container to move the contents to
    """

    def __init__(self, volume: float, dest: Container = None, name: str = ""):
        super().__init__()

        self.name = name
        self.volume = volume
        self.dest = dest

    def execute(self, protocol: Protocol):
        """Create the pipette instructions for the move."""

        new_containers: List[Container] = []
        for container in protocol.containers:
            if self.dest:
                new_container = self.dest.create(container.contents)
            else:
                new_container = container.create(container.contents)
            new_containers.append(new_container)
            transfer = Transfer(src=container, dest=new_container, volume=self.volume)
            self.transfers.append(transfer)

        protocol.add_instruction(Instruction(name=self.name, transfers=self.transfers))
        protocol.containers = new_containers


class Add(Step):
    """Add contents to the existing containers.

    Arguments:
        add {Union[Container, Content]} -- the src container or content to add
        volume {float} -- the volume of the new content to add

    Keyword Arguments:
        name {str} -- the name of this step in the protocol
        instructions {List[str]} -- extra instructions to add to this step
    """

    def __init__(
        self,
        add: Union[Container, Content],
        volume: float,
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        if isinstance(add, Container):
            self.add = add
        else:
            self.add = Fridge(add)
        self.volume = volume
        self.name = name
        self.instructions = instructions if instructions else []

    def execute(self, protocol: Protocol):
        """Add new content from a src container to each other container."""

        for container in protocol.containers:
            transfer = Transfer(src=self.add, dest=container, volume=self.volume)
            self.transfers.append(transfer)

        protocol.add_instruction(
            Instruction(
                name=self.name, transfers=self.transfers, instructions=self.instructions
            )
        )


class ThermoCycle(Step):
    """Thermo cycle for PCR, digestion/ligation, etc.

    Make a list of temperature instructions with the temperature
    to set it at, the length of time and the number of cycles.

    For example, a ThermoCycle for a typical PCR would be expressed as:

    ```
    ThermoCycle(cycles=30, temps=[
        Temperature(temp=97, time=5 * 60),  # denature
        Temperature(temp=55, time=30),  # annealing
        Temperature(temp=72, time=60),  # extension
    ])
    ```

    Arguments:
        temps {List[Temperature]} -- list of temperatures in a gradient

    Keyword Arguments:
        name {str} -- the name of this step in the protocol
        cycles {int} -- the number of thermo cycles (default: {1})
        mutate {Optional[Callable[Container], List[Content]]]} --
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
        instructions {List[str]} -- list of additional instructions to add
    """

    def __init__(
        self,
        temps: List[Temperature],
        name: str = "",
        cycles: int = 1,
        mutate: Optional[Callable[[Container], Container]] = None,
        instructions: List[str] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = temps
        self.cycles = cycles
        self.mutate = mutate
        self.instructions = instructions if instructions else []

    def execute(self, protocol: Protocol):
        """Create an instruction for temps and their durations.

        If a mutate function was passed, mutate the contents
        of the containers appropriately.

        Arguments:
            protocol {Protocol} -- the protocol to add a thermoCycle step to
        """

        if self.mutate:
            new_containers: List[Container] = []
            for container in protocol.containers:
                new_containers.append(self.mutate(container))
            protocol.containers = new_containers

        if self.cycles > 1:
            self.instructions += [f"For {self.cycles} cycles"]

        protocol.add_instruction(
            Instruction(
                name=self.name, temps=self.temps, instructions=self.instructions
            )
        )


class Incubate(Step):
    """Incubate contents for some time.

    Incubate the contents of containers and in a fridge or some other incubator.

    Keyword Arguments:
        name {str} -- the name of this step in the protocol
        temp {Temperature} -- the target the container should be in (default: {None})
        mutate {Optional[Callable[Container], List[Content]]]} --
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
    """

    def __init__(
        self,
        temp: Temperature,
        name: str = "",
        mutate: Optional[Callable[[Container], Container]] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = [temp]
        self.mutate = mutate

    def execute(self, protocol: Protocol):
        """Create an instruction for temperatures and their durations.

        Arguments:
            protocol {Protocol} -- the protocol to add an incubtion step to
        """

        if self.mutate:
            new_containers: List[Container] = []
            for container in protocol.containers:
                new_containers.append(self.mutate(container))
            protocol.containers = new_containers

        protocol.add_instruction(Instruction(name=self.name, temps=self.temps))

