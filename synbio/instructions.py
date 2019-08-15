"""Instructions for biologists and robots.

Created by steps within a Protocol.
"""

import math
from typing import List
from uuid import uuid1


class Transfer:
    """Transfer contents from one container to another.

    Arguments:
        src {Container} -- the source Container
        dest {Container} -- the destination Container
        volume {float} -- the volume to transfer in microliters
    """

    def __init__(self, src, dest, volume: float = 10.0):
        self.src = src
        self.dest = dest
        self.volume = volume

    def split(self, max_volume: float, multiple_of: float) -> List["Transfer"]:
        """Split a Transfer into multiple other Transfers based on max volume.

        This was necessitated by the Labcyte Echo that only transfers up to 10 uL
        and with transfers that are a multiple of 2.5 nL.

        Arguments;
            max_volume {float} -- the max volume of a single transfer in uL
            multiple_of {float} -- each transfer has to be a multiple of this (in uL)

        Returns:
            List[Transfer] -- list of transfers meeting restraints
        """

        transfer_count = math.ceil(self.volume / max_volume)
        volume_per_transfer = self.volume / transfer_count
        volume_per_transfer = round(volume_per_transfer / multiple_of) * multiple_of

        split_transfers: List["Transfer"] = []
        for _ in range(transfer_count):
            split_transfers.append(Transfer(self.src, self.dest, volume_per_transfer))
        return split_transfers

    def __hash__(self):
        return hash(self.src) + hash(self.dest) + hash(self.volume)


class Temperature:
    """A temperature instruction, has a temperature and time component.

    Arguments:
        temp {float} -- temperature in degrees Celcius
        time {float} -- length of time in seconds
    """

    def __init__(self, temp: float = 27.0, time: float = 10.0):
        self.temp = temp
        self.time = time


class Instruction:
    """A single instruction set generated by a single Step."""

    def __init__(
        self,
        name: str = "",
        transfers: List[Transfer] = None,
        temps: List[Temperature] = None,
        instructions: List[str] = None,
    ):
        self.id = uuid1()
        self.name = name
        self.transfers = transfers
        self.temps = temps
        self.instructions = instructions or []

    def to_txt(self, index: int = -1) -> str:
        """Create a string representation of this instruction."""

        if self.temps:
            self.instructions += [_temp_instructions(self.temps)]

        if not self.name and not self.instructions:
            return ""

        if not self.name and len(self.instructions) > 1:
            self.name = self.instructions.pop(0)

        txt = ""
        if self.name and self.instructions:
            if index > 0:
                txt += f"{index}. {self.name}:\n"
                for j, sub_instruction in enumerate(self.instructions):
                    txt += f"\t{index}.{j+1}. {sub_instruction}\n"
            else:
                txt += f"{self.name}:\n"
                for j, sub_instruction in enumerate(self.instructions):
                    txt += f"\t{sub_instruction}\n"
        elif self.name:
            if index > 0:
                txt += f"{index}. {self.name}\n"
            else:
                txt += f"{self.name}\n"
        elif self.instructions:
            if index > 0:
                for sub_instruction in self.instructions:
                    txt += f"{index}. {sub_instruction}\n"
            else:
                for sub_instruction in self.instructions:
                    txt += f"{sub_instruction}\n"
        return txt

    def __hash__(self):
        """Create a unique hash for this Instruction."""

        return hash(self.id)


def to_txt(name: str, instructions: List[Instruction]) -> str:
    """Return a text representation of the instructions for a protocol.

    ```txt
    Combinatorial MoClo:
    1. Mix the Assembly mix
        1.1. Create 200 µL 'assembly-mix' from 1:1 T4 ligase buffer (10X) and NEB Golden Gate Assembly Mix
    2. Setup the Setup Plate as specified
    ```

    Arguments:
        name {str} -- the protocol's name
        instructions {str} -- the list of protocol instructions

    Returns:
        str -- a string representation of the protocol
    """

    txt = name + ":\n" if name else ""
    i = 1
    for instruction in instructions:
        instruction_txt = instruction.to_txt(i)
        if instruction_txt:
            txt += instruction_txt
            i += 1

    return txt


def _temp_instructions(temps: List[Temperature] = None) -> str:
    """Create an instruction for a thermocycle or incubate step."""

    temps = temps or []
    temps_str = [
        f"{t.temp} °C {_display_time(t.time)}" if t.time > 0 else f"Hold at {t.temp} °C"
        for t in temps
    ]
    return "Heat at " + ", ".join(temps_str)


def _display_time(seconds: float, granularity=2) -> str:
    """Return a human readable instruction for some time duration.

    Copy and pasted from:
    https://stackoverflow.com/questions/4048651/python-function-to-convert-seconds-into-minutes-hours-and-days/4048773

    If the number is positive, it returns the time with a granularity of 2.
    Example:

    ```python
    >>>display_time(3601)
    ..."for 1 hours 1 seconds"
    ```

    If the number is negative, hold the temperature.
    Example:

    ```python
    >>>display_time(-1)
    ..."hold"
    ```

    Arguments:
        seconds {int} -- the number of seconds of a step

    Keyword Arguments:
        granularity {int} -- the granularity of the displayed time

    Returns:
        str -- the steps total time in human readable string format
    """

    if seconds < 1:
        return "hold"

    seconds = round(seconds)
    intervals = (
        ("weeks", 604_800),  # 60 * 60 * 24 * 7
        ("days", 86400),  # 60 * 60 * 24
        ("hours", 3600),  # 60 * 60
        ("minutes", 60),
        ("seconds", 1),
    )
    result = []
    for name, count in intervals:
        value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip("s")
            result.append("{} {}".format(value, name))
    return "for " + ", ".join(result[:granularity])
