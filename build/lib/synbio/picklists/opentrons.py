"""Output of a python script for the OpenTrons platform."""

from typing import List

from ..instructions import Instruction


def to_opentrons(instructions: List[Instruction]) -> str:
    """Given a list of instructions, generate a python Opentrons script

    API documentation for the Opentrons platform is available at:
    https://docs.opentrons.com/writing.html

    Use of this script requires the `opentrons` package that is not
    installed with this library. You can install it with `pip install opentrons`

    Arguments:
        instructions {List[Instruction]} -- list of assembly instructions

    Returns:
        str -- the file contents of the Opentrons automation script
    """

    pass
