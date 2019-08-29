"""Script generation for the opentrons platform."""

from enum import Enum

from ..instructions import Instruction


class Robot(Enum):
    """Opentrons has an OT-2 and an older OT-One.

    https://opentrons.com/
    """

    OT1 = "OT-One"
    OT2 = "OT-2"


def to_opentrons(
    instruction: Instruction,
    existing_plates: int,
    robot: Robot = Robot.OT2,
    separate_reagents: bool = False,
) -> str:
    """Given an instruction from a Protocol (from some Step), create an Opentrons script.

    Depends on the user having the 'opentrons' package installed

    API reference: https://docs.opentrons.com/api.html

    Example protocols: https://protocols.opentrons.com/protocol/dinosaur

    Args:
        instruction: a single step's instruction
        existing_plates: number of plates before these in protocol

    Keyword Args:
        robot: which of the Opentrons robots to make the script for
        separate_reagents: Whether to separate reagent plate from other wells

    Returns:
        the python script
    """

    pass
