"""Output of a python script for the OpenTrons platform."""

from typing import Dict, List, Union

from ..containers import Plates
from ..instructions import Instruction

TRANSFER_MAX_VOLUME = 10.0
"""Max transfer of Echo is 10,000 nL or 10 uL.
https://www.labcyte.com/media/pdf/SPC-Echo-55X-Liquid-Handler.pdf
"""

TRANSFER_MULTIPLE = 0.0025
"""Transfers have to be a multiple of 2.5 nL.
https://github.com/Edinburgh-Genome-Foundry/Plateo/blob/master/plateo/exporters/picklist_to_labcyte_echo_picklist_file.py
"""


def to_labcyte(instruction: Instruction, existing_plates: int) -> str:
    """Given an instruction from a Protocol (from some Step), create a Labcyte picklist.

    Transfers have to be in multiples of 2.5 nL and cannot be more than 10 uL:
    https://www.labcyte.com/media/pdf/SPC-Echo-55X-Liquid-Handler.pdf

    There is really great documentation for this (in a relative sense) at:
    https://www.labcyte.com/documentation/ECP_HTML5/Content/PROJECTS/ECP_UG/CONTENT/c_CreatingPickLists.htm

    Arguments:
        instruction {Instruction} -- a single step's instruction
        existing_plates {int} -- number of plates before these in protocol

    Returns:
        str -- the picklist in CSV string format
    """

    if not instruction.transfers:
        raise ValueError(f"no transfers in Instruction: {instruction}")

    columns = [
        "Source Plate Barcode",
        "Source Well",
        "Destination Plate Barcode",
        "Destination Well",
        "Transfer Volume",
    ]

    # set each wells location in setup and destination plates
    plates = Plates.from_instruction(
        instruction, src_containers=True, existing_plates=existing_plates
    )

    rows: List[Dict[str, Union[str, float]]] = []
    for transfer in instruction.transfers:
        for split_transfer in transfer.split(TRANSFER_MAX_VOLUME, TRANSFER_MULTIPLE):
            src = split_transfer.src
            dest = split_transfer.dest

            rows += [
                {
                    "Source Plate Barcode": plates.container_to_plate_name[src],
                    "Source Well": plates.container_to_well_name[src],
                    "Destination Plate Barcode": plates.container_to_plate_name[dest],
                    "Destination Well": plates.container_to_well_name[dest],
                    "Transfer Volume": split_transfer.volume * 1000,  # UGLY: uL to nL
                }
            ]

    def csv_row(row: Dict[str, Union[str, float]]) -> str:
        """Make each row ',' separated."""

        return ",".join(str(row[c]) for c in columns)

    header_row = ",".join(columns) + "\n"
    return header_row + "\n".join(csv_row(row) for row in rows).strip()
