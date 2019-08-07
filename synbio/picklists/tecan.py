"""Generate Tecan EVO picklists."""

from typing import Dict, List, Union

from ..containers import Plates
from ..instructions import Instruction


def to_tecan(instruction: Instruction, existing_plates: int) -> str:
    """Given a picklist for the Freedom EVO platform

    https://lifesciences.tecan.com/freedom-evo-platform

    This code is based entirely on Zulko's Plateo:
    https://github.com/Edinburgh-Genome-Foundry/Plateo/blob/master/plateo/exporters/picklist_to_tecan_evo_picklist_file.py

    Steps start with "A" for aspirate, "D" for dispense and "W" for a wash (new tip)
    between each pipette step

    Optimizations:
        TODO: reagents like water don't need "wash" steps between each dispense
        TODO: reagents like water can be done with one larger aspirate step
        TODO: reagents that ALWAYS mix don't need "wash" steps between dispenses

    Arguments:
        instruction {Instruction} -- a single step's instruction
        existing_plates {int} -- number of plates before these in protocol

    Returns:
        str -- the picklist in CSV string format
    """

    if not instruction.transfers:
        raise ValueError(f"no transfers in Instruction: {instruction}")

    columns = [
        "Action",
        "RackLabel",
        "RackID",
        "RackType",
        "Position",
        "TubeID",
        "Volume",
        "LiquidClass",
        "TipType",
        "TipMask",
    ]

    # set each wells location in setup and destination plates
    plates = Plates.from_instruction(
        instruction, src_containers=True, existing_plates=existing_plates
    )

    rows: List[Dict[str, Union[str, float]]] = []
    for transfer in instruction.transfers:
        src = transfer.src
        dest = transfer.dest
        volume = round(transfer.volume, 2)

        rows += [
            {
                "Action": "A",
                "RackLabel": plates.container_to_plate_name[src],
                "Position": plates.container_to_well_index[src],
                "Volume": volume,
            },
            {
                "Action": "D",
                "RackLabel": plates.container_to_plate_name[dest],
                "Position": plates.container_to_well_index[dest],
                "Volume": volume,
            },
            {"Action": "W"},
        ]

    def csv_row(row: Dict[str, Union[str, float]]) -> str:
        """Make each row ';' separated."""

        return ";".join(str(row.get(c, "")) for c in columns)

    return "\n".join(csv_row(row) for row in rows).strip()
