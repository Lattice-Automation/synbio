"""Generate Hamilton picklists."""

from typing import Dict, List, Union

from ..containers import Layout
from ..instructions import Instruction


def to_hamilton(
    instruction: Instruction, existing_plates: int, separate_reagents: bool = False
) -> str:
    """Create a picklist for a Hamilton robot

    Args:
        instruction: A single step's instruction
        existing_plates: The number of plates before these in protocol

    Keyword Args:
        separate_reagents: Whether to separate reagent plate from other wells

    Returns:
        The picklist in CSV string format
    """

    if not instruction.transfers:
        raise ValueError(f"no transfers in Instruction: {instruction}")

    # set each wells location in setup and destination plates
    plates = Layout.from_instruction(
        instruction,
        src_containers=True,
        existing_plates=existing_plates,
        separate_reagents=separate_reagents,
    )

    columns = ["LabID", "SourceID", "TargetID"]

    rows: List[Dict[str, Union[str, float]]] = []
    for transfer in instruction.transfers:
        src = transfer.src
        dest = transfer.dest

        rows.append(
            {
                "LabID": plates.container_to_plate_name[src],
                "SourceID": plates.container_to_well_name[src],
                "TargetID": plates.container_to_well_name[dest],
            }
        )

    def csv_row(row: Dict[str, Union[str, float]]) -> str:

        return ",".join(str(row.get(c, "")) for c in columns)

    csv = ",".join(columns) + "\n"
    csv += "\n".join(csv_row(row) for row in rows).strip()
    return csv
