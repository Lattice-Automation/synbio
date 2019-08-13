"""Test Tecan automation script generation."""

import unittest

from synbio.containers import Well
from synbio.instructions import Instruction, Transfer
from synbio.picklists import to_tecan


class TestTecan(unittest.TestCase):
    """Tecan picklist generation."""

    def test_to_tecan(self):
        """Create a Tecan picklist."""

        src = Well()
        dest = Well()
        instruction = Instruction(transfers=[Transfer(src, dest, volume=5.0)])
        picklist = to_tecan(instruction, 0)
        lines = picklist.split("\n")

        self.assertEqual(3, len(lines))
        self.assertEqual("A;Plate:1;;;1;;5.0;;;", lines[0])
        self.assertEqual("D;Plate:2;;;1;;5.0;;;", lines[1])
        self.assertEqual("W;;;;;;;;;", lines[2])

        with self.assertRaises(ValueError):
            to_tecan(Instruction(), 0)  # needs to have transfers to use in picklist
