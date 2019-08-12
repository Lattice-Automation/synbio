"""Test Labcyte automation script generation."""

import unittest

from synbio.containers import Well
from synbio.instructions import Instruction, Transfer
from synbio.picklists import to_labcyte


class TestLabcyte(unittest.TestCase):
    """Labcyte picklist generation."""

    def test_to_labcyte(self):
        """Create a labcyte picklist."""

        src = Well()
        dest = Well()
        instruction = Instruction(transfers=[Transfer(src, dest, volume=5.0)])
        picklist = to_labcyte(instruction, 0)
        lines = picklist.split("\n")

        self.assertIn("Source Plate Barcode", lines[0])
        self.assertIn("Source Well", lines[0])
        self.assertIn("Destination Plate Barcode", lines[0])
        self.assertIn("Destination Well", lines[0])
        self.assertEqual("Plate:1,A1,Plate:2,A1,5000.0", lines[1])

        # max is 10 uL, so a 20 uL needs to be split in half
        instruction = Instruction(transfers=[Transfer(src, dest, volume=20.0)])
        picklist = to_labcyte(instruction, 0)
        lines = picklist.split("\n")

        self.assertEqual(3, len(lines))  # one transfer split into two
        self.assertEqual("Plate:1,A1,Plate:2,A1,10000.0", lines[1])
        self.assertEqual("Plate:1,A1,Plate:2,A1,10000.0", lines[2])

        with self.assertRaises(ValueError):
            to_labcyte(Instruction(), 0)  # needs to have transfers to use in picklist
