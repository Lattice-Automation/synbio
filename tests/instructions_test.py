"""Test Instructions."""

import unittest

from synbio.containers import Well
from synbio.instructions import _display_time, Transfer


class InstructionsTest(unittest.TestCase):
    """Test instructions."""

    def test_display_time(self):
        """Create time messages from step durations in seconds."""

        self.assertEqual("for 1 hour, 1 second", _display_time(3601))
        self.assertEqual("hold", _display_time(0))
        self.assertEqual("hold", _display_time(-1))

    def test_transfer_split(self):
        """Splitup a transfer based on max-volume."""

        transfer = Transfer(src=Well(), dest=Well(), volume=1000)  # millileter

        max_volume = 50.0
        split = transfer.split(max_volume, 1.0)
        self.assertEqual(20, len(split))
        self.assertTrue(all(t.volume and t.volume <= max_volume for t in split))
