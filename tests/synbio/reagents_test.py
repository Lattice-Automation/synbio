"""Test Reagents"""

import unittest

from synbio.containers import Reagent


class TestReagents(unittest.TestCase):
    """Test Reagent class."""

    def test_hash(self):
        """Hash reagent."""

        self.assertEqual(Reagent("water"), Reagent("water"))
