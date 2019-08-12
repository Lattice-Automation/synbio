"""Lab species and organisms."""

import unittest

from synbio.species import Species


class TestSpecies(unittest.TestCase):
    """Test Species"""

    def test_hash(self):
        """Hash species."""

        self.assertEqual(Species("e coli"), Species("e coli"))
