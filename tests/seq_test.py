"""Test sequence utility functions."""


import unittest

from synbio.seq import mutate


class TestSeq(unittest.TestCase):
    """Test sequence specific methods."""

    def test_mutate(self):
        """Create a set of mutated sequences off by 2 bp."""

        seq = "ATGC"

        mutants_1_edit = mutate(seq, edit_distance=1)
        mutants_2_edit = mutate(seq, edit_distance=2)

        self.assertIn(seq, mutants_1_edit)
        self.assertIn("TTGC", mutants_1_edit)
        self.assertNotIn("TAGC", mutants_1_edit)
        self.assertIn(seq, mutants_2_edit)
        self.assertIn("TTGC", mutants_2_edit)
        self.assertIn("TAGC", mutants_2_edit)  # 2 edit distance
        self.assertNotIn("TACC", mutants_2_edit)  # 3 edit distance
