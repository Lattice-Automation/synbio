"""Mix testing."""

import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from synbio import Mix, Reagent, Species
from synbio.containers import content_id


class TestMix(unittest.TestCase):
    """Test Mix class"""

    def test_mix_init(self):
        """Fail to create a Mix with nonsensical params."""

        with self.assertRaises(ValueError):
            Mix({}, fill_with=Reagent("water"))  # no fill volume

        with self.assertRaises(ValueError):
            Mix({}, fill_to=50.0)  # no fill content

    def test_call(self):
        """Call Mix on a list of a contents."""

        water = Reagent("water")
        master_mix = Reagent("master mix")

        mix = Mix({master_mix: 4.0, SeqRecord: 2.0}, fill_with=water, fill_to=20.0)

        f1 = SeqRecord(Seq("GGAGttgac"))
        f2 = SeqRecord(Seq("acagtctca"))
        contents, volumes = mix([f1, f2])

        self.assertEqual(
            self._content_ids([f1, f2, master_mix, water]), self._content_ids(contents)
        )
        self.assertEqual([2.0, 2.0, 4.0, 12.0], volumes)

    def _content_ids(self, contents):
        return [content_id(c) for c in contents]
