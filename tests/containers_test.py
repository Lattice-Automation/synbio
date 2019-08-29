"""Test Conatiner classes and functions
"""

import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Restriction import EcoRI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.containers import content_id, Layout, Well, Species, Reagent


class TestContainers(unittest.TestCase):
    """Container tests"""

    def test_content_id(self):
        """Create/find the unique IDs of container contents"""

        enzyme = EcoRI
        record = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna))

        self.assertEqual(str(EcoRI), content_id(enzyme))
        self.assertEqual(str(record.seq), content_id(record))
        self.assertEqual("e_coli", content_id(Species("e_coli")))
        self.assertEqual("water", content_id(Reagent("water")))

        with self.assertRaises(TypeError):
            content_id("")

    def test_layout(self):
        """Layout creation and serialization to CSV format."""

        plate_contents = [Well(SeqRecord(Seq("ATGATAGAT"))) for i in range(140)]
        plates = Layout(plate_contents)

        plate_csv = plates.to_csv()

        self.assertIn("Plate:1", plate_csv)
        self.assertIn(",,Plate:2", plate_csv)
        self.assertIn("B,", plate_csv)  # row label
        self.assertIn(",".join([str(n) for n in range(1, 13)]), plate_csv)  # col labels
        self.assertGreaterEqual(len(plate_csv.split("\n")), 9)

        # test plates name, well index (within plates) and well name (within plates)
        self.assertEqual("Plate:1", plates.container_to_plate_name[plate_contents[1]])
        self.assertEqual("Plate:1", plates.container_to_plate_name[plate_contents[95]])
        self.assertEqual("Plate:2", plates.container_to_plate_name[plate_contents[96]])

        self.assertEqual("Plate:2", plates.container_to_plate_name[plate_contents[105]])
        self.assertEqual(2, plates.container_to_well_index[plate_contents[1]])
        self.assertEqual("B1", plates.container_to_well_name[plate_contents[1]])
        self.assertEqual(2, len(plates))  # two plates

        # STILL just two plates, no reagents in contents
        plates = Layout(plate_contents, separate_reagents=True)
        self.assertEqual(2, len(plates))

        # Three plates, reagents are in a separate plate
        plate_reagent_contents = plate_contents + [
            Well(Reagent("water")),
            Well(Reagent("mix")),
        ]
        plates = Layout(plate_reagent_contents, separate_reagents=True)
        self.assertEqual(3, len(plates))

    def test_hash(self):
        """Hash containers."""

        well_set = {Well(), Well()}
        self.assertEqual(2, len(well_set))

    def test_withdraw(self):
        """Withdraw from a well, know when 'empty'"""

        c1 = Well([Reagent("e coli"), Reagent("water")], [50.0, 30.0])
        c1.volume_max = 100.0

        self.assertFalse(c1.empty())
        self.assertEqual(80.0, c1.volume())

        c1.withdraw(20.0)

        self.assertFalse(c1.empty())

        self.assertEqual(80.0, c1.volume())  # didn't change
        self.assertEqual(20.0, c1.withdrawn)
        self.assertTrue(c1.empty(70.0))  # there isn't 70 uL left

        with self.assertRaises(RuntimeError):
            c1.withdraw(100.0)

    def test_init_err(self):
        """Throw an error if init arguments are invalid."""

        self.assertEqual(
            Well(
                contents=[Reagent("water"), Reagent("water"), Reagent("water")],
                volumes=[2.0, 2.0, 0.0],
            ).volumes,
            [2.0, 2.0, 0.0],
        )

    def test_lt(self):
        """Compare two containers to see which should come first."""

        c1 = Well([Reagent("water")])
        c2 = Well([SeqRecord(Seq("ATGATAGAT"))])
        c3 = Well([Reagent("assembly-mix")])

        self.assertLess(c2, c1)
        self.assertLess(c3, c1)
