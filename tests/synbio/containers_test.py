"""Test Conatiner classes and functions
"""

import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Emboss.Primer3 import Primers
from Bio.Restriction import EcoRI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.containers import content_id, Plates, Well, Species, Reagent, Fridge


class TestContainers(unittest.TestCase):
    """Container tests"""

    def test_content_id(self):
        """Create/find the unique IDs of container contents"""

        enzyme = EcoRI
        record = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna))
        primers = Primers()
        primers.forward_seq = "AGATGATAGAT"
        primers.reverse_seq = "AGATGAGCCAG"

        expected_primer_id = f"primers:{primers.forward_seq};{primers.reverse_seq}"

        self.assertEqual(str(EcoRI), content_id(enzyme))
        self.assertEqual(str(record.seq), content_id(record))
        self.assertEqual(expected_primer_id, content_id(primers))
        self.assertEqual("e_coli", content_id(Species("e_coli")))
        self.assertEqual("water", content_id(Reagent("water")))

        with self.assertRaises(TypeError):
            content_id("")

    def test_plates(self):
        """Plate creation and serialization to CSV format."""

        plate_contents = [
            Well(SeqRecord(Seq("ATGATAGAT"), id=str(i))) for i in range(140)
        ]
        plates = Plates(plate_contents, rows=8, cols=12)

        plate_csv = plates.to_csv()

        self.assertIn("Plate:1", plate_csv)
        self.assertIn(",,Plate:2", plate_csv)
        self.assertIn("B,", plate_csv)  # row label
        self.assertIn(",".join([str(n) for n in range(1, 13)]), plate_csv)  # col labels
        self.assertGreaterEqual(len(plate_csv.split("\n")), 9)

        # test plates name, well index (within plates) and well name (within plates)
        content_2 = plate_contents[1]

        self.assertEqual("Plate:1", plates.container_to_plate_name[content_2])
        self.assertEqual("Plate:1", plates.container_to_plate_name[plate_contents[95]])
        self.assertEqual("Plate:2", plates.container_to_plate_name[plate_contents[96]])
        self.assertEqual("Plate:2", plates.container_to_plate_name[plate_contents[105]])
        self.assertEqual(2, plates.container_to_well_index[content_2])
        self.assertEqual("B1", plates.container_to_well_name[content_2])
        self.assertEqual(2, len(plates))  # two plates

    def test_hash(self):
        """Hash containers."""

        self.assertEqual(hash(Fridge()), hash(Fridge()))
        self.assertNotEqual(Well(), Well())

    def test_withdraw(self):
        """Withdraw from a well, know when 'empty'"""

        c1 = Well([Reagent("e coli"), Reagent("water")], [50.0, 30.0], volume_max=100.0)

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

        with self.assertRaises(ValueError):
            Well(volumes=[40, 80, 90], volume_max=50.0)

        with self.assertRaises(ValueError):  # only one reagent, two volumes
            Well(volumes=[2.0, 5.0], contents=Reagent("e coli"))

    def test_lt(self):
        """Compare two containers to see which should come first."""

        c1 = Well([Reagent("water")])
        c2 = Well([SeqRecord(Seq("ATGATAGAT"))])
        c3 = Well([Reagent("assembly-mix")])

        self.assertLess(c2, c1)
        self.assertLess(c3, c1)
