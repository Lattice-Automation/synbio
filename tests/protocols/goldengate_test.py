"""Test GoldenGate design method."""

import os
import unittest

from Bio import SeqIO
from Bio.Restriction import BsaI
from Bio.SeqIO import parse

from synbio.designs import Combinatorial
from synbio.protocols import GoldenGate

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "data", "goldengate")
OUT_DIR = os.path.join(DIR_NAME, "..", "output")


class TestGoldenGate(unittest.TestCase):
    """Test GoldenGate class methods."""

    def test_goldengate(self):
        """End to end GoldenGate test."""

        # read in all the records
        records = []
        for (_, _, filenames) in os.walk(TEST_DIR):
            for file in filenames:
                gb = os.path.join(TEST_DIR, file)
                if not gb.endswith(".gb"):
                    continue
                for record in parse(gb, "genbank"):
                    records.append(record)

        record_sets = []
        for f_type in ["promoter", "RBS", "CDS", "terminator"]:

            def test(r):
                return any(
                    f.type == f_type and f.location.start < 50 for f in r.features
                )

            new_bin = [r for r in records if test(r)][:5]
            record_sets.append(new_bin)  # add a new bin

        records = [r for record_set in record_sets for r in record_set] + [
            self.read("DVK_AE.gb")
        ]

        # create a protocol, add GoldenGate as the sole protocol step, and run
        protocol = GoldenGate(
            name="Combinatorial GoldenGate",
            design=Combinatorial(records, linear=False),  # circular plasmid records
            include=["KanR"],
        )
        protocol.run()

        # export human and robotic instructions
        protocol.to_csv(os.path.join(OUT_DIR, "gg.csv"))
        protocol.to_txt(os.path.join(OUT_DIR, "gg_protocol.txt"))
        protocol.to_fasta(os.path.join(OUT_DIR, "gg.fasta"))
        protocol.to_genbank(os.path.join(OUT_DIR, "gg.gb"))
        protocol.to_picklists(os.path.join(OUT_DIR, "gg.tecan.gwl"), platform="tecan")
        protocol.to_picklists(
            os.path.join(OUT_DIR, "gg.labcyte.gwl"), platform="labcyte"
        )

    def test_lincoln(self):
        """Create a set of plasmids from a Combinatorial assembly."""

        records = []
        test_dir = os.path.join(TEST_DIR, "lincoln")
        for (_, _, filenames) in os.walk(test_dir):
            for file in filenames:
                if not file.endswith(".fa"):
                    continue
                test_file = os.path.join(test_dir, file)
                for record in parse(test_file, "fasta"):
                    records.append(record)

        design = Combinatorial(records)
        protocol = GoldenGate(design, enzymes=[BsaI], separate_reagents=True)
        protocol.run()

        csv_output = protocol.to_csv(os.path.join(OUT_DIR, "lincoln.layout.csv"))
        protocol.to_picklists(
            os.path.join(OUT_DIR, "lincoln.hamilton.csv"), platform="hamilton"
        )

        self.assertTrue(protocol.output)
        self.assertEqual(1, csv_output.count("Plate:1"))
        self.assertIn(
            ",,,,,,,,,,,,,A", csv_output
        )  # empty wells in plate 1 before plate2

    def read(self, filename):
        """Read in a single Genbank file from the test directory."""

        return next(SeqIO.parse(os.path.join(TEST_DIR, filename), "genbank"))
