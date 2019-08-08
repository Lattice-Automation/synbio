"""Test GoldenGate design method."""

import os
import unittest

from Bio import SeqIO
from Bio.SeqIO import parse

from synbio import Protocol, Combinatorial
from synbio.composite import GoldenGate

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "goldengate")
OUT_DIR = os.path.join(DIR_NAME, "..", "..", "output")


class TestGoldenGate(unittest.TestCase):
    """Test GoldenGate class methods."""

    def test_goldengate(self):
        """End to end GoldenGate test."""

        # create the test directory for file output
        if not os.path.exists(OUT_DIR):
            os.mkdir(OUT_DIR)

        # create a library design with multiple "bins"
        design = Combinatorial()

        # read in all the records
        records = []
        for (_, _, filenames) in os.walk(TEST_DIR):
            for file in filenames:
                gb_filename = os.path.join(TEST_DIR, file)
                for record in parse(gb_filename, "genbank"):
                    records.append(record)

        for f_type in ["promoter", "RBS", "CDS", "terminator"]:

            def test(r):
                return any(
                    f.type == f_type and f.location.start < 50 for f in r.features
                )

            new_bin = [r for r in records if test(r)][:10]

            design.append(new_bin)  # add a new bin

        # add a backbone
        design.append(self.read("DVK_AE.gb"))

        # create a protocol, add GoldenGate as the sole composite step, and run
        protocol = Protocol(name="Combinatorial GoldenGate", design=design)
        protocol.add(GoldenGate(resistance="KanR"))
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

    def read(self, filename):
        """Read in a single Genbank file from the test directory."""

        return next(SeqIO.parse(os.path.join(TEST_DIR, filename), "genbank"))
