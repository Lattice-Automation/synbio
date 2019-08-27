"""Test Gibson Assembly design steps."""

import os
import unittest

from Bio.SeqIO import parse

from synbio.designs import Plasmid
from synbio.protocols import Gibson

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "data", "gibson")
OUT_DIR = os.path.join(DIR_NAME, "..", "output")


class TestGibsonProtocol(unittest.TestCase):
    """Test Gibson protocol."""

    def test_gibson_protocol(self):
        """Run a Gibson Assembly protocol."""

        insert = next(parse(os.path.join(TEST_DIR, "pdsred2.fa"), "fasta"))
        backbone = next(parse(os.path.join(TEST_DIR, "pDusk.fa"), "fasta"))

        protocol = Gibson(
            design=Plasmid([insert, backbone]), name="Gibson plasmid assembly"
        )
        protocol.run()

        # export human and robotic instructions
        protocol.to_csv(os.path.join(OUT_DIR, "gibson.csv"))
        protocol.to_txt(os.path.join(OUT_DIR, "gibson_protocol.txt"))
        protocol.to_fasta(os.path.join(OUT_DIR, "gibson.fasta"))
        protocol.to_genbank(os.path.join(OUT_DIR, "gibson.gb"), split=True)
        protocol.to_picklists(
            os.path.join(OUT_DIR, "gibson.tecan.gwl"), platform="tecan"
        )
        protocol.to_picklists(
            os.path.join(OUT_DIR, "gibson.labcyte.gwl"), platform="labcyte"
        )
