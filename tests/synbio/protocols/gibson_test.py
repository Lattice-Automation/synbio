"""Test Gibson Assembly design steps."""

import os
import unittest

from Bio import SeqIO
from Bio.SeqIO import parse

from synbio import Plasmid
from synbio.protocols import Gibson

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "gibson")
OUT_DIR = os.path.join(DIR_NAME, "..", "..", "output")


class TestGibsonProtocol(unittest.TestCase):
    """Test Gibson protocol."""

    def test_gibson_protocol(self):
        """Run a Gibson Assembly protocol."""

        insert = next(parse(os.path.join(TEST_DIR, "BBa_K1085023.fa"), "fasta"))
        backbone = next(parse(os.path.join(TEST_DIR, "pSB1C3.fa"), "fasta"))

        protocol = Gibson("Gibson plasmid assembly", design=Plasmid([insert, backbone]))
        protocol.run()

        # export human and robotic instructions
        protocol.to_csv(os.path.join(OUT_DIR, "gibson.csv"))
        protocol.to_txt(os.path.join(OUT_DIR, "gibson_protocol.txt"))
        protocol.to_fasta(os.path.join(OUT_DIR, "gibson.fasta"))
        protocol.to_genbank(os.path.join(OUT_DIR, "gibson.gb"))
        protocol.to_picklists(
            os.path.join(OUT_DIR, "gibson.tecan.gwl"), platform="tecan"
        )
        protocol.to_picklists(
            os.path.join(OUT_DIR, "gibson.labcyte.gwl"), platform="labcyte"
        )
