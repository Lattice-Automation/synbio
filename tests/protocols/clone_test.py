"""Test Clone protocol."""

import os
import unittest

from Bio.Restriction import BamHI, NotI
from Bio.SeqIO import parse

from synbio.designs import Plasmid
from synbio.protocols import Clone

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "data", "cloning")
OUT_DIR = os.path.join(DIR_NAME, "..", "output")


class TestCloneProtocol(unittest.TestCase):
    """Test Clone protocol."""

    def test_clone_protocol(self):
        """Run Clone protocol.
        
        Using sequences from https://benchling.com/tutorials/31/digest-and-ligate"""

        backbone = next(
            parse(os.path.join(TEST_DIR, "pdusk.gb"), "genbank")
        )  # has KanR in backbone
        insert = next(parse(os.path.join(TEST_DIR, "pdsred2.gb"), "genbank"))
        expected = next(parse(os.path.join(TEST_DIR, "pdusk-dsred2.gb"), "genbank"))

        protocol = Clone(
            "Cloning Protocol",
            enzymes=[BamHI, NotI],
            design=Plasmid([insert, backbone], linear=False),  # circular plasmids
            include=["KanR"],
        )
        protocol.run()

        # test output sequence
        self.assertEqual(3, len(protocol.output))
        self.assertTrue(
            any(
                a.seq in (expected + expected).seq
                or a.seq in (expected + expected).seq.reverse_complement()
                for a in protocol.output
            )
        )

        # export human and robotic instructions
        protocol.to_csv(os.path.join(OUT_DIR, "clone.csv"))
        protocol.to_txt(os.path.join(OUT_DIR, "clone_protocol.txt"))
        protocol.to_fasta(os.path.join(OUT_DIR, "clone.fasta"))
        protocol.to_genbank(os.path.join(OUT_DIR, "clone.gb"))
        protocol.to_picklists(
            os.path.join(OUT_DIR, "clone.tecan.gwl"), platform="tecan"
        )
        protocol.to_picklists(
            os.path.join(OUT_DIR, "clone.labcyte.gwl"), platform="labcyte"
        )
