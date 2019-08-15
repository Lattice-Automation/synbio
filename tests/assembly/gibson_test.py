"""Gibson testing."""

import logging
import os
import unittest

from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from synbio.assembly import gibson
from synbio.assembly.gibson import _hamming_set
from synbio.primers import MIN_PRIMER_LEN, MAX_PRIMER_LEN

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "data", "gibson")


class TestGibson(unittest.TestCase):
    """Test Gibson class"""

    def test_gibson(self):
        """Create primers for a Gibson assembly."""

        s1 = Seq(
            "GCGTTTTATAGAAGCCTAGGGGAACAGATTGGTCTAATTAGCTTAAGAGAGTAAATTCTGGGATCATTCAGTAGTAATCACAAATTTACGGTGGGGCTTTTTTGGCGGATCTTTACAGATACTAACCAGGTGATTTCAACTAATTTAGTTGACGATTTAGGCGCGCTATCCCGTAATCTTCAAATTAAAACATAGCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGTTATTCGCGCCCACTCTCCCATTTATCCGCGCAAGCGGATGCGATGCGATTGCCCGCT"
        )
        s2 = Seq(
            "AAGATATTCTTACGTGTAACGTAGCTAAGTATTCTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTCATTGGGGGCTTATACAGGCGTAGACTACAATGGGCCCAACTCAATCACAGCTCGAGCGCCTTGAATAACATACTCATCTCTATACATTCTCGACAATCTATCGAGCGAGTCGATTATCAACGGGTGTGTTGCAGTTTTAATCTCTTGCCAGCATTGTAATAGCCACCAAGAGATTGATGATAGTCATGGGTGCTGAGCTGAGA"
        )
        s3 = Seq(
            "CGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTGTGCGCATCACTCTCAATGTACAAGCAACCCAAGAAGGCTGAGCCTGGACTCAACCGGTTGCTGGGTGAACTCCAGACTCGGGGCGACAACTCTTCATACATAGAGCAAGGGCGTCGAACGGTCGTGAAAGTCTTAGTACCGCACGTACCAACTTACTGAGGATATTGCCTGAAGCTGTACCGTTTTAGGGGGGGAAGGTTGAAGATCTCCTCTTCTCATGACTGAACTCGCGAGGGCCG"
        )
        r1 = SeqRecord(s1)
        r2 = SeqRecord(s2)
        r3 = SeqRecord(s3)

        self._run_and_verify([r1, r2, r3])

    def test_gibson_offtarget_primer(self):
        """Create Gibson primers when there's offtarget in one's end (1)."""

        s1 = Seq(  # GCGTTTTATAGAAGCCTAGGGGAAC shows up twice
            "GCGTTTTATAGAAGCCTAGGGGAACAGATTGGTCTAATTAGCTTAAGAGAGTAAATGGCGTTTTATAGAAGCCTAGGGGAACCTGGGATCATTCAGTAGTAATCACAAATTTACGGTGGGGCTTTTTTGGCGGATCTTTACAGATACTAACCAGGTGATTTCAACTAATTTAGTTGACGATTTAGGCGCGCTATCCCGTAATCTTCAAATTAAAACATAGCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGTTATTCGCGCCCACTCTCCCATTTATCCGCGCAAGCGGATGCGATGCGATTGCCCGCT"
        )
        s2 = Seq(
            "AAGATATTCTTACGTGTAACGTAGCTAAGTATTCTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTCATTGGGGGCTTATACAGGCGTAGACTACAATGGGCCCAACTCAATCACAGCTCGAGCGCCTTGAATAACATACTCATCTCTATACATTCTCGACAATCTATCGAGCGAGTCGATTATCAACGGGTGTGTTGCAGTTTTAATCTCTTGCCAGCATTGTAATAGCCACCAAGAGATTGATGATAGTCATGGGTGCTGAGCTGAGA"
        )
        s3 = Seq(
            "CGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTGTGCGCATCACTCTCAATGTACAAGCAACCCAAGAAGGCTGAGCCTGGACTCAACCGGTTGCTGGGTGAACTCCAGACTCGGGGCGACAACTCTTCATACATAGAGCAAGGGCGTCGAACGGTCGTGAAAGTCTTAGTACCGCACGTACCAACTTACTGAGGATATTGCCTGAAGCTGTACCGTTTTAGGGGGGGAAGGTTGAAGATCTCCTCTTCTCATGACTGAACTCGCGAGGGCCG"
        )
        r1 = SeqRecord(s1)
        r2 = SeqRecord(s2)
        r3 = SeqRecord(s3)

        self._run_and_verify([r1, r2, r3])

    def test_gibson_offtarget_primer2(self):
        """Create Gibson primers when there's offtarget in one's end (2)."""

        insert = next(parse(os.path.join(TEST_DIR, "BBa_K1649003.fa"), "fasta"))
        backbone = next(parse(os.path.join(TEST_DIR, "pDusk.fa"), "fasta"))

        plasmid, primer_pairs = gibson([insert, backbone])

        self.assertTrue(plasmid and primer_pairs)

    def test_hamming_set(self):
        """Create a set of off-by-one DNA sequences."""

        hset = _hamming_set("ATG")

        self.assertEqual(10, len(hset))
        self.assertTrue(all(len(s) == 3 for s in hset))
        self.assertIn("TTG", hset)
        self.assertNotIn("TAG", hset)

    def _run_and_verify(self, records):
        """Verify the primers and plasmid sequence are correct."""

        plasmid, primer_pairs = gibson(records)

        doubled_seq = "".join(str(r.seq) for r in records + records)
        self.assertIsInstance(plasmid, SeqRecord)
        self.assertIn(str(plasmid.seq), doubled_seq)

        for i, primers in enumerate(primer_pairs):
            self.assertGreater(len(primers.fwd), MIN_PRIMER_LEN)
            self.assertGreater(len(primers.rev), MIN_PRIMER_LEN)
            self.assertLess(len(primers.fwd), MAX_PRIMER_LEN)
            self.assertLess(len(primers.rev), MAX_PRIMER_LEN)

            # primers' sequences are in the final plasmid
            self.assertIn(primers.fwd, plasmid.seq)
            self.assertIn(primers.rev.reverse_complement(), plasmid.seq + plasmid.seq)

            record = records[i]
            self.assertIn(primers.fwd[-MIN_PRIMER_LEN:], record.seq)
            self.assertIn(
                primers.rev[-MIN_PRIMER_LEN:].reverse_complement(), record.seq
            )

        return plasmid, primer_pairs
