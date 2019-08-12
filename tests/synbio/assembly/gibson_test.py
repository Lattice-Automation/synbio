"""Gibson testing."""

import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from synbio.assembly import gibson
from synbio.assembly.gibson import _hamming_set


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
        """Create primers extending in 3' direction when there's an offtarget."""

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

    def test_hamming_set(self):
        """Create a set of off-by-one DNA sequences."""

        hset = _hamming_set("ATG")

        self.assertEqual(10, len(hset))
        self.assertTrue(all(len(s) == 3 for s in hset))
        self.assertIn("TTG", hset)
        self.assertNotIn("TAG", hset)

    def _run_and_verify(self, records):
        """Verify the primers and plasmid sequence are correct."""

        plasmid, primers = gibson(records)

        doubled_seq = "".join(str(r.seq) for r in records + records)
        self.assertIsInstance(plasmid, SeqRecord)
        self.assertIn(str(plasmid.seq), doubled_seq)

        for i, primer in enumerate(primers):
            fwd = primer.fwd
            rev = primer.rev

            record = records[i]

            # primers' sequences are in the final plasmid
            self.assertIn(fwd, plasmid.seq)
            self.assertIn(rev.reverse_complement(), plasmid.seq + plasmid.seq)

            self.assertIn(fwd[-10:], record.seq)
            self.assertIn(rev[-10:].reverse_complement(), record.seq)

        return plasmid, primers
