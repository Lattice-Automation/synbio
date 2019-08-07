"""Test Design classes
"""

from typing import Iterable, List
import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.design import Plasmid, Combinatorial


class TestDesign(unittest.TestCase):
    def setUp(self):
        """Instantiate test records.
        """

        self.r1 = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna))
        self.r2 = SeqRecord(Seq("GGGGGGAGAGAAACCCACAATAT", unambiguous_dna))
        self.r3 = SeqRecord(Seq("TTTTGAGAGAGATTTAGAGATA", unambiguous_dna))
        self.r4 = SeqRecord(Seq("GAGATAGATAGACCAGATAGA", unambiguous_dna))

        # same as r3
        self.r5 = SeqRecord(Seq("GGGGGGAGAGAAACCCACAATAT", unambiguous_dna))

    def test_plasmid(self):
        """Instantiate and iterate over a Plasmid design."""

        plasmid = Plasmid([self.r1, self.r2, self.r3])
        expected = [[self.r1, self.r2, self.r3]]

        self.assertEqual(self.seq_list(expected), self.seq_list(plasmid))
        self.assertEqual("plasmid", plasmid.__name__)

    def test_combinatorial(self):
        """Test instantiation and iteration over Combinatorial design."""

        comb = Combinatorial([[self.r1, self.r2], [self.r3, self.r4]])
        expected = [
            [self.r1, self.r3],
            [self.r1, self.r4],
            [self.r2, self.r3],
            [self.r2, self.r4],
        ]

        self.assertEqual(self.seq_list(expected), self.seq_list(list(comb)))
        self.assertEqual("combinatorial", comb.__name__)

    def test_get_all_records(self):
        """Get all unique records out of a design."""

        def seq_only(records: List[SeqRecord]) -> List[str]:
            return [str(r.seq) for r in records]

        plasmid = Plasmid([self.r1, self.r2, self.r3, self.r4])
        comb = Combinatorial([[self.r1, self.r2], [self.r3, self.r4]])

        plasmid_records = set(seq_only(plasmid.get_all_records()))
        comb_records = set(seq_only(comb.get_all_records()))

        self.assertEqual(4, len(comb_records))
        self.assertEqual(plasmid_records, comb_records)

    def test_plasmid_append(self):
        """Add a new SeqRecord to a plasmid design."""

        plasmid = Plasmid()
        self.assertEqual(0, len(plasmid))
        plasmid.append(self.r1)
        plasmid.append(self.r2)
        self.assertEqual(2, len(plasmid))
        plasmid.append([self.r3, self.r4])
        self.assertEqual(4, len(plasmid))

    def test_combinatorial_append(self):
        """Add a new bin and record to a combinatorial design."""

        comb = Combinatorial()
        self.assertEqual(0, len(comb))
        comb.append([self.r1, self.r3])
        self.assertEqual(1, len(comb))
        comb.append(self.r2)
        self.assertEqual(2, len(comb))
        self.assertIsInstance(comb.bins[1], list)

    def seq_list(self, records: Iterable[List[SeqRecord]]):
        """Map a list of SeqRecords to a list of Seq strings (comparable)."""

        return [r.seq for l in records for r in l]

