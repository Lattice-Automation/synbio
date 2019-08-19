"""Test Design classes."""

from typing import Iterable, List
import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.designs import Plasmid, Combinatorial, CombinatorialBins, PlasmidLibrary


class TestDesign(unittest.TestCase):
    """Design specification."""

    def setUp(self):
        """Instantiate test records."""

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
        """Create a combinatorial design."""

        combinations = Combinatorial([self.r1, self.r2, self.r3])
        expected = [[self.r1, self.r2, self.r3]]

        self.assertEqual(self.seq_list(expected), self.seq_list(combinations))
        self.assertEqual("combinatorial", combinations.__name__)

    def test_combinatorial_bins(self):
        """Test instantiation and iteration over CombinatorialBins design."""

        bins = CombinatorialBins([[self.r1, self.r2], [self.r3, self.r4]])
        expected = [
            [self.r1, self.r3],
            [self.r1, self.r4],
            [self.r2, self.r3],
            [self.r2, self.r4],
        ]

        self.assertEqual(self.seq_list(expected), self.seq_list(list(bins)))
        self.assertEqual("combinatorial_bins", bins.__name__)

    def test_get_all_records(self):
        """Get all unique records out of a design."""

        def seq_only(records: List[SeqRecord]) -> List[str]:
            return [str(r.seq) for r in records]

        plasmid = Plasmid([self.r1, self.r2, self.r3, self.r4])
        comb = CombinatorialBins([[self.r1, self.r2], [self.r3, self.r4]])

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

        comb = CombinatorialBins()
        self.assertEqual(0, len(comb))
        comb.append([self.r1, self.r3])
        self.assertEqual(1, len(comb))
        comb.append(self.r2)
        self.assertEqual(2, len(comb))
        self.assertIsInstance(comb.bins[1], list)

    def test_library(self):
        """Create and traverse a library of SeqRecords/SeqRecord-lists."""

        lib = PlasmidLibrary([[self.r1, self.r2]])
        lib.append([self.r3, self.r4])
        lib.append(self.r5)

        self.assertEqual("[3 x [SeqRecord]]", str(lib))

        lib_list = list(lib)
        self.assertEqual(3, len(lib_list))
        self.assertTrue(all(isinstance(l, list) for l in lib_list))
        self.assertEqual(
            self.seq_list([[self.r1, self.r2]]), self.seq_list([lib_list[0]])
        )

    def seq_list(self, records: Iterable[List[SeqRecord]]):
        """Map a list of list of SeqRecords to a list of Seq strings (comparable)."""

        return [r.seq for l in records for r in l]
