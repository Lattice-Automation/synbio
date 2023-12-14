from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from synbio.designs import Plasmid
from Bio.Restriction import BbsI

from synbio.protocols.clone import Clone
import unittest

class TestFlag(unittest.TestCase):

    def test_strict(self):
        with open('seq_records.fasta', 'r') as file:
            design = list(SeqIO.parse(file, 'fasta'))

        protocol = Clone(
                    name="test clone",
                    enzymes=[BbsI],
                    design=Plasmid(design, linear=False),
                    min_count=4
                )

        with self.assertRaises(Exception):
            protocol.run()

    def test_without_strict(self):
        with open('seq_records.fasta', 'r') as file:
            design = list(SeqIO.parse(file, 'fasta'))

        protocol = Clone(
                    name="test clone",
                    enzymes=[BbsI],
                    design=Plasmid(design, linear=False),
                    min_count=4,
                    strict=False
                )
                
        protocol.run()
        self.assertEqual(type(protocol.output[0]), SeqRecord)

if __name__ == '__main__':
    unittest.main()