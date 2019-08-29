"""Test Protocol Steps."""

import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio import Protocol
from synbio.designs import Combinatorial
from synbio.containers import Well, Reagent
from synbio.designs import Plasmid
from synbio.steps import Pipette, Setup


class TestSteps(unittest.TestCase):
    """Test individual Protocol steps."""

    def setUp(self):
        """Instantiate a test Protocol with SeqRecords in containers."""

        self.r1 = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna), id="r1")
        self.r2 = SeqRecord(Seq("GGGGGGAGAGAAACCCACAATAT", unambiguous_dna), id="r2")
        self.r3 = SeqRecord(Seq("TTTTGAGAGAGATTTAGAGATA", unambiguous_dna), id="r3")
        self.r4 = SeqRecord(Seq("GAGATAGATAGACCAGATAGA", unambiguous_dna), id="r4")

        self.protocol = Protocol("", Combinatorial())

        self.protocol.containers = [
            Well([self.r1, self.r2], volumes=[50.0, 50.0]),
            Well([self.r3, self.r4], volumes=[50.0, 50.0]),
        ]

    def test_setup(self):
        """Create a setup plate with multiple wells."""

        # the water here requires 300.0 uL
        # this exceeds the max for a well
        d1 = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna), id="r1")
        r1 = Reagent("water")

        w1 = Well([d1, r1], [50.0, 100.0])
        w2 = Well([d1, r1], [50.0, 100.0])

        target = [w1, w2]

        protocol = Protocol("", Plasmid())

        protocol = Setup(target)(protocol)

        setup_wells = protocol.containers

        self.assertEqual(3, len(setup_wells))

        water_wells = [w for w in setup_wells if r1 in w]
        self.assertEqual(2, len(water_wells))

        dna_wells = [w for w in setup_wells if d1 in w]
        self.assertEqual(1, len(dna_wells))

        self.assertTrue(all(w.volumes and w.volumes[0] for w in setup_wells))

    def test_pipette(self):
        """Test Pipette step to move SeqRecords into new wells."""

        def contents(containers):  # get lists of ids (SeqRecords not comparable)
            return [c.id for container in containers for c in container]

        target = [
            Well([self.r1, self.r3], volumes=[20.0, 20.0]),
            Well([self.r2, self.r4], volumes=[20.0, 20.0]),
        ]

        pipette_step = Pipette(target=target)

        self.assertNotEqual(contents(self.protocol.containers), contents(target))

        pipette_step(self.protocol)

        self.assertEqual(contents(self.protocol.containers), contents(target))
