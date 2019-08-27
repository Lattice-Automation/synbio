"""Test Protocol methods and inspection
"""

import unittest

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.containers import Reagent, Well, Fridge
from synbio.designs import Plasmid
from synbio.instructions import Instruction, Transfer
from synbio.protocol import Protocol
from synbio.steps import Pipette


class TestProtocol(unittest.TestCase):
    """Test Protocol class methods."""

    def setUp(self):
        """Create a mock Protocol."""

        self.r1 = SeqRecord(Seq("AGATAGACCAGAAGATAGA", unambiguous_dna))
        self.r2 = SeqRecord(Seq("GGGGGGAGAGAAACCCACAATAT", unambiguous_dna))

        self.protocol = Protocol(
            name="Mock Plasmid Design", design=Plasmid(records=[self.r1, self.r2])
        )

    def test_init(self):
        """Create a new Protocol."""

        protocol = Protocol("test plasmid", Plasmid([self.r1, self.r2]))

        self.assertEqual("test plasmid", protocol.name)
        self.assertEqual(0, len(protocol))  # no steps yet

    def test_add(self):
        """Add a new step to a Protocol."""

        curr_steps = len(self.protocol)

        self.protocol.add(Pipette([]))
        self.assertEqual(curr_steps + 1, len(self.protocol))

        with self.assertRaises(TypeError):
            self.protocol.add([])

    def test_str(self):
        """Create a string representation of a protocol."""

        p_name = self.protocol.name
        p_type = "plasmid"

        self.assertIn(p_name, str(self.protocol))
        self.assertIn(p_type, str(self.protocol))
        self.assertIn("design:", str(self.protocol))

    def test_inputs(self):
        """Get a name,value map for protocol inputs."""

        protocol = Protocol("input test", Plasmid([]))
        protocol.instructions.append(
            Instruction(
                transfers=[
                    Transfer(src=Fridge(Reagent("e coli")), dest=Well(), volume=30.0)
                ]
            )
        )
        protocol.instructions.append(
            Instruction(
                transfers=[
                    Transfer(
                        src=Fridge(
                            SeqRecord(
                                Seq("GGGGGGAGAGAAACCCACAATAT", unambiguous_dna),
                                id="mock_part",
                            )
                        ),
                        dest=Well(),
                        volume=50.0,
                    )
                ]
            )
        )

        inputs = protocol.input
        self.assertIn("mock_part", inputs.keys())
        self.assertIn("e coli", inputs.keys())
