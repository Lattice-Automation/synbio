"""Test Gibson Assembly design steps."""

import os
import unittest

from Bio import SeqIO
from Bio.SeqIO import parse

from synbio import Protocol, Combinatorial
from synbio.composite import Gibson

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.join(DIR_NAME, "..", "..", "goldengate")
OUT_DIR = os.path.join(DIR_NAME, "..", "..", "output")


class TestGibson(unittest.TestCase):
    """Test Gibson class methods."""

    pass
