"""Test oligo functions."""

import unittest

from synbio.oligos import calc_tm


class TestOligos(unittest.TestCase):
    """Test oligo functions"""

    def test_mix_init(self):
        """Test oligo tm calculation."""

        # values are from Table 1 of IDT's:
        # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
        # with a 1.5mM Mg concentration which looks typical according to NEB
        experimental_tms = {
            "GGGACCGCCT": 51.9,
            "CCATTGCTACC": 42.7,
            "GCAGTGGATGTGAGA": 55.1,
            "CTGGTCTGGATCTGAGAACTTCAGG": 67.7,
            "CTTAAGATATGAGAACTTCAACTAATGTGT": 59.7,
            "AGTCTGGTCTGGATCTGAGAACTTCAGGCT": 71.6,
        }

        for seq, actual in experimental_tms.items():
            calc = calc_tm(seq)
            actual = actual
            if abs(calc - actual) > 3:
                raise AssertionError(
                    f"large tm diff: got {calc} expected {actual} for {seq}"
                )
