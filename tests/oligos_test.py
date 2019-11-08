"""Test oligo functions."""

import unittest

from synbio.oligos import calc_tm, fold, _bulge, _pair, _hairpin


class TestOligos(unittest.TestCase):
    """Test oligo functions"""

    def test_calc_tm(self):
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
            self.assertAlmostEqual(calc, actual, delta=3)  # within 3 deg tm difference

    def test_fold(self):
        """Test DNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = {
            # "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,
            "AAGGGGTTGGTCGCCTCGACTAAGCGGCTGGATTCC": -2.5,
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.1,
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
            # the below is a three branched structure
            # "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,
            # "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC": -9.35
        }

        for seq, unafold_est in unafold_dgs.items():
            calc_dg = fold(seq)

            # accepting a 25% difference
            delta = abs(0.25 * unafold_est)
            self.assertAlmostEqual(unafold_est, calc_dg, delta=delta)

    def test_bulge(self):
        """Test delta G calc of a bulge."""

        # mock bulge of CAT on one side and AG on other
        # from pg 429 of SantaLucia, 2004
        pair = "CT/GA"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"  # nonsense sequence

        pair_dg = _bulge(pair, seq, 5, 7, 310.15)
        self.assertAlmostEqual(3.22, pair_dg, delta=0.4)

    def test_pair(self):
        """Test delta G of pairs with and without mismatches."""

        # from pg 429 of SantaLucia, 2004
        pair = "CT/GA"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"

        pair_dg = _pair(pair, seq, 5, 27, 310.15)

        self.assertAlmostEqual(-1.28, pair_dg, delta=0.1)

    def test_hairpin(self):
        """Test delta G of a hairpin structure."""

        # from page 428 of SantaLucia, 2004
        # hairpin = "CGCAAG"
        seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        k = 3
        j = 8
        hairpin_dg = _hairpin(seq, k, j, 310.15)
        self.assertAlmostEqual(0.67, hairpin_dg, delta=0.1)
