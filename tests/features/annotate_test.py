"""Test SeqRecord annotation."""

import unittest

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from synbio.features.annotate import (
    _cull,
    _filter_hits,
    _get_features,
    _get_hits,
    _get_matches,
    _reduce_hits,
    annotate,
    Hit,
)


class TestAnnotate(unittest.TestCase):
    """Test feature annotate on a sequence."""

    def test_annotate(self):
        """Annotate a plasmid with known features: pUC"""

        seq = "gtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctgcaagcttaatgtagtcttatgcaatactcttgtagtcttgcaacatggtaacgatgagttagcaacatgccttacaaggagagaaaaagcaccgtgcatgccgattggtggaagtaaggtggtacgatcgtgccttattaggaaggcaacagacgggtctgacatggattggacgaaccactgaattgccgcattgcagagatattgtatttaagtgcctagctcgatacataaacgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagggacttgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcggcttgctgaagcgcgcacggcaagaggcgaggggcggcgactggtgagtacgccaaaaattttgactagcggaggctagaaggagagagatgggtgcgagagcgtcagtattaagcgggggagaattagatcgcgatgggaaaaaattcggttaaggccagggggaaagaaaaaatataaattaaaacatatagtatgggcaagcagggagctagaacgattcgcagttaatcctggcctgttagaaacatcagaaggctgtagacaaatactgggacagctacaaccatcccttcagacaggatcagaagaacttagatcattatataatacagtagcaaccctctattgtgtgcatcaaaggatagagataaaagacaccaaggaagctttagacaagatagaggaagagcaaaacaaaagtaagaccaccgcacagcaagcggccgctgatcttcagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagtggtgcagagagaaaaaagagcagtgggaataggagctttgttccttgggttcttgggagcagcaggaagcactatgggcgcagcgtcaatgacgctgacggtacaggccagacaattattgtctggtatagtgcagcagcagaacaatttgctgagggctattgaggcgcaacagcatctgttgcaactcacagtctggggcatcaagcagctccaggcaagaatcctggctgtggaaagatacctaaaggatcaacagctcctggggatttggggttgctctggaaaactcatttgcaccactgctgtgccttggaatgctagttggagtaataaatctctggaacagatttggaatcacacgacctggatggagtgggacagagaaattaacaattacacaagcttaatacactccttaattgaagaatcgcaaaaccagcaagaaaagaatgaacaagaattattggaattagataaatgggcaagtttgtggaattggtttaacataacaaattggctgtggtatataaaattattcataatgatagtaggaggcttggtaggtttaagaatagtttttgctgtactttctatagtgaatagagttaggcagggatattcaccattatcgtttcagacccacctcccaaccccgaggggacccagagagggcctatttcccatgattccttcatatttgcatatacgatacaagcctgttagagagataattagaattaatttgactgtaaacacaaagatattagtacaaaatacgtgacgtagaaagtaataatttcttgggtagtttgcagttttaaaattatgttttaaaatggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaaCACCGGCCCGCTTTGCATACGCCGTgtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgctttttttctcctatatattcattgtcctaattttaattcttgcctaatttcgtctactttaactttagcgttttgaacagattcaccaacacctataatccgtagcctaggttcagttccacttgggcgaacagcaaatcatgacttatcttctagataacggggagggcctatttcccatgattccttcatatttgcatatacgatacaaggctgttagagagataattagaattaatttgactgtaaacacaaagatattagtacaaaatacgtgacgtagaaagtaataatttcttgggtagtttgcagttttaaaattatgttttaaaatggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgcgcaactccatcgaagccgagtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttgaattctagatcttgagacaaatggcagtattcatccacaattttaaaagaaaaggggggattggggggtacagtgcaggggaaagaatagtagacataatagcaacagacatacaaactaaagaattacaaaaacaaattacaaaaattcaaaattttcgggtttattacagggacagcagagatccactttggcgccggctcgagtggctccggtgcccgtcagtgggcagagcgcacatcgcccacagtccccgagaagttggggggaggggtcggcaattgaaccggtgcctagagaaggtggcgcggggtaaactgggaaagtgatgtcgtgtactggctccgcctttttcccgagggtgggggagaaccgtatataagtgcagtagtcgccgtgaacgttctttttcgcaacgggtttgccgccagaacacaggtgtcgtgacgcgggatccgccaccatgaccgagtacaagcccacggtgcgcctcgccacccgcgacgacgtccccagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgacatcggcaaggtgtgggtcgcggacgacggcgccgcggtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggagtctcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaggcggccgagcgcgccggggtgcccgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgcctgaacgcgttaagtcgacaatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcgtcgactttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatctgctttttgcttgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtacgtatagtagttcatgtcatcttattattcagtatttataacttgcaaagaaatgaatatcagagagtgagaggaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggctctagctatcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctagggacgtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatgggacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaatttaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagacccc"
        record = SeqRecord(seq)

        annotated_record = annotate(record)

        self.assertEqual(0, len(record.features))
        self.assertTrue(len(annotated_record.features) > 1)

    def test_get_features_dna(self):
        """Get a list of DNA features for a SeqRecord."""

        mock_feature1 = SeqRecord(
            "TCCTCCCGGcagcaaaaaaGGGctcaagacccgttta".upper(),
            id="s1",
            name="mock1",
            annotations={"type": "promoter"},
        )
        mock_feature2 = SeqRecord(
            "taaacgggtcttgaggggttttttgc".upper(),
            id="s2",
            name="mock2",
            annotations={"type": "terminator"},
        )
        seq = "ggatatagtTCCTCCCGGcagcaaaaaacccctcaagacccgtttagaggccccaaggggttatgctagttattgctcagcggtggcagcagccaactcagcttcctttcgggctttgttagcagccggatctcagt".upper()
        kmer_map = {
            "TCCTCCCGG": [("s1", 0)],
            "CCTCCCGGC": [("s1", 1)],
            "CTCCCGGCA": [("s1", 2)],
            # taaacgggtcttgaggggttttttgc feature
            # rc: gcaaaaaacccctcaagacccgttta
            "TAAACGGGT": [("s2", 0)],
            "AAACGGGTC": [("s2", 1)],
            "AACGGGTCT": [("s2", 2)],
        }
        subject_map = {"s1": mock_feature1, "s2": mock_feature2}
        identity = 0.9

        features = _get_features(seq, kmer_map, subject_map, identity, True, False)
        feature1 = features[0]
        feature2 = features[1]

        self.assertEqual(2, len(features))
        self.assertEqual(SeqFeature, type(feature1))
        self.assertEqual(feature1.id, mock_feature1.name)
        self.assertEqual(
            "TCCTCCCGGcagcaaaaaacccctcaagacccgttta".upper(), feature1.extract(seq)
        )
        self.assertEqual(feature1.location.start, 9)
        self.assertEqual(feature1.location.end, 46)
        self.assertEqual(feature1.location.strand, 1)

        self.assertEqual(SeqFeature, type(feature2))
        self.assertEqual(feature2.id, mock_feature2.name)
        self.assertEqual("taaacgggtcttgaggggttttttgc".upper(), feature2.extract(seq))
        self.assertEqual(feature2.location.start, 20)
        self.assertEqual(feature2.location.end, 46)
        self.assertEqual(feature2.location.strand, -1)

    def test_get_features_protein(self):
        """Get a list of protein features for a SeqRecord."""

        # corresponds to agtTCCTCCCGGcagcaaaaaacc in the query seq
        mock_feature1 = SeqRecord(
            "SSSRQQKT", id="s1", name="s1", annotations={"type": "CDS"}
        )
        # corresponds to acgggtcttgaggggttttttgct in the query seq (reverse complement)
        # agcaaaaaacccctcaagacccgt in the template sequence
        mock_feature2 = SeqRecord(
            "TGLEGFFA", id="s2", name="s2", annotations={"type": "CDS"}
        )
        seq = "ggatatagtTCCTCCCGGcagcaaaaaacccctcaagacccgtttagaggccccaaggggttatgctagttattgctcagcggtggcagcagccaactcagcttcctttcgggctttgttagcagccggatctcagt".upper()
        kmer_map = {
            "SSS": [("s1", 0)],
            "SSR": [("s1", 1)],
            "SRQ": [("s1", 2)],
            "TGL": [("s2", 0)],
            "GLE": [("s2", 1)],
            "LEG": [("s2", 2)],
        }
        subject_map = {"s1": mock_feature1, "s2": mock_feature2}
        identity = 1.0

        features = _get_features(seq, kmer_map, subject_map, identity, True, True)

        self.assertEqual(2, len(features))
        self.assertEqual(mock_feature1.name, features[0].id)
        self.assertEqual(
            "agtTCCTCCCGGcagcaaaaaacc".upper(), features[0].extract(seq.upper())
        )
        self.assertEqual(mock_feature2.name, features[1].id)
        self.assertEqual(
            "acgggtcttgaggggttttttgct".upper(), features[1].extract(seq.upper())
        )

    def test_get_hits(self):
        """Get kmer hits against a query sequence."""

        seq = "ATGATACAGATACGAAAGTAT"
        rid = "asdf"
        id_map = {rid: SeqRecord("ATGGAT")}
        kmers = {"ATG": [(rid, 0)], "GAT": [(rid, 3)]}

        hits = _get_hits(seq, kmers, id_map, circular=False)

        self.assertEqual(2, len(hits))

    def test_filter_hits(self):
        """Filter out hits that don't show up enough."""

        seq1 = "ATGATAGACAGATAGAGATAGATGGGGAGA"
        seq2 = "GGGAC"
        subject1 = SeqRecord("ATGATAGACAGATAG", id="s1")
        subject2 = SeqRecord("AGATAGATGGGGAGA", id="s2")
        hit1 = Hit("ATGATA", 0, subject1, 0)
        hit2 = Hit("AGATAG", 0, subject1, 9)
        hit3 = Hit("AGATAG", 0, subject2, 0)
        hits = [hit1, hit2, hit3]

        filtered_hits = _filter_hits(seq1, hits)
        filtered_hits_small_seq = _filter_hits(seq2, hits)

        self.assertEqual(2, len(filtered_hits))
        self.assertIn(hit1, filtered_hits)
        self.assertIn(hit2, filtered_hits)
        self.assertEqual([], filtered_hits_small_seq)

    def test_reduce_hits(self):
        """Reduce the total number of hits when they correspond to the same feature+range."""

        subject1 = SeqRecord("ATGATAGACAGATAG", id="s1")
        subject2 = SeqRecord("AGATAGATGGGGAGA", id="s2")
        hit1 = Hit("ATGATA", 0, subject1, 0)
        hit2 = Hit("TGATAG", 1, subject1, 1)
        hit3 = Hit("GATAGA", 2, subject1, 2)
        hit4 = Hit("TGGGGC", 50, subject2, 51)

        reduced_hits = _reduce_hits([hit1, hit2, hit3, hit4])

        self.assertEqual([hit1, hit4], reduced_hits)

    def test_get_matches(self):
        """Get matches from hits against a query sequence."""

        query = "TCTCATGTGATATC"
        subject1 = SeqRecord("ACTCATGTGATAT", id="s1", name="feature1")
        hit1 = Hit("TCATG", 1, subject1, 1)
        hits = [hit1]

        matches = _get_matches(query, hits, 0.8)

        self.assertEqual(1, len(matches))
        match = matches[0]
        self.assertEqual(match.subject.id, subject1.id)
        self.assertEqual(match.query_start, 1)
        self.assertEqual(match.query_end, 12)
        self.assertEqual(match.subject_start, 1)
        self.assertEqual(match.subject_end, 12)

    def test_cull(self):
        """Remove features that have high overlap with one another."""

        seq = "ATGATAGACAGATAGAGATAGATGGGGAGA"
        feature1 = SeqFeature(FeatureLocation(1, 10, strand=1), id="1")
        feature2 = SeqFeature(FeatureLocation(2, 9, strand=1), id="2")
        feature3 = SeqFeature(FeatureLocation(12, 21, strand=1), id="3")
        feature4 = SeqFeature(FeatureLocation(14, 22, strand=1), id="4")
        identity = 0.8

        features = _cull(seq, [feature1, feature2, feature3, feature4], identity)

        self.assertEqual(2, len(features))
        self.assertTrue(all(f.id in ("1", "3") for f in features))

