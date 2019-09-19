"""Test clustering functions."""

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from synbio.features.cluster import _parse_row, _consensus_name


class TestCluster(unittest.TestCase):
    """Test feature database creation."""

    def test_parse_row(self):
        """Parse a row of iGEM XML to a SeqRecord."""

        row = """<row>
                <field name="part_id">2557</field>
                <field name="ok">0</field>
                <field name="part_name">BBa_S01288</field>
                <field name="short_desc">Intermediate part from assembly 236</field>
                <field name="description" xsi:nil="true" />
                <field name="part_type">Intermediate</field>
                <field name="author">Randy Rettberg</field>
                <field name="owning_group_id">7</field>
                <field name="status">Deleted</field>
                <field name="dominant">0</field>
                <field name="informational">0</field>
                <field name="discontinued">1</field>
                <field name="part_status"></field>
                <field name="sample_status">Discontinued</field>
                <field name="p_status_cache"></field>
                <field name="s_status_cache"></field>
                <field name="creation_date">2003-12-03</field>
                <field name="m_datetime">2018-11-06 15:23:57</field>
                <field name="m_user_id">0</field>
                <field name="uses">0</field>
                <field name="doc_size">533</field>
                <field name="works"></field>
                <field name="favorite">0</field>
                <field name="specified_u_list">_149_156_603_145_193_147_161_603_145_</field>
                <field name="deep_u_list">_149_156_603_145_193_147_161_603_145_</field>
                <field name="deep_count">9</field>
                <field name="ps_string" xsi:nil="true" />
                <field name="scars"></field>
                <field name="default_scars"></field>
                <field name="owner_id">24</field>
                <field name="group_u_list">_1_</field>
                <field name="has_barcode">0</field>
                <field name="notes" xsi:nil="true" />
                <field name="source"></field>
                <field name="nickname"></field>
                <field name="categories">//classic/intermediate/uncategorized</field>
                <field name="sequence">tcacacaggaaagtactagatgagcacaaaaaagaaaccattaacacaagagcagcttgaggacgcacgtcgccttaaagcaatttatgaaaaaaagaaaaatgaacttggcttatcccaggaatctgtcgcagacaagatggggatggggcagtcaggcgttggtgctttatttaatggcatcaatgcattaaatgcttataacgccgcattgcttgcaaaaattctcaaagttagcgttgaagaatttagcccttcaatcgccagagaaatctacgagatgtatgaagcggttagtatgcagccgtcacttagaagtgagtatgagtaccctgttttttctcatgttcaggcagggatgttctcacctgagcttagaacctttaccaaaggtgatgcggagagatgggtaagcacaaccaaaaaagccagtgattctgcattctggcttgaggttgaaggtaattccatgaccgcaccaacaggctccaagccaagctttcctgacggaatgttaattctcgttgaccctgagcaggctgttgagccaggtgatttctgcatagccagacttgggggtgatgagtttaccttcaagaaactgatcagggatagcggtcaggtgtttttacaaccactaaacccacagtacccaatgatcccatgcaatgagagttgttccgttgtggggaaagttatcgctagtcagtggcctgaagagacgtttggcgctgcaaacgacgaaaactacgctttagtagcttaataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatatactagagacctgtacgatcctacaggtgcttatgttaagtaattgtattcccagcgatacaatagtgtgacaaaaatccaatttattagaatcaaatgtcaatccattaccgttttaatgatatataacacgcaaaacttgcgacaaacaataggtaatactagagattaaagaggagaaatactagatgaaaaacataaatgccgacgacacatacagaataattaataaaattaaagcttgtagaagcaataatgatattaatcaatgcttatctgatatgactaaaatggtacattgtgaatattatttactcgcgatcatttatcctcattctatggttaaatctgatatttcaatcctagataattaccctaaaaaatggaggcaatattatgatgacgctaatttaataaaatatgatcctatagtagattattctaactccaatcattcaccaattaattggaatatatttgaaaacaatgctgtaaataaaaaatctccaaatgtaattaaagaagcgaaaacatcaggtcttatcactgggtttagtttccctattcatacggctaacaatggcttcggaatgcttagttttgcacattcagaaaaagacaactatatagatagtttatttttacatgcgtgtatgaacataccattaattgttccttctctagttgataattatcgaaaaataaatatagcaaataataaatcaaacaacgatttaaccaaaagagaaaaagaatgtttagcgtgggcatgcgaaggaaaaagctcttgggatatttcaaaaatattaggttgcagtgagcgtactgtcactttccatttaaccaatgcgcaaatgaaactcaatacaacaaaccgctgccaaagtatttctaaagcaattttaacaggagcaattgattgcccatactttaaaaattaataacactgatagtgctagtgtagatcactactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata</field>
                <field name="sequence_sha1">?^XT]^P^?f ?]M</field>
                <field name="sequence_update">5</field>
                <field name="seq_edit_cache">&lt;script src='http://parts.igem.org/cgi/partsdb/seq_edit/all.js'&gt;&lt;/script&gt; &lt;DIV id='sequencePaneDiv' style='clear:both'&gt; &lt;INPUT type='hidden' id='new_dna_format' name='new_dna_format' value='' /&gt; &lt;INPUT type='hidden' id='selection_start' name='selection_start' value='14' /&gt; &lt;INPUT type='hidden' id='selection_end'   name='selection_end'   value='0' /&gt;&lt;/DIV&gt; &lt;script&gt; var sequence = new String('ttaggtattgactgtactatcagttccgtcataatatgaaccataagttcaccac');        var seqFeatures = new Array( ['brick',1,55,'R0050', 0], ['reg',8,13,'-35', 0], ['reg',30,35,'-10', 0], ['rbs',13,27,'OR2', 0], ['rbs',37,51,'OR1', 0], ['rarrow_p',43,43,'putative', 0]); var subParts = null; var Format = '_ruler_'; var PrimaryPartName = 'BBa_R0050'; var PrimaryPartID = '188'; var Selection_Start = 0; var Selection_End = 0; showSeqFeatures(false);  &lt;/script&gt;&lt;div style='position:relative;clear:both;width:100%'&gt;&lt;div style=''&gt;
&lt;STYLE type='text/css'&gt;
.compatibility_div ul,
.compatibility_div li {
display: inline;
}
.compatibility_div li {
position: relative;
padding-top: 2px;
padding-left:4px;
padding-right:3px;
margin-right:2px;
margin-bottom: 5px;
        """

        record = _parse_row(row)

        self.assertEqual("BBa_S01288", record.id)
        self.assertEqual(
            "tcacacaggaaagtactagatgagcacaaaaaagaaaccattaacacaagagcagcttgaggacgcacgtcgccttaaagcaatttatgaaaaaaagaaaaatgaacttggcttatcccaggaatctgtcgcagacaagatggggatggggcagtcaggcgttggtgctttatttaatggcatcaatgcattaaatgcttataacgccgcattgcttgcaaaaattctcaaagttagcgttgaagaatttagcccttcaatcgccagagaaatctacgagatgtatgaagcggttagtatgcagccgtcacttagaagtgagtatgagtaccctgttttttctcatgttcaggcagggatgttctcacctgagcttagaacctttaccaaaggtgatgcggagagatgggtaagcacaaccaaaaaagccagtgattctgcattctggcttgaggttgaaggtaattccatgaccgcaccaacaggctccaagccaagctttcctgacggaatgttaattctcgttgaccctgagcaggctgttgagccaggtgatttctgcatagccagacttgggggtgatgagtttaccttcaagaaactgatcagggatagcggtcaggtgtttttacaaccactaaacccacagtacccaatgatcccatgcaatgagagttgttccgttgtggggaaagttatcgctagtcagtggcctgaagagacgtttggcgctgcaaacgacgaaaactacgctttagtagcttaataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatatactagagacctgtacgatcctacaggtgcttatgttaagtaattgtattcccagcgatacaatagtgtgacaaaaatccaatttattagaatcaaatgtcaatccattaccgttttaatgatatataacacgcaaaacttgcgacaaacaataggtaatactagagattaaagaggagaaatactagatgaaaaacataaatgccgacgacacatacagaataattaataaaattaaagcttgtagaagcaataatgatattaatcaatgcttatctgatatgactaaaatggtacattgtgaatattatttactcgcgatcatttatcctcattctatggttaaatctgatatttcaatcctagataattaccctaaaaaatggaggcaatattatgatgacgctaatttaataaaatatgatcctatagtagattattctaactccaatcattcaccaattaattggaatatatttgaaaacaatgctgtaaataaaaaatctccaaatgtaattaaagaagcgaaaacatcaggtcttatcactgggtttagtttccctattcatacggctaacaatggcttcggaatgcttagttttgcacattcagaaaaagacaactatatagatagtttatttttacatgcgtgtatgaacataccattaattgttccttctctagttgataattatcgaaaaataaatatagcaaataataaatcaaacaacgatttaaccaaaagagaaaaagaatgtttagcgtgggcatgcgaaggaaaaagctcttgggatatttcaaaaatattaggttgcagtgagcgtactgtcactttccatttaaccaatgcgcaaatgaaactcaatacaacaaaccgctgccaaagtatttctaaagcaattttaacaggagcaattgattgcccatactttaaaaattaataacactgatagtgctagtgtagatcactactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata",
            str(record.seq),
        )
        self.assertEqual(
            len(record.features), 3
        )  # only three are greater than the default word size

    def test_consensus_name(self):
        """Find a consensus name for a list of SeqRecords."""

        frt1 = SeqRecord(
            Seq(""),
            id="BBa_J61020",
            annotations={
                "short_desc": "[FRT]",
                "description": "Site for recombination by flp recombinase.  Genes flanked by FRT sites (oriented in the same direction) in the genome can be excised with the introduction of flp recombinase-expressing helper plasmid pCP20.  Note that the original FRT sequence contained an XbaI site that has been removed with a point mutation for compatibility with standard assembly.  See Datsenko and Wanner for details of its use in markerless knockouts and knockins.",
                "categories": "//function/recombination/flp",
                "nickname": "FRT",
            },
        )
        frt2 = SeqRecord(
            Seq(""),
            id="BBa_J72001",
            annotations={
                "short_desc": "{FRT} recombination site for flp recombinase in BBb",
                "description": "Later",
                "categories": "//function/recombination/flp",
                "nickname": "FRT",
            },
        )
        frt3 = SeqRecord(
            Seq(""),
            id="BBa_K337019",
            annotations={
                "short_desc": "FRT site with spacing",
                "description": "FRT sites can be used to stably integrate DNA parts into a host cell genome via the Flip-in System.",
                "categories": "",
                "nickname": "",
            },
        )

        feature_name = _consensus_name([frt1, frt2, frt3])

        self.assertEqual("FRT", feature_name)

