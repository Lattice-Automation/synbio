"""Primers for PCR."""

from collections import defaultdict
from typing import Dict, Tuple, Union

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3.bindings import designPrimers

from .seq import mutate


MIN_PRIMER_LEN = 18
"""The minimum number of bp in a primer."""

MAX_PRIMER_LEN = 50
"""The maximum number of bp in a primer."""


class Primers:
    """Primers created by Primer3 for a SeqRecord.

    Args:
        fwd: the FWD primer
        fwd_tm: the tm of the FWD primer
        rev: the REV primer
        rev_tm: the tm of the REV primer
    """

    def __init__(self, fwd: Seq, fwd_tm: float, rev: Seq, rev_tm: float):
        self.fwd = fwd
        self.fwd_tm = fwd_tm
        self.rev = rev
        self.rev_tm = rev_tm

    @classmethod
    def pcr(
        cls,
        seq: Union[str, Seq, SeqRecord],
        fwd_padding: str = "",
        rev_padding: str = "",
    ) -> "Primers":
        """Create Primers to amplify a sequence-like object.

        Args:
            record: the sequence to amplify via primers

        Keyword Args:
            fwd_padding: Additional bp that are added to
                the 5' end of the FWD primer. Example use case:
                a restriction enzyme
            rev_padding: Additional bp that are added to
                the 5' end of the REV primer. Keep in mind that
                these are added to the 5' end of the rev primer,
                so they're the reverse complement bp of the
                template sequence

        Returns:
            A Primers object to amplify the SeqRecord
        """

        template = _get_seq(seq)

        p3_output = designPrimers(
            {
                "SEQUENCE_TEMPLATE": template,
                "SEQUENCE_INCLUDED_REGION": [0, len(template)],
            },
            {
                "PRIMER_PRODUCT_SIZE_RANGE": [len(template), len(template)],
                "PRIMER_TASK": "pick_cloning_primers",
                "PRIMER_NUM_RETURN": 1,
                "PRIMER_PICK_ANYWAY": 1,
                "PRIMER_MIN_SIZE": MIN_PRIMER_LEN,
            },
        )

        primers = cls.from_p3(p3_output)

        # add additional bp to the primers
        # the tm of the primers aren't updated: they refer just
        # to the binding bp of the primers to the template sequence
        if fwd_padding:
            primers.fwd = fwd_padding + primers.fwd
        if rev_padding:
            primers.rev = rev_padding + primers.rev

        return primers

    @classmethod
    def from_p3(cls, p3_output: Dict[str, str]) -> "Primers":
        """Create Primers from primer3 output (a Dict)

        Documentation for primer3 output (IO format) is available at:
        http://primer3.ut.ee/primer3web_help.htm

        Returns:
            a primer pair from Primer3 output
        """

        assert p3_output

        if "PRIMER_WARNING" in p3_output:
            raise RuntimeWarning(f"Primer3 warning: {p3_output['PRIMER_WARNING']}")

        if "PRIMER_ERROR" in p3_output:
            raise RuntimeError(f"Primer3 error: {p3_output['PRIMER_ERROR']}")

        if p3_output["PRIMER_PAIR_NUM_RETURNED"] == 0:
            raise RuntimeError(f"Primer3: no primer pair created")

        # based on my REPP version of same thing:
        # https://github.com/Lattice-Automation/repp/blob/master/internal/repp/primer3.go
        def primer(side: str) -> Tuple[Seq, float]:
            """Get the sequence and tm of the FWD or REV primer."""

            assert side in ("LEFT", "RIGHT")
            seq = p3_output[f"PRIMER_{side}_0_SEQUENCE"]
            seq_tm = p3_output[f"PRIMER_{side}_0_TM"]
            return Seq(seq, alphabet=IUPACUnambiguousDNA()), float(seq_tm)

        fwd_seq, fwd_tm = primer("LEFT")
        rev_seq, rev_tm = primer("RIGHT")

        return cls(fwd_seq, fwd_tm, rev_seq, rev_tm)

    def __eq__(self, other) -> bool:
        """Primers equality checking.

        Primers are the same if their sequences are same in FWD/REV direction.

        Args:
            other: Other Primers object

        Returns:
            Whether the two Primers objects are the same
        """

        if not isinstance(other, Primers):
            return False

        return self.fwd == other.fwd and self.rev == other.rev

    def specific(
        self,
        seq: Union[str, Seq, SeqRecord],
        check_len: int = 12,
        edit_distance: int = 1,
    ) -> bool:
        """Return whether the primers will produce a single, unambiguous PCR product.

        This returns false if either:
            1. The primers don't anneal to the template sequence and won't amplify
            2. The primers bind in multiple locations and may create multiple
                PCR product via ectopic binding

        Ectopic/off-target binding is checked by looking at the last 10 bp of the primer.
        If there is a binding site anywhere in the template sequence with only
        one bp of difference between either of the primers' 3' end and the template
        sequence, we classify that as being sufficient for an ectopic binding site.
        This assumption is based on the following paper, which looked at the effects
        of single bp mutations in qPCR rates. Single bp changes had a deleterious effect:

        "Quantitative effects of position and type of single mismatch on
        single base primer extension". Journal of Microbiological Methods, 2009.

        Args:
            seq: A template sequence like object (str, Seq or SeqRecord)

        Keyword Args:
            check_len: The number of bp from 3' end of primers to check for
                binding in the template sequence. Larger number bp mean this
                is less likely to find ectopic binding sites. Conversely,
                smaller check_len's mean this is more likely to find
                off-target, ectopic sites and not count the primers as specific
            edit_distance: The number of edits to check for within the
                ends of each primer when looking for ectopic binding sites

        Returns:
            Whether the primers will create a single PCR product
        """

        template = _get_seq(seq).upper()
        template_rc = str(Seq(template).reverse_complement())

        kmer_count: Dict[str, int] = defaultdict(int)
        for i in range(len(template) - check_len + 1):
            kmer = template[i : i + check_len]
            kmer_rc = template_rc[i : i + check_len]

            kmer_count[kmer] += 1
            kmer_count[kmer_rc] += 1

        end_fwd = self.fwd[-check_len:]
        end_rev = self.rev[-check_len:]

        # return False if the fwd and rev primer don't bind
        # to opposite strands
        if not (
            (end_fwd in template and end_rev in template_rc)
            or (end_fwd in template_rc and end_rev in template)
        ):
            return False

        # return False if the fwd and or rev primer bind more than once
        if kmer_count[end_fwd] > 1 or kmer_count[end_rev] > 1:
            return False

        fwd_mutants = mutate(end_fwd, edit_distance=edit_distance)
        rev_mutants = mutate(end_rev, edit_distance=edit_distance)
        mutants = fwd_mutants.union(rev_mutants)

        # remove the original binding sites, those were checked above
        mutants.remove(end_fwd)
        mutants.remove(end_rev)

        # return False if any of the off-by one sequences will bind
        # more than once
        if any(kmer_count[m] > 0 for m in mutants):
            return False

        return True


def _get_seq(seq: Union[str, Seq, SeqRecord]) -> str:
    """Get the sequence as a str from seq-like object."""

    if isinstance(seq, str):
        return seq
    if isinstance(seq, Seq):
        return str(seq)
    if isinstance(seq, SeqRecord):
        return str(seq.seq)
    raise TypeError(seq)
