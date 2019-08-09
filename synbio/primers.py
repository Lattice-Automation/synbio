"""Primers for PCR."""

from typing import Dict, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3.bindings import designPrimers

MIN_PRIMER_LEN = 18
"""The minimum number of bp for a primer."""


class Primers:
    """Primers created by Primer3 for a SeqRecord.

    Arguments:
        fwd {Seq} -- the FWD primer
        fwd_tm {float} -- the tm of the FWD primer
        rev {Seq} -- the REV primer
        rev_tm {float} -- the tm of the REV primer
    """

    def __init__(self, fwd: Seq, fwd_tm: float, rev: Seq, rev_tm: float):
        self.fwd = fwd
        self.fwd_tm = fwd_tm
        self.rev = rev
        self.rev_tm = rev_tm

    @classmethod
    def for_record(cls, record: SeqRecord) -> "Primers":
        """Create Primers to amplify a SeqRecord.

        Arguments:
            record {SeqRecord} -- the sequence to amplify via primers

        Returns:
            Primers -- to PCR the record
        """

        p3_output = designPrimers(
            {
                "SEQUENCE_TEMPLATE": str(record.seq),
                "SEQUENCE_INCLUDED_REGION": [0, len(record.seq)],
            },
            {
                "PRIMER_PRODUCT_SIZE_RANGE": [len(record.seq), len(record.seq)],
                "PRIMER_TASK": "pick_cloning_primers",
                "PRIMER_NUM_RETURN": 1,
                "PRIMER_PICK_ANYWAY": 1,
                "PRIMER_MIN_SIZE": MIN_PRIMER_LEN,
            },
        )

        return cls.from_p3(p3_output)

    @classmethod
    def from_p3(cls, p3_output: Dict[str, str]) -> "Primers":
        """Create Primers from primer3 output (a Dict)

        Documentation for primer3 output (IO format) is available at:
        http://primer3.ut.ee/primer3web_help.htm

        Returns:
            Primers -- a primer pair from Primer3 output
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

            return Seq(seq), float(seq_tm)

        fwd_seq, fwd_tm = primer("LEFT")
        rev_seq, rev_tm = primer("RIGHT")

        return cls(fwd_seq, fwd_tm, rev_seq, rev_tm)
