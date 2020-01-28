"""Primers for PCR."""

from typing import Union

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primers import primers


class Primers:
    """Primers created by Primer3 for a SeqRecord.

    Attributes:
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
            record: the sequence-like object to amplify via primers

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

        fwd, rev = primers(template, add_fwd=fwd_padding, add_rev=rev_padding)

        return Primers(
            Seq(fwd.seq, alphabet=IUPACUnambiguousDNA()),
            fwd.tm,
            Seq(rev.seq, alphabet=IUPACUnambiguousDNA()),
            rev.tm,
        )

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


def _get_seq(seq: Union[str, Seq, SeqRecord]) -> str:
    """Get the sequence as a str from seq-like object."""

    if isinstance(seq, str):
        return seq
    if isinstance(seq, Seq):
        return str(seq)
    if isinstance(seq, SeqRecord):
        return str(seq.seq)
    raise TypeError(
        f"invalid non-sequence-like argument: '{seq}'. "
        + "Expecting str, Bio.Seq, or Bio.SeqRecord"
    )

