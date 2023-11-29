"""Create primers to assembly SeqRecords via Gibson Assembly."""

from typing import Tuple, List, Optional, Iterable

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..designs import Combinatorial
from ..primers import Primers


MIN_HOMOLOGY = 20
"""The minimum homology (bp) between two adjacent SeqRecords."""

MAX_HOMOLOGY = 100
"""The maximum homology (bp) between two adjacent SeqRecords."""

OFFTARGET_CHECK_LEN = 10
"""The number of bp from the 3' end of a primer to use in an offtarget check.

Any other sequence spans in the template or reverse complement sequence
within 1 hamming distance of those bp is considered an off-target.
"""


def gibson_many(
    design: Iterable[List[SeqRecord]], hifi: bool = False
) -> List[Tuple[SeqRecord, List[Primers]]]:
    """Create many Gibson assemblies for each combination of SeqRecords

    Args:
        design: the list of seqrecord combinations to circularize

    Keyword Args:
        hifi: whether to use NEB's HiFi DNA ssembly (default: {False})

    Returns:
        A list of assembled Plasmids and Primers
    """

    if isinstance(design, Combinatorial):
        raise RuntimeError(
            f"Cannot create a Gibson library with Combinatorial design. See 'CombinatorialBins'."
        )

    return [gibson(rs, hifi=hifi) for rs in design]


def gibson(
    records: List[SeqRecord], hifi: bool = False
) -> Tuple[SeqRecord, List[Primers]]:
    """Create primers for records for a single Gibson Assembly.

    Create primers to mutate the records' sequences (after PCR) so they
    anneal to their neighboring records

    Args:
        records: list of records to assemble

    Keyword Args:
        hifi: whether to use HiFi DNA assembly

    Returns:
        1. assembled plasmid (SeqRecord)
        2. list of primer pairs (same length as records)
    """

    assert records

    records = [r.upper() for r in records]
    plasmid = records[0].upper()
    primers: List[Primers] = [Primers.pcr(records[0])]

    for i, f1 in enumerate(records):
        j = (i + 1) % len(records)
        f2 = records[j]

        if j != 0:
            primers.append(Primers.pcr(f2))

        # if hifi is false, mismatches is 0
        homology, homology_length, mismatch_lengths = _record_homology(
            str(f1.seq), str(f2.seq), hifi
        )

        # remove mismatches up to 10 bp if doing HiFi DNA assembly
        if hifi and mismatch_lengths:
            f1_mm_length, f2_mm_length = mismatch_lengths
            f1.seq = f1.seq[: len(f1.seq) - MAX_HOMOLOGY + f1_mm_length - 1]
            f2.seq = f2.seq[f2_mm_length:]

        if homology:  # homology already exists between records.
            plasmid += f2[homology_length:].upper()
            plasmid.id += f"|{f2.id}"
        else:
            # homology does not exist between records, introduce it to primers
            plasmid += f2.upper()
            plasmid.id += f"|{f2.id}"
            _mutate_primers(primers[i], primers[j], MIN_HOMOLOGY // 2)

    plasmid.id = "+".join(r.id for r in records if r.id != "<unknown id>")
    plasmid.seq = Seq(str(plasmid.seq.upper()), annotations={"molecule_type": "DNA"})

    # extend primers in 5' direction to avoid duplicate junctions
    _fix_duplicate_junctions(records, primers)

    return plasmid, primers


def _record_homology(
    f1: str, f2: str, hifi: bool = False
) -> Tuple[bool, int, Tuple[int, int]]:
    """Checks homology between two sequences.

    Args:
        f1: first sequence
        f2: second sequence

    Returns:
        1. true if homologous, false if not
        2. length of homologous sequence
        3. number of mismatches on either side
    """

    f1 = f1[-MAX_HOMOLOGY:]
    f2 = f2[:MAX_HOMOLOGY]
    common_substrings = _common_substring_table(f1, f2, hifi)

    # A row in the matrix corresponds to a letter in the first word
    # Eg: cat has three rows. The third row corresponds to t

    # Gets the last 11 rows of the matrix since
    # we need to account for up to 10 mismatches
    if hifi:
        last_11 = common_substrings[-11:]
        f2_length = len(last_11[0])
        # flattens last 11 rows into one list
        indices = [j for i in last_11 for j in i]
    else:  # otherwise the last row is enough
        indices = common_substrings[-1]
    while indices:
        # the length of the longest substring in the selected
        # region (indices) of matrix and the end of the substring
        end_value = max(indices)
        # 0 indicates no common substrings
        if end_value < MIN_HOMOLOGY or end_value == 0:
            break
        if hifi:
            # the index of the last row before
            # the 11 rows mentioned above
            last_start = len(common_substrings) - 11
            # list that contains end_value
            end_list = (indices.index(end_value) // f2_length) + last_start
            # index of end_value in end_list
            end_index = indices.index(end_value) % f2_length
            # index of the first letter of the substring
            start_index = end_index - end_value + 1
            start_value = common_substrings[end_list - end_value + 1][start_index]
            # the start_index must be less than 12 which accounts for up to 10
            # mismatches and the start_value must be 1 which means it is the start
            # of the substring
            if start_index < 12 and start_value == 1:
                return True, end_value, (end_list + 1, start_index - 1)
        else:
            if indices.index(end_value) == end_value:
                return True, end_value, (0, 0)
        indices[indices.index(end_value)] = 0
    return False, 0, (0, 0)


def _common_substring_table(f1: str, f2: str, hifi: bool = False) -> List[List[int]]:
    """Builds matrix of common substring lengths."""

    matrix = [[0] * (1 + len(f2)) for i in range(1 + len(f1))]
    for x in range(1, 1 + len(f1)):
        for y in range(1, 1 + len(f2)):
            # Only need lower diagonal half of matrix if not hifi
            if not hifi:
                if y > x:
                    continue
            if f1[x - 1] == f2[y - 1]:
                matrix[x][y] = matrix[x - 1][y - 1] + 1
            else:
                matrix[x][y] = 0
    return matrix


def _mutate_primers(p1: Primers, p2: Primers, homology_to_add: int):
    """Mutates primers of records if they're not homologous.

    Args:
        p1: primers of the left-most SeqRecord
        p2: primers of the right-most SeqRecord
        homology_to_add: the number of bp to add to both sides' primers from the other
    """

    # add the first homology_to_add bp of the next record (p2)
    # to the beginning of p1's reverse primer
    new_p1_rev = p2.fwd[:homology_to_add].reverse_complement() + p1.rev

    # add the last homology_to_add bp of the previous record (p1)
    # to the beginning of p2's forward primer
    new_p2_fwd = p1.rev[:homology_to_add].reverse_complement() + p2.fwd

    p1.rev = new_p1_rev
    p2.fwd = new_p2_fwd


def _fix_duplicate_junctions(records: List[SeqRecord], primers: List[Primers]):
    """Mutates primers to avoid unintentional annealing between records that shouldn't be

    Args:
        records: list of records in this assembly
        primers: list of primers to PCR each records
    """

    for i, record in enumerate(records):
        next_frag = records[(i + 1) % len(records)]
        prev_frag = records[(i - 1) % len(records)]

        # first MAX_HOMOLOGY bp of the record
        r_end_one = record.seq[:MAX_HOMOLOGY]

        # last MAX_HOMOLOGY bp of the record
        r_end_two = record.seq[-MAX_HOMOLOGY:]

        r_next_one = next_frag.seq[:MAX_HOMOLOGY]
        r_prev_two = prev_frag.seq[-MAX_HOMOLOGY:]

        ends = _ends(r_end_two, r_next_one, records)
        _mutate_junction(record, primers[i], next_frag, r_end_two, ends, True)

        ends = _ends(r_end_one, r_prev_two, records)
        _mutate_junction(record, primers[i], prev_frag, r_end_one, ends, False)


def _ends(end: Seq, neighbor_end: Seq, records: List[SeqRecord]) -> List[Seq]:
    """ Returns a list of all record ends except for
        the record end and neighbor record end passed in."""

    ends: List[Seq] = []
    for r in records:
        r_end_one = r.seq[:MAX_HOMOLOGY]
        r_end_two = r.seq[-MAX_HOMOLOGY:]
        if r_end_one != end and r_end_one != neighbor_end:
            ends.append(r_end_one)
        if r_end_two != end and r_end_two != neighbor_end:
            ends.append(r_end_two)
    return ends


def _mutate_junction(
    record: SeqRecord,
    primers: Primers,
    neighbor: SeqRecord,
    f_end: Seq,
    ends: List[Seq],
    end_of_record: bool,
) -> Optional[int]:
    """ Extends record's primer and pcr sequence if homology exists between
        the input record and any other record other than the record end
        it should anneal to. """

    i = 1
    while _has_offtarget_junction(f_end, ends, end_of_record):
        if i > max(len(primers.fwd), len(primers.rev)):
            return 0

        n_seq = neighbor.seq
        f_seq = record.seq
        if end_of_record:
            primers.rev = n_seq[:i].reverse_complement() + primers.rev
            f_seq += n_seq[:i]
            f_end = f_seq[-MAX_HOMOLOGY:]
        else:
            primers.fwd = n_seq[-i:] + primers.fwd
            f_seq = n_seq[-i:] + f_seq
            f_end = f_seq[:MAX_HOMOLOGY]
        record.seq = f_seq
        i += 1
    return None


def _has_offtarget_junction(f_end: Seq, ends: List[Seq], end_of_record: bool) -> bool:
    """Check whether f_end has offtarget junction with any other end.

    Args:
        f_end: the end of the SeqRecord we're checking
        ends: the ends of the other SeqRecords
        end_of_record: whether we're checking the 5' or 3'. 3' if end

    Returns:
        whether there's an offtarget junction
    """

    revcomp = False
    for end in ends:
        if revcomp:
            end = end.reverse_complement()
        if end_of_record:
            if _record_homology(str(f_end), str(end))[0]:
                return True
        else:
            if _record_homology(str(end), str(f_end))[0]:
                return True
        revcomp = not revcomp
    return False
