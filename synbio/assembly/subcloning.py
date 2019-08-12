"""Subcloning for digestion and ligation of SeqRecords together.

see: https://en.wikipedia.org/wiki/Subcloning
"""

import logging
from typing import Dict, List, Set, Tuple, Iterable, Union

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Restriction import RestrictionBatch, BsaI, BpiI
from Bio.Restriction.Restriction import RestrictionType
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import networkx as nx
from networkx.algorithms.cycles import simple_cycles
from networkx.exception import NetworkXNoCycle

from ..designs import Plasmid


CATALYZE_CACHE: Dict[str, List[Tuple[str, SeqRecord, str]]] = {}
"""Store the catalyze results of each SeqRecord. Avoid lots of string searches."""


def goldengate(
    design: Union[List[SeqRecord], Iterable[List[SeqRecord]]],
    resistance: str = "",
    min_count: int = -1,
) -> List[Tuple[SeqRecord, List[SeqRecord]]]:
    """Simulate a digestion and ligation using BsaI and BpiI.

    Accepts lists of SeqRecords that may combine and returns the lists
    of SeqRecord design that may circularize into new vectors.

    Arguments:
        design {Iterable[List[SeqRecord]]} -- possible combinations of fragments

    Keyword Arguments:
        resistance {str} -- the feature to filter assemblies on (default: {""})
        min_count {int} -- minimum number of SeqRecords for an assembly to be considered

    Returns:
        List[List[SeqRecord]] -- list of combinations of digested SeqRecords
    """

    return subclone(design, [BsaI, BpiI], resistance, min_count)


def subclone(
    design: Union[List[SeqRecord], Iterable[List[SeqRecord]]],
    enzymes: List[RestrictionType],
    resistance: str = "",
    min_count: int = -1,
) -> List[Tuple[SeqRecord, List[SeqRecord]]]:
    """Return lists of combinations of fragments that will circularize. One to many.

    Given a list of SeqRecord combinations return combinations that will circularize
    into new plasmids. Simulate a digest of each fragment with all enzymes,
    check overhangs with other SeqRecords in the same well, and check what may
    form during ligation.

    Arguments:
        design {List[List[SeqRecord]]} -- list of possible SeqRecord combinations
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with

    Keyword Arguments:
        resistance {str} -- the resistance to filter assemblies against (ex: 'KanR')
        min_count {int} -- mininum number of SeqRecords for an assembly to be considered

    Returns:
        List[Tuple[SeqRecord, List[SeqRecord]]] -- list of tuples with:
            1. composite plasmids expected after subcloning
            2. list of linearized SeqRecords that went into the design
    """

    assert design, f"Cannot subclone. List of SeqRecords required."

    design = list(design)
    if not isinstance(design[0], list):
        design = [design]  # make a list of lists

    valid_assemblies: List[Tuple[SeqRecord, List[SeqRecord]]] = []
    for combination in design:
        for assembly in _valid_assemblies(combination, enzymes, resistance, min_count):
            plasmid = SeqRecord(Seq("", IUPACUnambiguousDNA()))
            for fragment in assembly:
                plasmid += fragment.upper()
            plasmid.id = "|".join(f.id for f in assembly if f.id != "<unknown id>")

            valid_assemblies.append((plasmid, assembly))

    if type(design) == Plasmid and len(valid_assemblies) > 1:
        warn = RuntimeWarning(
            f"Plasmid design specified using only first valid assembly"
        )
        logging.warning(warn)
        valid_assemblies = valid_assemblies[:1]

    return valid_assemblies


def _valid_assemblies(
    record_set: List[SeqRecord],
    enzymes: List[RestrictionType],
    resistance: str,
    min_count: int,
) -> List[List[SeqRecord]]:
    """Parse a single list of SeqRecords to find all circularizable plasmids.

    Turn each SeqRecord's post-digest seqs into a graph where the nodes are
    the overhangs and the edges are the linear fragments
    post-digest/catalyzing with BsaI/BpiI.

    Arguments:
        record_set {List[SeqRecord]} -- single record set that might circularize
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with
        resistance {str} -- the resistance to filter assemblies
        min_count {int} -- mininum number of SeqRecords for an assembly to be considered

    Returns:
        List[List[SeqRecord]] -- a list of SeqRecords/fragment combos to mix for MoClo
    """

    graph = nx.DiGraph()

    input_seqs: Set[str] = set()  # stored list of input seqs (not new combinations)
    for record in record_set:
        input_seqs.add(str(record.seq + record.seq))
        for left, frag, right in _catalyze(record, enzymes):
            graph.add_node(left)
            graph.add_node(right)
            graph.add_edge(left, right, frag=frag)

    try:  # find all circularizable cycles
        cycles = simple_cycles(graph)
    except NetworkXNoCycle:
        return []

    # get the fragments, enzymes back out of the cycle
    assemblies: List[List[SeqRecord]] = []
    for cycle in cycles:
        if min_count > 0 and len(cycle) < min_count:
            continue

        fragments = []
        for i, overhang in enumerate(cycle):
            next_overhang = cycle[(i + 1) % len(cycle)]
            fragments.append(graph.edges[overhang, next_overhang]["frag"])

        # make sure it's not just a re-ligation of insert + backbone
        new_plasmid = "".join([str(f.seq) for f in fragments])
        if any(new_plasmid in seq for seq in input_seqs):
            continue

        # filter for plasmids that have the resistance feature
        if resistance:
            has_resistance = False
            for fragment in fragments:
                if _has_resistance(fragment, resistance):
                    has_resistance = True
                    break

            if not has_resistance:
                continue

        # try and re-order the fragments to match the input order
        record_set_ids = [r.id for r in record_set]
        fragment_indexes = [record_set_ids.index(f.id) for f in fragments]
        fragment_min_index = min(fragment_indexes)
        fragment_first = fragment_indexes.index(fragment_min_index)
        fragments = fragments[fragment_first:] + fragments[:fragment_first]

        assemblies.append(fragments)
    return assemblies


def _catalyze(
    record: SeqRecord, enzymes: List[RestrictionType]
) -> List[Tuple[str, SeqRecord, str]]:
    """Catalyze a SeqRecord and return all post-digest SeqRecords with overhangs.

    Overhangs are returned as the overhang plus the position of the cut
    in the 5' end (^) and 3' end (_). So a 5' overhang may be:
    ^AAAA_. But a 3' overhang may be: _AAAA^.

    Arguments:
        record {SeqRecord} -- the SeqRecord to digest with enzymes
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with

    Returns:
        List[Tuple[str, SeqRecord, str]] --
            left overhang, cut fragment, right overhang
    """

    record_id = _record_id(record)
    if record_id in CATALYZE_CACHE:
        return CATALYZE_CACHE[record_id]

    batch = RestrictionBatch(enzymes)
    batch_sites = batch.search(record.seq, linear=False)

    # order all cuts with enzymes based on index
    cuts_seen: Set[int] = set()
    enzyme_cuts: List[Tuple[RestrictionType, int]] = []
    for enzyme, cuts in batch_sites.items():
        for cut in cuts:
            if cut in cuts_seen:
                continue
            cuts_seen.add(cut)
            enzyme_cuts.append((enzyme, cut - 1))  # revert to 0-based
    enzyme_cuts = sorted(enzyme_cuts, key=lambda x: x[1])

    # list of left/right overhangs for each fragment
    frag_w_overhangs: List[Tuple[str, SeqRecord, str]] = []
    for i, (enzyme, cut) in enumerate(enzyme_cuts):
        next_enzyme, next_cut = enzyme_cuts[(i + 1) % len(enzyme_cuts)]
        frag = record[cut:next_cut]

        # calculate the junction on the left-hand side
        len_ovhg = len(enzyme.ovhgseq)
        left_overhang = (
            record[cut - len_ovhg : cut] + "^"
            if enzyme.is_3overhang()
            else "^" + record[cut : cut + len_ovhg]
        )

        # calculate the junction on the right-hand side
        right_overhang = (
            record[next_cut - len_ovhg : next_cut] + "^"
            if next_enzyme.is_3overhang()
            else "^" + record[next_cut : next_cut + len_ovhg]
        )

        if next_cut < cut:  # wraps around the zero-index
            frag = (record + record)[cut : next_cut + len(record)]

        frag_w_overhangs.append((str(left_overhang.seq), frag, str(right_overhang.seq)))

    CATALYZE_CACHE[record_id] = frag_w_overhangs  # store for future look-ups

    return frag_w_overhangs


def _has_resistance(record: SeqRecord, resistance: str) -> bool:
    """Return whether any of a record's features/qualifiers match the resistance specified.

    Arguments:
        record {SeqRecord} -- the record being checked for resistance
        resistance {str} -- the resistance to filter for

    Returns:
        bool -- whether the record has any features or qualifiers with specified resistance
    """

    for feature in record.features:
        if feature.id is resistance:
            return True
        for _, value in feature.qualifiers.items():
            if resistance in value:
                return True
    return False


def _record_id(record: SeqRecord) -> str:
    """Get the record id, unique, for a seqRecord."""

    return record.id if record.id != "<unknown id>" else str(record.seq)
