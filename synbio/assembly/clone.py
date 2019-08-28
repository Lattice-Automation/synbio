"""Cloning for digestion and ligation of fragments."""

from collections import defaultdict
from typing import Dict, List, Set, Tuple, Iterable, Optional

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Restriction import RestrictionBatch, BsaI, BpiI
from Bio.Restriction.Restriction import RestrictionType
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import networkx as nx
from networkx.algorithms.cycles import simple_cycles
from networkx.exception import NetworkXNoCycle

from ..containers import content_id
from ..designs import CombinatorialBins


CATALYZE_CACHE: Dict[str, List[Tuple[str, SeqRecord, str]]] = {}
"""Store the catalyze results of each SeqRecord. Avoid lots of string searches."""


def clone(
    record_set: Iterable[SeqRecord],
    enzymes: List[RestrictionType],
    include: List[str] = None,
    min_count: int = -1,
    linear: bool = True,
):
    """Simulate a digestion/ligation with a list of enzymes to get expected plasmids.

    Each SeqRecord is digested with all enzymes and a graph is created with
    nodes from each digested records overhangs. Circularizable assemblies/plasmids
    are returned as a list of ligated plasmids.

    Note: this does not return, explicitly, which enzymes and SeqRecords were
    used to create each plasmid. But that information is implicitly available
    in each SeqRecord's id.

    Args:
        record_set: SeqRecords to combine into a plasmid
        enzymes: Enzymes to digest each SeqRecord with

    Keyword Args:
        include: feature names to filter plasmids against (default: {None})
        min_count: mininum # of SeqRecords in each output plasmid (default: {-1})
        linear: Whether the individual SeqRecords are assumed to be linear
    """

    cloned_plasmids: List[SeqRecord] = []
    for plasmids, _ in clone_combinatorial(
        list(record_set), enzymes, include=include, min_count=min_count, linear=linear
    ):
        cloned_plasmids.extend(plasmids)
    return cloned_plasmids


def goldengate(
    record_set: Iterable[List[SeqRecord]],
    enzymes: List[RestrictionType] = [BsaI, BpiI],
    include: List[str] = None,
    min_count: int = -1,
    linear: bool = True,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Simulate a digestion and ligation using BsaI and BpiI.

    Accepts lists of SeqRecords that may combine and returns the lists
    of SeqRecords that may circularize into new vectors.

    Args:
        record_set: possible combinations of fragments

    Keyword Args:
        include: the feature to filter assemblies on (default: {""})
        min_count: minimum number of SeqRecords for an assembly to be considered
        linear: Whether the individual SeqRecords are assumed to be linear

    Returns:
        A list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    return clone_many_combinatorial(
        record_set, enzymes, include=include, min_count=min_count, linear=linear
    )


def clone_many_combinatorial(
    design: Iterable[List[SeqRecord]],
    enzymes: List[RestrictionType],
    include: List[str] = None,
    min_count: int = -1,
    linear: bool = True,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Parse a single list of SeqRecords to find all circularizable plasmids.

    Turn each SeqRecord's post-digest seqs into a graph where the nodes are
    the overhangs and the edges are the linear fragments post-digest/catalyzing
    with BsaI/BpiI.

    Args:
        record_set: single record set that might circularize
        enzymes: list of enzymes to digest the input records with

    Keyword Args:
        include: List of strings to filter assemblies against
        min_count: The mininum number of SeqRecords for an assembly to be considered
        linear: Whether the individual SeqRecords are assumed to be linear

    Returns:
        List[Tuple[List[SeqRecord], List[SeqRecord]]] -- list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    seen_fragment_ids: Set[str] = set()
    all_plasmids_and_fragments: List[Tuple[List[SeqRecord], List[SeqRecord]]] = []
    for record_set in design:
        for plasmids, fragments in clone_combinatorial(
            record_set, enzymes, include=include, min_count=min_count, linear=linear
        ):

            # we don't want to re-use the fragment combination more than once
            fragment_ids = _hash_fragments(fragments)
            if fragment_ids in seen_fragment_ids:
                continue
            seen_fragment_ids.add(fragment_ids)

            all_plasmids_and_fragments.append((plasmids, fragments))
    return all_plasmids_and_fragments


def clone_combinatorial(
    record_set: List[SeqRecord],
    enzymes: List[RestrictionType],
    include: List[str] = None,
    min_count: int = -1,
    linear: bool = True,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Parse a single list of SeqRecords to find all circularizable plasmids.

    Turn each SeqRecord's post-digest seqs into a graph where the nodes are
    the overhangs and the edges are the linear fragments
    post-digest/catalyzing with BsaI/BpiI.

    Args:
        record_set: single record set that might circularize
        enzymes: list of enzymes to digest the input records with

    Keyword Args:
        include: the include to filter assemblies
        min_count: mininum number of SeqRecords for an assembly to be considered
        linear: Whether the individual SeqRecords are assumed to be linear

    Returns:
        A list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    graph = nx.MultiDiGraph()

    seen_seqs: Set[str] = set()  # stored list of input seqs (not new combinations)
    for record in record_set:
        seen_seqs.add(str(record.seq + record.seq).upper())
        seen_seqs.add(str((record.seq + record.seq).reverse_complement().upper()))

        for left, frag, right in _catalyze(record, enzymes, linear):
            graph.add_node(left)
            graph.add_node(right)
            graph.add_edge(left, right, frag=frag)

    try:  # find all circularizable cycles
        cycles = simple_cycles(graph)
    except NetworkXNoCycle:
        return []

    # get the fragments, enzymes back out of the cycle
    ids_to_fragments: Dict[str, List[SeqRecord]] = defaultdict(list)
    ids_to_plasmids: Dict[str, List[SeqRecord]] = defaultdict(list)
    for cycle in cycles:
        # filter for the minimum number of SeqRecords
        if min_count > 0 and len(cycle) < min_count:
            continue

        combinations = CombinatorialBins()
        for i, overhang in enumerate(cycle):
            next_overhang = cycle[(i + 1) % len(cycle)]
            record_bin = []
            for out_edge in graph.out_edges(keys=True):
                src, dest, index = out_edge
                if src != overhang or dest != next_overhang:
                    continue
                record_bin.append(graph.edges[src, dest, index]["frag"])
            combinations.append(record_bin)

        for fragments in combinations:
            # create the composite plasmid
            plasmid = SeqRecord(Seq("", IUPACUnambiguousDNA()))
            for fragment in fragments:
                plasmid += fragment.upper()

            # make sure it's not just a re-ligation of insert + backbone
            plasmid_seq = str(plasmid.seq)
            if any(plasmid_seq in seq for seq in seen_seqs):
                continue

            # filter for plasmids that have an 'include' feature
            if not _has_features(plasmid, include):
                continue

            # re-order the fragments to try and match the input order
            fragments = _reorder_fragments(record_set, fragments)

            seen_seqs.add(str(plasmid.seq + plasmid.seq))
            seen_seqs.add(str((plasmid.seq + plasmid.seq).reverse_complement()))

            # make a unique id for the fragments
            fragments_id = _hash_fragments(fragments)
            ids_to_fragments[fragments_id] = fragments
            ids_to_plasmids[fragments_id].append(plasmid)

    plasmids_and_fragments: List[Tuple[List[SeqRecord], List[SeqRecord]]] = []
    for ids, fragments in ids_to_fragments.items():
        plasmids = ids_to_plasmids[ids]
        for i, plasmid in enumerate(plasmids):
            plasmid.id = "+".join(f.id for f in fragments if f.id != "<unknown id>")
            plasmid.description = f"cloned from {', '.join(str(e) for e in enzymes)}"

            if len(plasmids) > 1:
                plasmid.id += f"({i + 1})"
        plasmids_and_fragments.append((plasmids, fragments))
    return plasmids_and_fragments


def _reorder_fragments(
    input_set: List[SeqRecord], output_set: List[SeqRecord]
) -> List[SeqRecord]:
    """Sort the output set of SeqRecords to try and match input order.

    If the first SeqRecord of the input_set is the third SeqRecord of the
    output_set, reorder the output_set so its first SeqRecord is the first again.

    This works because the SeqRecords are going to anneal to one another in a plasmid.

    Args:
        input_set: SeqRecords set to clone_combinatorial
        output_set: SeqRecords after circularization

    Returns:
        List[SeqRecord] -- re-ordered SeqRecords
    """

    fragment_ids = [content_id(r) for r in input_set]
    fragment_indexes = [fragment_ids.index(f.id) for f in output_set]
    fragment_min_index = min(fragment_indexes)
    fragment_first = fragment_indexes.index(fragment_min_index)
    return output_set[fragment_first:] + output_set[:fragment_first]


def _hash_fragments(record_set: List[SeqRecord]) -> str:
    """Create a unique ID for a list of records

    Args:
        record_set: set of Records to make unique ID for

    Returns:
        A unique ID concatenating records IDs/Seqs
    """

    fragment_ids = [content_id(f) for f in record_set]
    fragment_ids = sorted(fragment_ids)
    return "".join(fragment_ids)


def _catalyze(
    record: SeqRecord, enzymes: List[RestrictionType], linear=True
) -> List[Tuple[str, SeqRecord, str]]:
    """Catalyze a SeqRecord and return all post-digest SeqRecords with overhangs.

    Overhangs are returned as the overhang plus the position of the cut
    in the 5' end (^) and 3' end (_). So a 5' overhang may be:
    ^AAAA_. But a 3' overhang may be: _AAAA^.

    Args:
        record: The SeqRecord to digest with enzymes
        enzymes: List of enzymes to digest the input records with

    Keyword Args:
        linear: Whether the record to catalyze is linear or circular

    Returns:
        Tuple with: (left overhang, cut fragment, right overhang)
    """

    record = record.upper()
    batch = RestrictionBatch(enzymes)
    batch_sites = batch.search(record.seq, linear=linear)

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
        if i == len(enzyme_cuts) - 1 and linear:
            continue

        next_enzyme, next_cut = enzyme_cuts[(i + 1) % len(enzyme_cuts)]

        enzyme_len = len(enzyme.ovhgseq)
        next_enzyme_len = len(next_enzyme.ovhgseq)

        # shift cuts left for 3overhang enzymes
        if enzyme.is_3overhang():
            cut -= enzyme_len
        if next_enzyme.is_3overhang():
            next_cut -= next_enzyme_len

        cut_rc = cut if enzyme.is_3overhang() else cut + enzyme_len
        next_cut_rc = (
            next_cut if next_enzyme.is_3overhang() else next_cut + next_enzyme_len
        )

        # find the cutsite sequences
        left = record[cut : cut + enzyme_len]
        right = record[next_cut : next_cut + next_enzyme_len]
        left_rc = right.reverse_complement()
        right_rc = left.reverse_complement()

        left = str(left.seq)
        right = str(right.seq)
        left_rc = str(left_rc.seq)
        right_rc = str(right_rc.seq)

        if enzyme.is_3overhang():
            left += "^"
            right_rc += "^"
        else:
            left = "^" + left
            right_rc = "^" + right_rc

        if next_enzyme.is_3overhang():
            right += "^"
            left_rc += "^"
        else:
            right = "^" + right
            left_rc = "^" + left_rc

        # shift cuts right again for 3overhang enzymes
        if enzyme.is_3overhang():
            cut += enzyme_len
        if next_enzyme.is_3overhang():
            next_cut += next_enzyme_len

        frag = record[cut:next_cut]
        frag_rc = record[cut_rc:next_cut_rc].reverse_complement()
        frag_rc.id = record.id

        if next_cut < cut:  # wraps around the zero-index
            frag = (record + record)[cut : next_cut + len(record)]
            frag.id = record.id
            frag_rc = (record + record)[
                cut_rc : next_cut_rc + len(record)
            ].reverse_complement()
            frag_rc.id = record.id

        frag_w_overhangs.append((left, frag, right))
        frag_w_overhangs.append((left_rc, frag_rc, right_rc))

    return frag_w_overhangs


def _has_features(record: SeqRecord, include: Optional[List[str]]) -> bool:
    """Return whether any of a record's features/qualifiers match the include specified.

    Args:
        record: the record being checked for include
        include: the include to filter for

    Returns:
        Whether the record has any features or qualifiers with specified include
    """

    if not include:
        return True

    assert isinstance(include, list)

    features: Set[str] = set()
    for feature in record.features:
        features.update(feature.id.lower().split())
        for _, values in feature.qualifiers.items():
            for value in values:
                features.update(value.lower().split())

    include_set = {i.lower() for i in include}
    intersect = features.intersection(include_set)

    return len(intersect) == len(include)
