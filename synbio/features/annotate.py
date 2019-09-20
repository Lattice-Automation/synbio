"""Align query records against the DNA and protein database."""

from collections import defaultdict
import os
import pickle
from math import floor
from typing import Dict, List, Set, Tuple

from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from .config import (
    DNA_ID_MAP_PICKLE,
    DNA_KMER_MAP_PICKLE,
    PROTEIN_ID_MAP_PICKLE,
    PROTEIN_KMER_MAP_PICKLE,
)


if os.path.isfile(DNA_ID_MAP_PICKLE):
    with open(DNA_ID_MAP_PICKLE, "r+b") as KMER_MAP:
        DNA_ID_MAP = pickle.load(KMER_MAP)

    with open(DNA_KMER_MAP_PICKLE, "r+b") as KMER_MAP:
        DNA_KMER_MAP = pickle.load(KMER_MAP)

    with open(PROTEIN_ID_MAP_PICKLE, "r+b") as ID_MAP:
        PROTEIN_ID_MAP = pickle.load(ID_MAP)

    with open(PROTEIN_KMER_MAP_PICKLE, "r+b") as ID_MAP:
        PROTEIN_KMER_MAP = pickle.load(ID_MAP)


class Hit:
    """Hits are kmer matches between the query sequence and a feature.

    Attributes:
        kmer: The kmer sequence of the hit
        query_loc: The index of the hit on the query seq
        subject: Reference to the SeqRecord the hit came from
        subject_loc: The index of the hit on the subject SeqRecord
    """

    def __init__(self, kmer: str, query_loc: int, subject: SeqRecord, subject_loc: int):
        self.kmer = kmer
        self.query_loc = query_loc
        self.subject = subject
        self.subject_loc = subject_loc


class Match:
    """Matches are final matches between features in the DBs and the query sequence.

    Attributes:
        subject: The matched feature (as a SeqRecord)
        query_start: The start index on the query sequence
        query_end: The end index on the query sequence
        subject_start: The start index on the subject sequence
        subject_end: The end index on the subject sequence
    """

    def __init__(
        self,
        subject: SeqRecord,
        query_start: int,
        query_end: int,
        subject_start: int,
        subject_end: int,
    ):
        self.subject = subject
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end


def annotate(
    record: SeqRecord, identity: float = 0.95, circular: bool = True, cull: bool = True
) -> SeqRecord:
    """Create a new SeqRecord with additional DNA and protein features.

    This annotation function replaces BLAST, which Lattice used to use
    for the same purpose. Downsides of using BLAST are that it introduces a
    binary that has to be pre-compiled for each major OS and is, honestly,
    a bit of overkill for the use-case described here: adding common plasmid
    features to a SeqRecord. This annotation function is based SUPER loosely
    on BLAST's approach where features' kmers are stored in a database and
    are used to seed alignments against the query sequence.

    Args:
        record: The record to annotate with features

    Keyword Args:
        identity: The identity ratio threshold
            below which feature hits are ignored
        circular: Whether the record to annotate is circular or
            linear. Affects annotation across the 1-index
        cull: Whether to remove features that are completely or nearly fully
            engulfed in others. Features are removed if they have either
            $identity bp in common with another feature or are fully engulfed
            by another feature

    Returns:
        A new SeqRecord with additional DNA and protein features
    """

    if not 0.0 < identity <= 1.0:
        raise ValueError(
            "identity must be a ratio greater than 0 and less than or equal to 1"
        )

    if isinstance(record.seq, Seq) and isinstance(record.seq.alphabet, generic_protein):
        raise ValueError("annotation only supports SeqRecords with DNA Alphabets")

    new_record = record.upper()
    seq = str(new_record.seq.upper())

    # get DNA and protein features from the currated database
    dna_features = _get_features(
        seq, DNA_KMER_MAP, DNA_ID_MAP, identity, circular, False
    )
    protein_features = _get_features(
        seq, PROTEIN_KMER_MAP, PROTEIN_ID_MAP, identity, circular, True
    )

    # cull the new features to avoid highly overlapping ones
    new_features = dna_features + protein_features
    if cull:
        new_features = _cull(seq, new_features, identity)

    # add to the new record
    new_record.features.extend(new_features)
    return new_record


def _get_features(
    seq: str,
    kmer_map: Dict[str, List[Tuple[str, int]]],
    subject_map: Dict[str, SeqRecord],
    identity: float,
    circular: bool,
    protein: bool,
) -> List[SeqFeature]:
    """Get a list of SeqFeatures for the sequence.

    Gather aligned features in both the FWD and the REV direction.

    Args:
        seq: The query sequence
        kmer_map: Map from kmer to a list of tuples with:
            1. the ID of the source feature
            2. the start index of the feature in the source feature
        subject_map: Map from subject ID to the source SeqRecord/feature
        float: The identity ratio that should be the threshold
        circular: Whether the seq is circular
        protein: Whether to gather coding sequence annotations, in which
            case the query sequence needs to be translated into each
            of its three open reading frames

    Returns:
        A list of SeqFeatures that align with the query sequence
    """

    features: List[SeqFeature] = []

    def add_features(strand: bool):
        query = seq.upper()
        if not strand:
            query = str(Seq(seq).reverse_complement())

        frames = [query]
        if protein:
            frames = [query[i:] for i in [0, 1, 2]]
            frames = [f[: (len(f) // 3) * 3] for f in frames]
            frames = [str(Seq(f).translate()) for f in frames]

        for i, frame in enumerate(frames):
            hits = _get_hits(frame, kmer_map, subject_map, circular)
            hits_filtered = _filter_hits(frame, hits)
            hits_reduced = _reduce_hits(hits_filtered)
            matches = _get_matches(frame, hits_reduced, identity)

            for match in matches:
                feature_subject = match.subject
                feature_type = feature_subject.annotations["type"]

                query_start = match.query_start
                query_end = match.query_end + 1
                if not strand:
                    query_start = len(frame) - match.query_end - 1
                    query_end = len(frame) - match.query_start
                if protein:
                    query_start = (query_start * 3) + i
                    query_end = (query_end * 3) + i

                if query_end >= len(query):
                    # TODO: why aren't features across the 1-index
                    # supported by BioPython
                    continue

                strand_int = 1 if strand else -1

                feature_loc = FeatureLocation(query_start, query_end, strand=strand_int)

                feature = SeqFeature(
                    id=feature_subject.name or feature_subject.id,
                    location=feature_loc,
                    type=feature_type,
                    strand=strand_int,
                    qualifiers=feature_subject.annotations,
                )
                features.append(feature)

    add_features(True)
    add_features(False)

    return features


def _get_hits(
    seq: str,
    kmer_map: Dict[str, List[Tuple[str, int]]],
    subject_map: Dict[str, SeqRecord],
    circular: bool,
) -> List[Hit]:
    """Gather all the matches from a seq's kmers in a kmer_map.

    Args:
        seq: The query sequence
        kmer_map: Map from kmer to a list of tuples with:
            1. the ID of the source feature
            2. the start index of the feature in the source feature
        subject_map: Map from subject ID to the source SeqRecord/feature
        circular: Whether the query sequence is circular

    Returns:
        A list of Hits
    """

    assert kmer_map and subject_map

    first_key = list(kmer_map.keys())[0]
    word_size = len(first_key)

    query = seq + seq if circular else seq

    hits: List[Hit] = []
    for i in range(len(query) - word_size + 1):
        if i >= len(seq):  # past wrapping around 1-index
            break

        kmer = query[i : i + word_size]
        if kmer not in kmer_map:
            continue

        for subject_id, subject_loc in kmer_map[kmer]:
            feature = subject_map[subject_id]

            # HACK: this will fail if, for example, all but 1bp
            # of a feature align with the query seq (and the
            # 1bp extends out to the left/right of the query seq)

            # it starts to the left of 1-index
            if i - subject_loc < 0:
                continue

            # it ends to the right of the end of query seq
            feature_left = len(feature) - subject_loc
            if feature_left + i >= len(query):
                continue

            hit = Hit(kmer, i, subject_map[subject_id], subject_loc)
            hits.append(hit)

    return hits


def _filter_hits(seq: str, hits: List[Hit]) -> List[Hit]:
    """Remove hits that don't show up enough to reach identity threshold.

    This is just a small heuristic method to avoid the hit expansion/DP
    part of alignment. The alignment/expansion part is expensive
    and this cuts down the number of hits that have to be expanded
    by looking ahead to check if there's any chance of there being
    enough bp for it to be worth it. For example, if there's ONLY one
    hit for GFP in a sequence but GFP has a 500bp sequence, it's not
    worth expanding it and we can skip that

    Args:
        seq: The length of the query sequence
        hits: The hits to filter

    Returns:
        Hits that meet the threshold we're filtering on
    """

    if not hits:
        return []

    word_size = len(hits[0].kmer)

    hit_map: Dict[str, SeqRecord] = {}
    hit_count: Dict[str, int] = defaultdict(int)

    for hit in hits:
        hit_map[hit.subject.id] = hit.subject
        hit_count[hit.subject.id] += 1

    hit_subjects: Set[str] = set()
    for subject_id, count in hit_count.items():
        subject_len = len(hit_map[subject_id])

        subject_threshold = floor(subject_len / float(word_size + 1))
        if count >= subject_threshold and subject_len < len(seq):
            hit_subjects.add(subject_id)

    hits_filtered: List[Hit] = []
    for hit in hits:
        if hit.subject.id in hit_subjects:
            hits_filtered.append(hit)

    return hits_filtered


def _reduce_hits(hits: List[Hit]) -> List[Hit]:
    """Get rid of redundant hits representing the same stretch.

    If hits all correspond to the same range on the query sequence
    the hits can be reduced to just the first hit since it will be
    expanded and engulf the other hits during expansion.

    Args:
        hits: Hits that we're reducing into a non-redundant set

    Returns:
        A new list of Hits reduced down to one hit per possible Match
    """

    hits_unique: Dict[str, Hit] = {}
    for hit in hits:
        # get the feature's probable start index on the subject
        query_start = hit.query_loc - hit.subject_loc
        hit_id = hit.subject.id + str(query_start)
        if hit_id not in hits_unique:
            hits_unique[hit_id] = hit

    return list(hits_unique.values())


def _get_matches(seq: str, hits: List[Hit], identity: float) -> List[Match]:
    """Expand Hits with an alignment against the query sequence, turn into Matches.

    Filter out hits/matches that fall beneath the identity threshold after expansion.

    Args:
        seq: The query sequence
        hits: Filtered hits to be expanded
        identity: The identity threshold below which features are ignored

    Returns:
        Matches with ranges on the query sequence
    """

    query = seq + seq  # double for if circular across the 1-index

    hits = sorted(hits, key=lambda hit: hit.query_loc)

    matches: List[Match] = []
    for hit in hits:
        i = hit.query_loc - hit.subject_loc  # index on query seq
        j = 0  # index on subject feature

        score = 0  # alignment bp matches
        threshold_hit = identity * len(hit.subject)
        threshold_miss = (1 - identity) * len(hit.subject)

        query_start = -1  # first bp in query seq
        query_end = 0  # last bp match on query seq
        subject_start = -1  # first bp on subject feature
        subject_end = 0  # last bp match on subject feature
        while j < len(hit.subject):
            # add to score if it's a matched bp
            if query[i] == hit.subject[j]:
                if query_start < 0:
                    query_start = i
                    subject_start = j
                score += 1
                query_end = i
                subject_end = j
            i += 1
            j += 1

            # bail if we're too low at this point in aligning
            if j % 10 == 0 and j - score >= threshold_miss:
                break

        if score >= threshold_hit:
            match = Match(
                hit.subject, query_start, query_end, subject_start, subject_end
            )
            matches.append(match)

    return matches


def _cull(seq: str, features: List[SeqFeature], identity: float) -> List[SeqFeature]:
    """Remove features that are nearly or completely engulfed by others.

    Features are removed if any of the following occurs:
        1. 100% engulfed by another feature
        2. percentage identity bp in common with another feature

    Args:
        seq: The query sequence
        features: Features to cull
        identity: The ratio identity used during aligning

    Returns:
        Culled features
    """

    features = sorted(features, key=len, reverse=True)

    hits: List[List[str]] = []
    for _ in range(len(seq)):
        hits.append([])

    features_unculled: List[SeqFeature] = []
    for feature in features:
        overlaps: Dict[str, int] = defaultdict(int)
        for i in range(int(feature.location.start), int(feature.location.end) + 1):
            for overlap in hits[i]:
                overlaps[overlap] += 1
            hits[i].append(feature.id)

        cull = False
        overlap_threshold = identity * len(feature)
        for overlap, count in overlaps.items():
            if count >= overlap_threshold:
                cull = True
                break

        if not cull:
            features_unculled.append(feature)

    return features_unculled
