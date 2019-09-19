"""Cluster iGEM annotated features with cd-hit to create DNA and protein feature databases."""

from collections import defaultdict
import functools
import os
import pickle
import re
import subprocess
from typing import Dict, Optional, List, Tuple, Set

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from fuzzywuzzy import fuzz

from .config import (
    DNA_WORD_SIZE,
    DNA_IDENTITY_THRESHOLD,
    DNA_LENGTH_DISTANCE_CUTOFF,
    IGEM_DIR,
    PROTEIN_WORD_SIZE,
    PROTEIN_IDENTITY_THRESHOLD,
    PROTEIN_LENGTH_DISTANCE_CUTOFF,
    CLUSTER_SIZE_THRESHOLD,
    CLUSTER_SOURCES_CUTOFF,
    RED_FLAG_NAMES,
    RED_FLAG_NAMES_INNER,
)


# File names
IGEM = os.path.join(IGEM_DIR, "igem.2017.xml")
ID_TO_RECORD = os.path.join(IGEM_DIR, "igem.pickle")
DNA = os.path.join(IGEM_DIR, "dna.fa")
DNA_CLSTR = os.path.join(IGEM_DIR, "dna")
DNA_DB = os.path.join(IGEM_DIR, "dna.db")
PROTEIN = os.path.join(IGEM_DIR, "protein.fa")
PROTEIN_CLSTR = os.path.join(IGEM_DIR, "protein")
PROTEIN_DB = os.path.join(IGEM_DIR, "protein.db")

# Regexes for parsing rows from the iGEM XML dump (SHA1s hinder XML parsing)
RE_NAME = re.compile(r"<field name=\"part_name\">(\w*)</field>")
RE_SEQ = re.compile(r"<field name=\"sequence\">(\w*)</field>")
RE_SHORT_DESC = re.compile(r"<field name=\"short_desc\">(.*)</field>")
RE_DESC = re.compile(r"<field name=\"description\">(.*)</field>")
RE_CATS = re.compile(r"<field name=\"categories\">(.*)</field>")
RE_CACHE = re.compile(r"<field name=\"seq_edit_cache\">(.*)")
RE_TYPE = re.compile(r"<field name=\"part_type\">(.*)</field>")
RE_NICKNAME = re.compile(r"<field name=\"nickname\">(.*)</field>")
RES = [RE_NAME, RE_SEQ, RE_SHORT_DESC, RE_DESC, RE_CATS, RE_CACHE, RE_TYPE, RE_NICKNAME]
RE_FEATURES = re.compile(r"seqFeatures = new Array\((.+?)\)")


Cluster = List[SeqRecord]
Clusters = List[Cluster]


def cluster():
    """Create feature databases.

    - Make a map from iGEM part ID to SeqRecord
        - Gather features out of each SeqRecord's "seq_edit_cache" script
        - Remove all features that are redundant w/ an iGEM part
    - Make a DNA/protein databases from each SeqRecord and feature in /data
    - Cluster the DNA/protein databases with cdhit
    - Gather all DNA/protein clusters above some threshold
    - Find common names for each feature in a cluster using their description
    """

    # parse the iGEM XML file
    if all(os.path.exists(f) for f in [DNA, PROTEIN, ID_TO_RECORD]):
        with open(ID_TO_RECORD, "r+b") as record_file:
            id_to_record = pickle.load(record_file)
    else:
        id_to_record = _parse_igem_data()
        _write_databases(id_to_record)
        with open(ID_TO_RECORD, "w+b") as record_file:
            pickle.dump(id_to_record, record_file)

    # run cd-hit
    if not os.path.exists(DNA_CLSTR + ".clstr"):
        _cluster_databases()

    # read in the SeqRecord/SeqFeature clusters
    dna_clusters, protein_clusters = _parse_clusters(id_to_record)

    # create a consensus SeqRecord for each cluster
    dna_features, protein_features = _create_features(dna_clusters, protein_clusters)

    # store
    with open(DNA_DB, "w") as dna_db:
        for feature in dna_features:
            dna_db.write(_fasta(feature))
    with open(PROTEIN_DB, "w") as protein_db:
        for feature in protein_features:
            protein_db.write(_fasta(feature))


def _parse_igem_data() -> Dict[str, SeqRecord]:
    """Parse the iGEM XML into a Dict with part_id -> SeqRecord.

    Args:
        test: Is a test, truncate the number of rows

    Returns:
        A map from each part's ID to a SeqRecord representing it.
    """

    id_to_record: Dict[str, SeqRecord] = {}
    with open(IGEM, "r") as data:
        for row in data.read().split("<row>"):
            record = _parse_row(row)
            if record:
                id_to_record[record.id] = record

    assert "BBa_J23104" in id_to_record

    return id_to_record


def _parse_row(row: str) -> Optional[SeqRecord]:
    """Parse a single row of the iGEM XML into a SeqRecord.

    I'm using the regex package here because there are characters
    in the XML's sha1 that throw the XML parsers that I tried. I know
    that this regex is slow...

    Args:
        row: A single 'row' element in the XML

    Returns:
        A SeqRecord with the part id, sequence, and description stored
    """

    def match(regex):
        row_match = regex.search(row)
        if row_match:
            return row_match[1]
        return ""

    matches = [match(r) for r in RES]
    name, seq, desc_short, desc_long, cats, cache, ftype, nickname = matches

    if not seq or len(seq) > 10_000:
        return None

    features: List[SeqFeature] = []
    feature_matches = RE_FEATURES.search(cache)
    if feature_matches:
        feature_match = feature_matches[1]
        for feature in feature_match.split("]"):
            feature = feature.replace(", [", "")

            if feature.count(",") < 4:
                continue

            f_type, f_start, f_end, f_name, f_strand, *_ = [
                f.replace("[", "").replace("'", "").replace("(", "").strip()
                for f in feature.split(",")
            ]

            f_start_int = int(f_start)
            f_end_int = int(f_end)

            if (
                f_start_int == f_end_int
                or f_start_int > f_end_int
                or f_end_int - f_start_int < DNA_WORD_SIZE
            ):
                continue

            features.append(
                SeqFeature(
                    id=f_name,
                    # have to -1 here. more CDS are of a length % 3 == 0 w/ this
                    # I don't think this was enforced on iGEM teams when making features
                    location=FeatureLocation(
                        f_start_int - 1, f_end_int, 1 if f_strand == "0" else -1
                    ),
                    type=_get_type(f_name, f_type.lower()),
                    strand=1 if f_strand == "0" else -1,
                )
            )

    return SeqRecord(
        Seq(seq, IUPACUnambiguousDNA()),
        id=name,
        dbxrefs=[name],
        annotations={
            "short_desc": desc_short,
            "description": desc_long,
            "categories": cats,
            "nickname": nickname,
            "type": _get_type(name, ftype),
        },
        features=features,
    )


def _get_type(name: str, ftype: str):
    """Get the feature type from its name and type
    
    Args:
        name: Feature's name
        ftype: Feature's type
    """

    type_map: Dict[str, str] = {
        "coding": "CDS",
        "protein": "CDS",
        "cds": "CDS",
        "reporter": "CDS",
        "dna": "misc_feature",
        "misc_p": "misc_feature",
        "rbs": "RBS",
        "h_pin": "misc_RNA",
        "stop": "terminator",
        "polyA_signal": "polya",
    }

    # a map from bad types to better ones
    if ftype in type_map:
        ftype = type_map[ftype]

    for t in ["promoter", "terminator"]:
        if t in name.lower():
            ftype = t
            break

    return ftype


def _write_databases(id_to_record: Dict[str, SeqRecord]):
    """Create the DNA and protein databases for each SeqRecord.

    Write both databases to a tmp directory, return paths

    Args:
        id_to_record: Dictionary from iGEM part ID to SeqRecord

    Returns:
        A tuple with two file names. The first is for the DNA database
        and the second is for the translated protein database
    """

    records = {r for r in id_to_record.keys()}

    with open(DNA, "w") as dna_file, open(PROTEIN, "w") as protein_file:
        for record in id_to_record.values():
            dna_file.write(f">{record.id}\n{str(record.seq)}\n")
            if _is_cds(record):
                protein_file.write(f">{record.id}\n{_translate(record.seq)}\n")

            for i, feature in enumerate(record.features):
                if (
                    feature.id in records or feature.type == "brick"
                ):  # would be redundant
                    continue

                feature_seq = feature.extract(record.seq)
                dna_file.write(f">{record.id}.{str(i)}\n{str(feature_seq)}\n")
                if (
                    feature.type in ("CDS", "protein", "coding")
                    and _translate(feature_seq)
                    and _translate(feature_seq).count("*") < 1
                ):
                    protein_file.write(
                        f">{record.id}.{str(i)}\n{_translate(feature_seq)}\n"
                    )


def _is_cds(record: SeqRecord) -> bool:
    """Returns whether a CDS should be translated for the protein database

    Args:
        record: The SeqRecord to check for being a CDS

    Returns:
        Whether it's a probable CDS
    """

    if record.annotations["type"] == "CDS":
        return True

    if (
        "CDS" in record.annotations["categories"]
        and "ATG" == record.seq[:3].upper()
        and record.seq.translate()
        and str(record.seq.translate()).count("*") < 2
    ):
        return True

    return False


def _translate(seq: Seq) -> str:
    """Translate a sequence from DNA to protein in way that cd-hit expects

    Args:
        seq: DNA sequence from a SeqRecord

    Returns:
        Protein sequence for a protein database
    """

    aas = seq.translate()
    while aas and aas[-1] == "*":  # cd-hit dislikes stop codons
        aas = aas[:-1]
    return str(aas)


def _cluster_databases():
    """Run mmseqs on the DNA and protein databases

    Args:
        dna_file: The name of the DNA file
        protein_file: The name of the protein file

    Returns:
        The name of the DNA and protein cluster file
    """

    # see: https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#Hierarchically_clustering
    subprocess.call(
        [
            "cd-hit-est",
            "-i",
            DNA,
            "-o",
            DNA_CLSTR,
            "-c",  # 90% sequence identity threshold
            str(DNA_IDENTITY_THRESHOLD),
            "-n",  # word size
            str(DNA_WORD_SIZE),
            "-s",
            str(DNA_LENGTH_DISTANCE_CUTOFF),
            "-d",  # use sequence name in fasta header for description
            "0",
            "-T",  # use all threads available
            "0",
            "-M",  # use 12GB of memory
            "12000",
            "-g",  # slower but more accurate
            "1",
        ]
    )

    subprocess.call(
        [
            "cd-hit",
            "-i",
            PROTEIN,
            "-o",
            PROTEIN_CLSTR,
            "-c",  # 85% sequence identity threshold
            str(PROTEIN_IDENTITY_THRESHOLD),
            "-n",  # word size of 5
            str(PROTEIN_WORD_SIZE),
            "-s",
            str(PROTEIN_LENGTH_DISTANCE_CUTOFF),
            "-d",  # use sequence name in fasta header for description
            "0",
            "-T",  # use all threads available
            "0",
            "-M",  # use 12GB of memory
            "12000",
        ]
    )


def _parse_clusters(id_to_record: Dict[str, SeqRecord]) -> Tuple[Clusters, Clusters]:
    """Parse clusters to gather reference seq and other meta

    - Filter clusters to those with at least $CLUSTER_SIZE_THRESHOLD members
    - Find the reference records/features for each cluster

    Example of a cluster from cd-hit:

    ```txt
        >Cluster 1865
        0	3738nt, >BBa_K1344008.0... at -/99.55%
        1	3707nt, >BBa_K2449015... at +/100.00%
        2	3792nt, >BBa_K2449018... *
        3	3707nt, >BBa_K2449018.0... at -/100.00%
        4	3707nt, >BBa_K2449019.0... at -/100.00%
        5	3792nt, >BBa_K2449020.0... at -/100.00%
        6	3707nt, >BBa_K2449020.1... at -/100.00%
    ```

    The first SeqRecord/SeqFeature in the returned cluster list is
    the reference sequence that's refered to in elsewhere

    Args:
        id_to_record: A dictionary from SeqRecord's id to the SeqRecord

    Returns:
        Two lists of clusters. The first list is for DNA, the second is for protein
    """

    def parse_cluster(src, protein=False):
        clusters: Clusters = []
        with open(src + ".clstr", "r") as dnafile:
            for clstr_str in dnafile.read().split(">Cluster"):
                if clstr_str.count("\n") - 1 < CLUSTER_SIZE_THRESHOLD:
                    continue

                cluster_group: Cluster = []
                for line in clstr_str.split("\n")[1:]:
                    if ">" not in line or "..." not in line:
                        continue

                    line_id = line[line.index(">") + 1 : line.index("...")]

                    rid = line_id  # record id
                    fid = ""  # feature id
                    findex = -1
                    if "." in rid:
                        rid, fid = rid.split(".")
                        findex = int(fid)  # index of feature on record

                    record = id_to_record[rid]
                    member = record.upper()  # copy member here

                    if findex > -1:
                        # turn the feature into a SeqRecord
                        feature = record.features[findex]
                        member = SeqRecord(
                            feature.extract(record.seq),
                            id=feature.id,
                            dbxrefs=[record.id],  # temporary
                            annotations={
                                "short_desc": "",
                                "description": "",
                                "categories": "",
                                "nickname": "",
                                "type": feature.type,
                            },
                        )

                    if protein:
                        member.seq = Seq(_translate(member.seq))
                    cluster_group.append(member)

                assert len(cluster_group) >= CLUSTER_SIZE_THRESHOLD
                clusters.append(cluster_group)
        assert any(len(c) == CLUSTER_SIZE_THRESHOLD for c in clusters)
        return clusters

    return (
        parse_cluster(DNA_CLSTR, protein=False),
        parse_cluster(PROTEIN_CLSTR, protein=True),
    )


def _create_features(
    dna_clusters: Clusters, protein_clusters: Clusters
) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Create the lists of consensus DNA and protein features

    1. Vote on feature name
    2. Gather categories associated with the SeqRecord
    3. Set 'dbxrefs' as a list of references in iGEM

    Args:
        dna_clusters: The list of DNA clusters
        protein_clusters: The list of protein clusters

    Returns:
        Two lists. One of DNA features and one of protein features
    """

    dna_features = [_consensus_record(c) for c in dna_clusters]
    protein_features = [_consensus_record(c) for c in protein_clusters]

    dna_features = [f for f in dna_features if f]
    protein_features = [f for f in protein_features if f]

    return dna_features, protein_features


def _consensus_record(records: Cluster) -> Optional[SeqRecord]:
    """Deduplicate a list of clustered SeqRecords into a single "consensus" SeqRecord

    Args:
        records: A list of SeqRecords that were clustered via cd-hit

    Returns:
        A single SeqRecord with a name, type, categories that reflect the cluster
    """

    name = _consensus_name(records)
    seq = records[0].seq.upper()
    record_type = _consensus_attribute(records, "type", 1) or "misc_feature"
    record_type = _get_type(name, record_type.lower())
    description = _consensus_attribute(records, "description", 1)
    if len(description) < 10:
        description = ""

    description = description.replace("|", " ")

    # avoid really short features
    if len(seq) < DNA_WORD_SIZE:
        return None

    # avoid BS names
    if len(name) > 50:
        return None

    if "00" in name:
        return None  # almost 100% iGEM

    for redflag in RED_FLAG_NAMES_INNER:
        if redflag.lower() in name.lower():
            return None

    if name.lower() in RED_FLAG_NAMES:
        return None

    cat_set = {
        cat
        for r in records
        for cat in r.annotations["categories"].split()
        if "categories" in r.annotations
    }

    dbxrefs: Set[str] = set()
    for record in records:
        for db in record.dbxrefs:
            dbxrefs.add(db)

    src_count = len({src[:7] for src in dbxrefs})
    if src_count < CLUSTER_SOURCES_CUTOFF:
        return None

    return SeqRecord(
        seq,
        id=records[0].id,
        name=name,
        annotations={
            "categories": ",".join(cat_set),
            "type": record_type,
            "description": description,
        },
        dbxrefs=sorted(list(dbxrefs)),
    )


def _consensus_name(records: Cluster) -> str:
    """Get the consensus name from a list of possible name sources.

    1. If there is a consensus nickname, use that
    2. Otherwise, gather all short_desc, description and nicknames,
        take the 6 shortest and use the one with the smallest
        token set ratio distance from the others

    Args:
        records: The Records whose feature we want a consensus name for

    Returns:
        The guessed consensus name of the feature
    """

    # nicknames take precedent
    consensus_nickname = _consensus_attribute(records, "nickname")
    if consensus_nickname and consensus_nickname not in RED_FLAG_NAMES_INNER:
        return consensus_nickname

    entries: List[str] = []
    for record in records:
        entries.append(record.id)
        for src in ["short_desc", "description", "nickname"]:
            if src in record.annotations and record.annotations[src].strip():
                entries.append(record.annotations[src])

    # filter entries for those with plausible names
    def keep(entry: str) -> bool:
        if len([e for e in entry if e.isalnum()]) < 3:
            return False
        if entry[:3] == "BBa" and entry.count(" ") == 0:
            return False
        return True

    entries_filtered = [e for e in entries if keep(e)]
    if not entries_filtered:
        entries = [e for e in entries if e]
        if not entries:
            return records[0].id.strip()  # give up, all empty, keep consensus id

    # heuristic so this doesn't take ages
    entry_max_count = 10
    if len(entries) > entry_max_count:
        # try to filter to just those under at or under a small number of words
        entries_short = [e for e in entries if len(e.split()) <= 5]
        if entries_short:
            entries = entries_short
        else:
            entries = sorted(entries, key=len)
            entries = entries[:entry_max_count]
    entries = sorted(entries)  # for functools.lru

    set_ratio_matrix: List[List[float]] = []
    for _ in range(len(entries)):
        set_ratio_matrix.append([100.0] * len(entries))

    # calc token set ratio of each pairing
    for i, entry in enumerate(entries):
        for j in range(i + 1, len(entries)):
            other = entries[j]
            set_ratio = _token_set_ratio(entry, other)
            set_ratio_matrix[i][j] = set_ratio
            set_ratio_matrix[j][i] = set_ratio

    # find the row with the max summed set ratio
    max_row = 0
    max_sum = sum(set_ratio_matrix[0])
    for i, row in enumerate(set_ratio_matrix):
        row_sum = sum(row)
        if row_sum > max_sum or (
            row_sum == max_sum and len(entries[i]) < len(entries[max_row])
        ):
            max_row = i
            max_sum = row_sum

    return entries[max_row].strip()


def _consensus_attribute(records: Cluster, field: str, threshold=2) -> str:
    """Gather a consensus nickname if there is one

    Args:
        records: the cluster whose nickname we want

    Returns:
        The consensus nickname, empty string if there is none
    """

    values: List[str] = []
    for record in [
        r for r in records if field in r.annotations and r.annotations[field]
    ]:
        values.append(record.annotations[field])

    value_to_count: Dict[str, int] = defaultdict(int)
    for value in values:
        value_to_count[value] += 1

    sorted_values = sorted(value_to_count.items(), key=lambda x: x[1], reverse=True)

    if sorted_values and sorted_values[0][1] >= threshold:
        return re.sub("<[^<]+?>", "", sorted_values[0][0], count=-1)
    return ""


@functools.lru_cache(maxsize=32768)
def _token_set_ratio(str1: str, str2: str) -> float:
    """Get the token set ratio. Memoize results."""

    return fuzz.token_set_ratio(str1, str2)


def _fasta(record: SeqRecord) -> str:
    """Create a FASTA from a SeqRecord for storing features to DB

    Args:
        record: A SeqRecord to store to a database file

    Returns:
        A string representation of the SeqRecord for storing
    """

    return f">{record.name}|{record.annotations['type']}|{record.annotations['description']}\n{str(record.seq)}\n"
