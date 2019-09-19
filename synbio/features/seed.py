"""Create a kmer map for each DNA and protein sequence."""

from collections import defaultdict
from hashlib import sha1
import pickle
from typing import Dict, List, Set, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .config import (
    DNA_DB,
    DNA_WORD_SIZE,
    DNA_ID_MAP_PICKLE,
    DNA_KMER_MAP_PICKLE,
    PROTEIN_DB,
    PROTEIN_WORD_SIZE,
    PROTEIN_ID_MAP_PICKLE,
    PROTEIN_KMER_MAP_PICKLE,
)


def seed():
    """Create kmer maps from ID to feature (SeqRecord) and kmers to feature IDs.

    This does two things:
        1. Parse the feature database (FASTA) and store each
            in a map from random unique IDs for each feature
            to the feature (SeqRecord). dna_id_map, protein_id_map.
            Saves this map to "dna.id.pickle", and protein
        2. Chop each feature into kmers and store in a map
            where each key is the kmer (str) and the value
            is a list of feature ids for features that have
            that kmer (which is what the feature adds its ID to).
            Saves this kmer map to "dna.kmermap.pickle", and protein
        
    This is SUPER loosely based on BLAST's initial word search approach.
    """

    dna_id_map = _id_map(DNA_DB, DNA_WORD_SIZE)
    protein_id_map = _id_map(PROTEIN_DB, PROTEIN_WORD_SIZE)

    with open(DNA_ID_MAP_PICKLE, "wb") as id_map_file:
        pickle.dump(dna_id_map, id_map_file)
    with open(PROTEIN_ID_MAP_PICKLE, "wb") as id_map_file:
        pickle.dump(protein_id_map, id_map_file)

    dna_kmer_map = _kmer_map(dna_id_map, DNA_WORD_SIZE)
    protein_kmer_map = _kmer_map(protein_id_map, PROTEIN_WORD_SIZE)

    with open(DNA_KMER_MAP_PICKLE, "wb") as kmer_map_file:
        pickle.dump(dna_kmer_map, kmer_map_file)
    with open(PROTEIN_KMER_MAP_PICKLE, "wb") as kmer_map_file:
        pickle.dump(protein_kmer_map, kmer_map_file)


def _id_map(filename: str, word_size: int) -> Dict[str, SeqRecord]:
    """Read in the database and create a map from unique ID to SeqRecord
    
    Args:
        filename: DNA or protein database filename
        word_size: Minumum length word size for a feature
    
    Returns:
        A map from unique ID (random) to the SeqRecord with a sequence
    """

    id_map: Dict[str, SeqRecord] = {}

    names_seen: Set[str] = set()

    def key(fname: str) -> str:
        fname = fname.lower()
        fname = fname.replace("origin", "ori")
        return fname

    with open(filename, "r") as db:
        db_lines = db.readlines()
        for i, line in enumerate(db_lines):
            if i % 2 == 1:
                continue

            line = line[1:].strip()  # starts with >
            name, ftype, desc = line.split("|")
            seq = db_lines[i + 1].strip()

            if len(seq) < word_size:
                continue

            name_key = key(name)
            if name_key in names_seen:
                continue
            names_seen.add(name_key)

            record = SeqRecord(
                Seq(seq),
                id=str(sha1(seq.encode()).hexdigest())[:8],
                name=name,
                description=desc,
                annotations={"type": ftype},
            )

            id_map[record.id] = record

    return id_map


def _kmer_map(
    id_map: Dict[str, SeqRecord], word_size: int
) -> Dict[str, List[Tuple[str, int]]]:
    """Create a map from kmer to the id of the matched feature
    and the index of the match within the matched feature.

    Args:
        id_map: Map from SeqRecord ID to SeqRecord from _id_map
        word_size: The length of each kmer word
    
    Returns:
        A map from kmer to a list of tuples with:
            1. the ID of the feature
            2. the 0-based index of feature's start index
    """

    kmer_map: Dict[str, List[Tuple[str, int]]] = defaultdict(list)

    for fid, feature in id_map.items():
        seq = str(feature.seq)
        for i in range(len(seq) - word_size + 1):
            sub_seq = seq[i : i + word_size]

            kmer_map[sub_seq].append((fid, i))

    return kmer_map
