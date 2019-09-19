"""Configuration for feature clustering and seeding."""

from os import path

# Dirnames and filenames
DIR_NAME = path.dirname(path.realpath(__file__))
FEATURE_DIR = path.join(DIR_NAME, "..", "..", "data", "features")
IGEM_DIR = path.join(FEATURE_DIR, "igem")
DNA_DB = path.join(FEATURE_DIR, "dna.db")
DNA_ID_MAP_PICKLE = path.join(FEATURE_DIR, "dna.id.pickle")
DNA_KMER_MAP_PICKLE = path.join(FEATURE_DIR, "dna.kmermap.pickle")
PROTEIN_DB = path.join(FEATURE_DIR, "protein.db")
PROTEIN_ID_MAP_PICKLE = path.join(FEATURE_DIR, "protein.id.pickle")
PROTEIN_KMER_MAP_PICKLE = path.join(FEATURE_DIR, "protein.kmermap.pickle")

# Parameters
DNA_WORD_SIZE = 11
DNA_IDENTITY_THRESHOLD = 0.95  # -c from cd-hit
DNA_LENGTH_DISTANCE_CUTOFF = 0.9  # -s from cd-hit
PROTEIN_WORD_SIZE = 5
PROTEIN_IDENTITY_THRESHOLD = 0.9
PROTEIN_LENGTH_DISTANCE_CUTOFF = 0.85  # -s from cd-hit
CLUSTER_SIZE_THRESHOLD = 3  # min number of cluster members to become feature
CLUSTER_SOURCES_CUTOFF = 2  # min number of part source to become a feature
RED_FLAG_NAMES = ["ori", "cds"]
RED_FLAG_NAMES_INNER = [
    "primer",
    "bba_",
    "spacer",
    "...",
    "attl",
    "attr",
    "attb",
    "attp",
    "coding sequence",
    "linker",
    "{&amp;#706;KILR&amp;#707;}",
    "tandem repeat",
    "inverted repeat",
    "a repeat",
    "no description",
    "][",
    "lacP\\",
    ";",
    "#",
]  # name substrings that are red-flags
