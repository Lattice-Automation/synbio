"""Parse SnapGene files to dna and protein FASTA files."""

import os
import sys

from Bio.Seq import Seq
import dgparse  # https://github.com/DeskGen/dgparse

from .config import FEATURE_DIR


# TODO: add dgparse and fuzzywuzzy to the environment

if sys.version_info[0] > 2:
    raise RuntimeError("dgparse requires Python 2.")


SNAPGENEDIR = os.path.join(FEATURE_DIR, "snapgene")
DNA_FILE = os.path.join(SNAPGENEDIR, "dna.db")
PROTEIN_FILE = os.path.join(SNAPGENEDIR, "protein.db")


def parse():
    """Parse all Snapgene features to dna and protein FASTA files.

    Snapgene files were retrieved from:
    https://www.snapgene.com/resources/plasmid-files/
    """

    feat_files = [os.path.join(d, f) for d, _, fs in os.walk(SNAPGENEDIR) for f in fs]
    feat_files = [f for f in feat_files if f.endswith(".dna")]

    # a dict with key = name and value = seq (1 for each feature)
    failures = []
    protein_features = {}
    dna_features = {}
    feature_categories = {}
    feature_description = {}
    for file in feat_files:
        with open(file, "rb") as snapgene_file:
            try:
                parsed_file = dgparse.snapgene.parse(snapgene_file)
                for feature in parsed_file["dnafeatures"]:
                    dnafeature = feature["dnafeature"]

                    seq = dnafeature["pattern"]["bases"]
                    if not seq:
                        continue

                    name = dnafeature["name"].encode("utf-8").strip()
                    category = dnafeature["category"].encode("utf-8").strip()
                    description = ""

                    if "description" in dnafeature and dnafeature["description"]:
                        description = (
                            dnafeature["description"]
                            .encode("utf-8")
                            .strip()
                            .replace("<html><body>", "")
                            .replace("</body></html>", "")
                            .replace("<br>", "")
                        )

                    feature_categories[name] = category
                    feature_description[name] = description

                    is_protein = (
                        "CDS" == dnafeature["category"]
                        and "translation" in dnafeature["properties"]
                    )

                    if is_protein:
                        seq = Seq(seq)
                        seq = seq.translate()

                        while seq[-1] == "*":
                            seq = seq[:-1]

                        protein_features[name] = str(seq)
                    else:
                        dna_features[name] = seq
            except Exception as err:
                # issue in ~30 files: ascii' codec can't encode character u'\xe4' in position 25: ordinal not in range(128)
                failures.append(file)
                continue

    with open(DNA_FILE, "w") as dnafile:
        for name, seq in dna_features.iteritems():
            cat = feature_categories[name]
            desc = feature_description[name]
            dnafile.write(">" + name + "|" + cat + "|" + desc + "\n" + seq + "\n")

    with open(PROTEIN_FILE, "w") as proteinfile:
        for name, seq in protein_features.iteritems():
            cat = feature_categories[name]
            desc = feature_description[name]
            proteinfile.write(">" + name + "|" + cat + "|" + desc + "\n" + seq + "\n")


if __name__ == "__main__":
    parse()
