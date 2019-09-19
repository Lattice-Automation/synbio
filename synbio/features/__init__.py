"""A subpackage for annotating SeqRecords with common DNA and protein features.

These features are derived from a mix of SnapGene curated features and
iGEM. When teams upload parts to iGEM, they annotate regions of each sequence.
We clustered those features from manual annotation by iGEM teams and used
the clusters consensus sequence and names as additional features to augment
the SnapGene feature dataset with.

Clustering on the iGEM feature dataset was performed with cd-hit with parameters
defined in synbio.features.config. Consensus names were extracted from the
descriptions and names for each iGEM part. Clustered features that show up multiple times
are stored as features in the database accessed during annotation.

`synbio.features.annotate()` is a function that accepts a Bio.SeqRecord and returns
a new one with additional Bio.SeqFeatures from the feature database. It uses a
kmer seeding to seed alignments and then filters for all DNA and protein features
that exceed the identity ratio threshold (0.95 by default).
"""

from .annotate import annotate
