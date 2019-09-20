"""Synbio designs.

Design encode plasmids, libraries (of plasmids) or combinatorial builds
for DNA assembly. They are the design specification used for the assembly step.
"""

from typing import Dict, Iterable, List, Union

from Bio.SeqRecord import SeqRecord

from .containers import content_id


class Design:
    """A SynBio design specification.

    Attributes:
        linear: Whether the SeqRecords of the design are linear
    """

    __name__ = "design"
    linear = True

    def __iter__(self):
        raise NotImplementedError

    def append(self, records: List[SeqRecord]):
        """Add to an existing design."""

        raise NotImplementedError

    def get_all_records(self) -> List[SeqRecord]:
        """Return all unique SeqRecords across all designs.

        Returns:
            The unique set of SeqRecord in this Design.
        """

        id_to_record: Dict[str, SeqRecord] = {}
        for record_set in self:
            for record in record_set:
                record_id = content_id(record)
                id_to_record[record_id] = record
        return list(id_to_record.values())

    def __len__(self) -> int:
        """Return the number of SeqRecord combinations."""

        return 0


class Plasmid(Design):
    """Create a plasmid from a list of SeqRecords

    Attributes:
        records: SeqRecords to concatenate into a plasmid
        linear: Whether the SeqRecords of the design are linear
    """

    __name__ = "plasmid"

    def __init__(self, records: Iterable[SeqRecord] = None, linear: bool = True):
        self.records = list(records) if records else []
        self.linear = linear

    def __iter__(self) -> Iterable[List[SeqRecord]]:
        return iter([self.records])

    def __str__(self) -> str:
        """Return the ids of each piece of DNA in this plasmid

        Returns:
            A high level description of this plasmid's design
        """

        if not self.records:
            return "[]"

        records_str: List[str] = []
        for record in self.records:
            if record.id == "<unknown id>":
                records_str.append(f"{len(record.seq)}bp")
            else:
                records_str.append(record.id)

        return "[" + ",".join(records_str) + "]"

    def __len__(self) -> int:
        """Return the number of inserts in the plasmid."""

        return len(self.records)

    def append(self, records: Union[SeqRecord, List[SeqRecord]]):
        """Add the SeqRecord to this plasmid's list of SeqRecords/fragments

        Args:
            record: the SeqRecord to add to the plasmid specification
        """

        if not records:
            raise ValueError(
                f"Can only add SeqRecords or non-empty lists of SeqRecords to Design"
            )

        if not isinstance(records, list):
            records = [records]
        for record in records:
            self.records.append(record)


class Combinatorial(Plasmid):
    """A list of SeqRecords. Find and use all valid combinations.  

    Structure-wise it's the same as a Plasmid but multiple assemblies
    are possible given a single set of SeqRecords. List of SeqRecords in,
    list of SeqRecords out. Versus List of SeqRecords in one SeqRecord
    out for a Plasmid design

    Attributes:
        records: list of SeqRecord
        linear: Whether the SeqRecords of the design are linear
    """

    __name__ = "combinatorial"


class CombinatorialBins(Design):
    """A list of bins of SeqRecords. All combinations between bins are attempted.

    Attributes:
        bins: list of SeqRecord bins
        linear: Whether the SeqRecords of the design are linear
    """

    __name__ = "combinatorial_bins"

    def __init__(self, bins: Iterable[List[SeqRecord]] = None, linear: bool = True):
        if bins:
            first_bin = list(bins)[0]

            if not isinstance(first_bin, list):
                raise ValueError(
                    "CombinatorialBins requires a list of SeqRecord bins (lists).\n"
                    + f"First bin is: {type(first_bin)}"
                )

        self.bins: List[List[SeqRecord]] = list(bins) if bins else []
        self.linear = linear

    def __iter__(self):
        """Create all combinations of bin SeqRecords. Return each combo as a list."""

        if not self.bins:
            return iter([])

        # separate each of the first bin into the start of new assemblies
        record_sets: List[List[SeqRecord]] = [[r] for r in self.bins[0]]
        for i in range(1, len(self.bins)):
            new_record_sets: List[List[SeqRecord]] = []
            for record_set in record_sets:
                for record in self.bins[i]:
                    new_record_sets.append(record_set + [record])

            record_sets = new_record_sets

        return iter(record_sets)

    def __str__(self) -> str:
        """Return high level description of the number of records at each bin

        Example:
        ```[[4 x SeqRecord], [500 x SeqRecord], [3 x SeqRecord]]```

        Returns:
            high level overview of the number of SeqRecords in each bin
        """

        return (
            f"[{', '.join(['[' + str(len(b)) + ' x SeqRecord]' for b in self.bins])}]"
        )

    def __len__(self) -> int:
        """Return the length of the bins."""

        return len(self.bins)

    def append(self, records: Union[SeqRecord, List[SeqRecord]]):
        """Add the SeqRecords as a bin for combinatorial design

        Args:
            record: the records to add as a bin
        """

        if not records:
            raise ValueError(
                f"Can only add SeqRecords or non-empty lists of SeqRecords to Design"
            )

        if not isinstance(records, list):
            records = [records]

        self.bins.append(records)


class PlasmidLibrary(Design):
    """A list of SeqRecord lists. Each nested list is combined into a Plasmid.

    Attributes:
        sets: Iterable of SeqRecord sets
        linear: Whether the SeqRecords of the design are linear
    """

    __name__ = "library"

    def __init__(self, sets: Iterable[List[SeqRecord]] = None, linear: bool = True):
        if sets:
            first_bin = list(sets)[0]

            if not isinstance(first_bin, list):
                raise ValueError(
                    "PlasmidLibrary requires a list of SeqRecord sets (lists).\n"
                    + f"First bin is: {type(first_bin)}"
                )

        self.sets: List[List[SeqRecord]] = list(sets) if sets else []
        self.linear = linear

    def __iter__(self) -> Iterable[List[SeqRecord]]:
        """Return an iterator over the sets"""

        if not self.sets:
            return iter([])
        return iter(self.sets)

    def __str__(self) -> str:
        """Return high level description of the number of records at each bin

        Example:
        ```[[4 x SeqRecord], [500 x SeqRecord], [3 x SeqRecord]]```

        Returns:
            high level overview of the number of SeqRecords in each bin
        """

        return f"[{len(self.sets)} x [SeqRecord]]"

    def __len__(self) -> int:
        """Return the number of library members."""

        return len(self.sets)

    def append(self, records: Union[SeqRecord, List[SeqRecord]]):
        """Add a new set or SeqRecord to the library

        Args:
            records: the records to add as a set
        """

        if not records:
            raise ValueError(
                f"Can only add SeqRecords or non-empty lists of SeqRecords to Design"
            )

        if not isinstance(records, list):
            records = [records]

        self.sets.append(records)
