"""Designs are build specifications that return lists of SeqRecord lists."""

from typing import List, Union, Dict

from Bio.SeqRecord import SeqRecord

from .containers import content_id


class Design:
    """A SynBio design specification."""

    __name__ = "design"

    def __iter__(self):
        raise NotImplementedError

    def append(self, records: Union[SeqRecord, List[SeqRecord]]):
        """Add to an existing design."""

        raise NotImplementedError

    def get_all_records(self) -> List[SeqRecord]:
        """Return all unique SeqRecords across all designs.

        Returns:
            {List[SeqRecord]} -- a unique list of SeqRecords used across designs
        """

        id_to_record: Dict[str, SeqRecord] = {}
        for record_set in self:
            for record in record_set:
                record_id = content_id(record)
                id_to_record[record_id] = record
        return list(id_to_record.values())


class Plasmid(Design):
    """Create a plasmid from a list of SeqRecords

    Arguments:
        records {List[SeqRecord]} -- SeqRecords to concatenate into a plasmid
    """

    __name__ = "plasmid"

    def __init__(self, records: List[SeqRecord] = None):
        self.records = records if records else []

    def __iter__(self):
        return iter([self.records])

    def __str__(self) -> str:
        """Return the ids of each piece of DNA in this plasmid

        Returns:
            [str] -- high level description of this plasmid's design
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

        Arguments:
            record {SeqRecord} -- the SeqRecord to add to the plasmid specification
        """

        if not isinstance(records, list):
            records = [records]
        for record in records:
            self.records.append(record)


class Combinatorial(Design):
    """Create a Plasmid design specification from a list of SeqRecords

    Arguments:
        bins {List[List[SeqRecord]]} -- list of SeqRecord bins
    """

    __name__ = "combinatorial"

    def __init__(self, bins: List[List[SeqRecord]] = None):
        self.bins = bins if bins else []

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
            str -- high level overview of the number of SeqRecords in each bin
        """

        return (
            f"[{', '.join(['[' + str(len(b)) + ' x SeqRecord]' for b in self.bins])}]"
        )

    def __len__(self) -> int:
        """Return the length of the bins."""

        return len(self.bins)

    def append(self, records: Union[SeqRecord, List[SeqRecord]]):
        """Add the SeqRecords as a bin for combinatorial design

        Arguments:
            record {Union[SeqRecord, List[SeqRecord]]} -- the records to add as a bin
        """

        if not isinstance(records, list):
            records = [records]

        self.bins.append(records)
