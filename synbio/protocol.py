"""A Protocol object for build assembly. Based on a Design and series of Steps."""

from collections import defaultdict
import logging
import os
import string
import unicodedata
from typing import Dict, List, Iterable, Tuple, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .containers import Container, Content, Fridge, Layout, content_id
from .design import Design
from .instructions import Transfer, Temperature, Instruction, to_txt
from .picklists import to_labcyte, to_tecan


class Step:
    """A single Step of a Protocol."""

    def __init__(self):
        self.transfers: List[Transfer] = []
        self.temps: List[Temperature] = []
        self.instructions: List[str] = []
        self.inputs: Dict[str, Tuple[Content, float]]

    def execute(self, protocol: "Protocol"):
        """Execute a step on a protocol, mutating containers and instructions."""

        raise NotImplementedError


class Protocol:
    """A Protocol for assembling DNA via steps in a laboratory.

    Arguments:
        name {str} -- the name of the protocol
        design {Design} -- the design specification for the build
        how {Union[Step, Iterable[Step]]} -- a single composite step for assembling the Design
    """

    def __init__(
        self,
        name: str = "",
        design: Design = Design(),
        how: Union[Step, Iterable[Step]] = None,
    ):
        self.name = name  # name of the protocol
        self.design = design  # the design specification
        self.steps: List[Step] = []  # list of steps for this assembly
        self.output: List[SeqRecord] = []

        # set steps from "how" if they were provided
        if how:
            if isinstance(how, Step):
                self.steps = [how]
            else:
                self.steps = list(how)

        # each of the below is set during a 'run()'
        self.containers: List[Container] = []
        self.instructions: List[Instruction] = []  # assembly instructions

        # map from each Instruction to the number of Layout before it
        self.plate_count = 0
        self.instruction_to_plate_count: Dict[Instruction, int] = defaultdict(int)

    def add(self, step: Step) -> "Protocol":
        """Add a step to the protocol.

        Arguments:
            step {Step} -- the Step to add to this protocol
        """

        if not isinstance(step, Step):
            raise TypeError

        self.steps.append(step)

        return self

    def run(self) -> "Protocol":
        """Run each step of the protcol. Build up the output records and instructions."""

        # all input records, each start out assigned to a single Fridge source
        records = self.design.get_all_records()
        self.containers = [Fridge(r) for r in records]  # first step, all fridge

        for step in self.steps:
            step.execute(self)

        return self

    @property
    def inputs(self) -> Dict[str, float]:
        """Return a map from protocol inputs to volumes in milliters.

        Returns:
            Dict[str, float] -- key, value map from input name to volume
        """

        # gather all transfers where the source was the Fridge
        id_to_volume: Dict[str, float] = defaultdict(float)
        for instruction in self.instructions:
            if not instruction.transfers:
                continue

            for transfer in instruction.transfers:
                if not isinstance(transfer.src, Fridge):
                    continue

                for content in transfer.src:
                    id_to_volume[content_id(content)] += transfer.volume

        # sort descending by volume
        id_volume_sorted = sorted(
            id_to_volume.items(), key=lambda x: x[1], reverse=True
        )
        return {cid: volume for cid, volume in id_volume_sorted}

    @property
    def outputs(self) -> List[SeqRecord]:
        """Gather the output SeqRecords from the final containers after all steps.

        Returns:
            List[SeqRecord] -- SeqRecords created from applying all steps
        """

        accumulator: List[SeqRecord] = []
        for container in self.containers:
            for content in container:
                if isinstance(content, SeqRecord):
                    accumulator.append(content)
        return accumulator

    def to_fasta(self, filename: str = "") -> int:
        """Write each output record to a FASTA file.

        Uses `SeqIO.write(records, filename, "fasta")`.

        Arguments:
            filename {str} -- the filename to write the FASTA file to

        Returns:
            int -- the number of records that were written
        """

        if not filename:
            filename = self._filename() + ".fasta"

        self._check_output()
        return SeqIO.write(self.outputs, filename, "fasta")

    def to_genbank(self, filename: str = "") -> int:
        """Write each output record to a Genbank file.

        Uses `SeqIO.write(records, filename, "genbank")`.

        TODO: add a separate flag for writing one GB per record

        Arguments:
            filename {str} -- the filename to write the Genbanks to

        Returns:
            int -- the number of records that were written
        """

        self._check_output()

        if not filename:
            filename = self._filename() + ".gb"

        # limit for id is 16 characters https://github.com/biopython/biopython/issues/747
        # shorten ids of those that exceed the limit
        output_renamed: List[SeqRecord] = []
        for i, output in enumerate(self.outputs):
            if len(output.id) > 16:
                output = output.upper()
                output.description = output.id
                output.id = "Seq" + str(i + 1)
                output_renamed.append(output)
            else:
                output_renamed.append(output)

        return SeqIO.write(output_renamed, filename, "genbank")

    def to_txt(self, filename: str = ""):
        """Write the protocol's instructions to a text file.

        Arguments:
            filename {str} -- the filename of the instructin file
        """

        if not filename:
            filename = self._filename() + ".csv"

        protocol_txt = to_txt(self.name, self.instructions)
        with open(filename, "w") as instruction_file:
            instruction_file.write(protocol_txt)

    def to_csv(self, filename: str = ""):
        """Write CSV file(s) describing the containers/Layout after each step.

        Rows of Layout/containers are written in CSV format with each step's
        name as its heading

        Example:
        ```py
        protocol.to_csv('instructions.csv')
        ```

        Arguments:
            filename {str} -- the name of the CSV(s) to write
        """

        if not filename:
            filename = self._filename() + ".csv"

        csv = ""
        row = 0
        for instruction in self.instructions:
            if not instruction.transfers:
                continue

            name = instruction.name
            row += 1
            csv += f"{name}:\n" if name else f"Setup step {row}:\n"
            csv += Layout.from_instruction(
                instruction,
                existing_plates=self.instruction_to_plate_count[instruction],
                log_volume=row == 1,
            ).to_csv()

        with open(filename, "w") as csvfile:
            csvfile.write(csv)

    def to_picklists(self, filename: str = "", platform: str = "tecan"):
        """Create picklists for robotic pipetting.

        For each step where there's plate to plate pipetting, create a
        robotic picklists. Steps where reagents or samples come from the Fridge
        are not written to a picklist right now.

        If there are multiple passable steps, each are saved with their
        index in their filename. Ex: picklist.1.gwl, picklist.2.gwl

        Keyword Arguments:
            filename {str} -- name of picklist file (default: {self.name})
            platform {str} -- picklist platform (default: {"tecan"})
        """

        picklist_generators = {
            "tecan": to_tecan,
            "hamilton": to_tecan,
            "labcyte": to_labcyte,
        }
        if platform not in picklist_generators:
            picklist_platforms = ", ".join(picklist_generators.keys())
            raise ValueError(
                f"'{platform}' is an unrecognized platform. Choose from: {picklist_platforms}"
            )

        self._check_output()

        if not filename:
            filename = self._filename()

            if platform == "tecan":
                filename += ".gwl"
            elif platform == "labcyte":
                filename += ".csv"

        # accumulate instructions from the protocol that are from plate to plate
        picklist_instructions: List[Instruction] = []
        for instruction in self.instructions:
            if not instruction.transfers:
                continue

            srcs = {t.src for t in instruction.transfers}
            dests = {t.dest for t in instruction.transfers}

            def no_fridge(containers: Iterable[Container]) -> bool:
                return all(not isinstance(s, Fridge) for s in containers)

            if no_fridge(srcs) and no_fridge(dests):
                picklist_instructions.append(instruction)

        if not picklist_instructions:
            raise RuntimeWarning(f"no picklist-capable steps in protocol")

        def picklist_name(index: int) -> str:
            if len(picklist_instructions) == 1:
                return filename

            fname, fext = os.path.splitext(filename)
            return fname + str(index + 1) + fext

        for i, instruction in enumerate(picklist_instructions):
            picklist = picklist_generators[platform](
                instruction, self.instruction_to_plate_count[instruction]
            )

            with open(picklist_name(i), "w") as picklist_file:
                picklist_file.write(picklist)

    def add_instruction(self, instruction: Instruction):
        """Add a single Instruction to this protocol's output.

        Instructions are generated by Steps. This is its own method
        because I'm also counting the number of "pipette-able" Layout after
        each step in instruction_to_plate_count.

        Arguments:
            instruction {Instruction} -- the instruction to add to Protocol
        """

        self.instructions.append(instruction)
        self.instruction_to_plate_count[instruction] = self.plate_count

        # add any pippete-able Layout
        if instruction.transfers:
            dest_containers = {t.dest for t in instruction.transfers}
            if all(not isinstance(c, Fridge) for c in dest_containers):
                self.plate_count += len(Layout.from_instruction(instruction))

    def _check_output(self):
        """Verify that the Protocol has steps and that they have been run.

        If there are no steps, or if there is no list of output Records
        in self.containers, there is nothing to write to a FASTA file.
        """
        if not self.steps and not self.outputs:
            raise RuntimeError(
                """No steps or containers specified in Protocol. Running `run()`"""
            )

        if self.steps and not self.outputs:
            logging.warning(
                """Steps specified but no output containers created. Running `run()`"""
            )
            self.run()

            if not self.outputs:
                raise RuntimeError(
                    """Failed to create valid assemblies after executing all steps."""
                )

    def _filename(self) -> str:
        """Gather a filename from this protocol's name.

        Convert this protocol's name to a name that's safe for a file.
        The below is from a Github Gist:
            https://gist.github.com/wassname/1393c4a57cfcbf03641dbc31886123b8
        """

        filename = self.name
        whitelist = "-_.() %s%s" % (string.ascii_letters, string.digits)
        filename = filename.replace(" ", "_")

        # keep only valid ascii chars
        cleaned_filename = (
            unicodedata.normalize("NFKD", filename).encode("ASCII", "ignore").decode()
        )

        # keep only whitelisted chars
        return "".join(c for c in cleaned_filename if c in whitelist)

    def __str__(self) -> str:
        """Return a human readable summary of the protocol.

        The python built in function str calls this. Should output a summary like:

        ```txt
        Combinatorial MoClo
            how: combinatorial
            design: [[15 x SeqRecord] [10 x SeqRecord] [2 x SeqRecord]]
            steps: []
        ```
        """

        name = self.name
        how = f"\thow: {self.design.__name__}"
        design = f"\tdesign: {str(self.design)}"

        return "\n".join([name, how, design])

    def __iter__(self) -> Iterable[Step]:
        """Iterate over the steps in this protocol."""

        return iter(self.steps)

    def __len__(self) -> int:
        """Return the number of steps in this protocol.

        Returns:
            int -- the number of steps
        """

        return len(self.steps)
