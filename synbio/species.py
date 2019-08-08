"""Lab species and organisms."""


class Species:
    """A Species. Ex E coli.

    Keyword Arguments:
        name {str} -- the species's name (default: {""})
    """

    def __init__(self, name: str = ""):
        assert name

        self.name = name

    def __hash__(self):
        """Hash species."""

        return hash(self.name)

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)
