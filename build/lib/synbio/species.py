"""Lab species and organisms."""


class Species:
    """A Species. Ex E coli.

    Keyword Arguments:
        name {str} -- the species's name (default: {""})
    """

    def __init__(self, name: str = ""):
        self.name = name
