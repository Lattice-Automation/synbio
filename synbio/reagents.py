"""Reagents."""


class Reagent:
    """A Reagent. Ex T4 Ligase, Buffer, etc.

    Keyword Arguments:
        name {str} -- the reagent's name (default: {""})
    """

    def __init__(self, name: str = ""):
        assert name

        self.name = name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)
