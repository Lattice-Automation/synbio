"""Reagents."""


class Reagent:
    """A Reagent. Ex T4 Ligase, Buffer, etc.

    Keyword Arguments:
        name {str} -- the reagent's name (default: {""})
    """

    def __init__(self, name: str = ""):
        self.name = name
