"""Sequence level utility functions."""

from typing import Set


def mutate(seq: str, edit_distance: int = 1) -> Set[str]:
    """Constructs all sequence within some edit distance of a starting seq

    Args:
        seq: The sequence to create mutations of

    Keyword Args:
        edit_distance: The max edit distance from the starting seq

    Returns:
        A set of possible off sequences within edit distance
    """

    assert edit_distance < len(seq)

    mutants: Set[str] = set([seq])
    for _ in range(edit_distance):
        mutants_new: Set[str] = set()
        for mutant in mutants:
            mutants_new.add(mutant)
            for i in range(len(mutant)):
                for base in "AGCT":
                    mutants_new.add(mutant[:i] + base + mutant[i + 1 :])
        mutants = mutants_new
    return mutants
