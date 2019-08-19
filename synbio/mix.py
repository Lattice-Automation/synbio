"""Mixtures of Contents."""

import inspect
import logging
from typing import Dict, Iterable, List, Tuple, Set, Any

from .containers import Content
from .reagents import Reagent
from .species import Species


class Mix:
    """An assembly Mix to specify how Reagents/Species/SeqRecords should come together

    'mix' is a parameter that accepts a Dict mapping each of the mix's ingredients
    to the volume in microliters of that ingredient to the mixture container.

    Below is an example for NEB's GoldenGate Assembly Mix (E1600) where
    there is 2 uL per DNA sample, 4 uL for assembly master mix, and the well
    is filled to 20 uL

    ```python
    Mix(
        mix={SeqRecord: 2.0, Reagent("master mix"): 4.0},
        fill_with=Reagent("water"),
        fill_to=20.0,
    )
    ```

    Keyword Args:
        mix: map from mix content to the
            volume required per well. Mix content can be a type (like SeqRecord)
            or an actual reagent (like Reagent("water")) (default: {None})
        fill_with: Content to fill rest of container to
            fill_to with (default: {None})
        fill_to: the volume to fill the container to (default: {0})
    """

    def __init__(
        self,
        mix: Dict[Any, float] = None,
        fill_with: Content = None,
        fill_to: float = 0,
    ):
        if fill_with and not fill_to:
            raise ValueError(
                f"Cannot specify what to 'fill_with' without specifying container's 'fill_to'"
            )

        if not fill_with and fill_to:
            raise ValueError(
                f"Cannot specify what to 'fill_to' without 'fill_with' Content (Reagent, Species or SeqRecord)"
            )

        self.mix = mix or {}
        self.fill_with = fill_with
        self.fill_to = fill_to

    def __call__(
        self, contents: Iterable[Content]
    ) -> Tuple[List[Content], List[float]]:
        """Call Mix on a list of contents to generate container's contents/volumes

        Figure out volume needed for each item in contents.
        Return those contents plus the additional contents needed in the mix (self.mix)
        Also return a second list with the volumes needed for each content in the
        first returned list

        Args:
            contents: the contents of a well to add to a mix

        Returns:
            Tuple[List[Content], List[float]] -- list of container contents and
                a list of volumes for each content item
        """

        contents_out: List[Content] = []
        volumes: List[float] = []

        seen: Set[Any] = set()

        for content in contents:
            if (
                isinstance(content, Reagent) or isinstance(content, Species)
            ) and content in self.mix:
                # reagent was explicitly specified
                contents_out.append(content)
                volumes.append(self.mix[content])
                seen.add(content)
            elif type(content) in self.mix:
                # example is SeqRecord
                contents_out.append(content)
                volumes.append(self.mix[type(content)])
                seen.add(type(content))
            elif any(
                isinstance(content, t) for t in self.mix.keys() if inspect.isclass(t)
            ):
                # example is RestrictionType class
                class_type = next(
                    t
                    for t in self.mix.keys()
                    if inspect.isclass(t) and isinstance(content, t)
                )
                contents_out.append(content)
                volumes.append(self.mix[class_type])
                seen.add(class_type)
            else:
                logging.warning(f"Content {content} not found in mix")

        for content, volume in self.mix.items():
            if content in seen or inspect.isclass(content):
                continue

            contents_out.append(content)
            volumes.append(volume)

        if self.fill_to and self.fill_with:
            volume_total = sum(volumes)
            volume_remaining = max([self.fill_to - volume_total, 0.0])
            contents_out.append(self.fill_with)
            volumes.append(volume_remaining)

        return (contents_out, volumes)
