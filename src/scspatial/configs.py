"""The grid of defensible preprocessing choices.

Every option here is one a competent analyst might pick and would not feel obliged to
justify in a methods section. That is the criterion for inclusion: not "plausible", but
"would pass review unremarked". The whole question is how much the answer moves inside that
space.

**Segmentation is absent from this grid, and that is a property of the data, not an
oversight.** Spot-based Visium has no segmentation step -- the spots are a fixed grid laid
down by the assay, and no analyst choice changes which transcripts land in which spot. The
design notes name segmentation as the stage whose errors most often become discoveries, and
that remains true; it is simply not a choice this modality offers. On imaging-based data it
would be the first axis in the grid.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterator
from dataclasses import dataclass

#: Normalisation strategies. `none` is included because analysts do skip it on already
#: size-factor-corrected data, not because it is recommended.
NORMALISATIONS = ("cpm_log1p", "cp10k_log1p", "none")

#: Highly variable gene selection. Both flavours are standard and neither is stated more
#: often than the other.
HVG_METHODS = ("seurat", "cell_ranger")
HVG_COUNTS = (2000, 4000)

#: Neighbourhood size for the kNN graph that clustering runs on.
NEIGHBOURS = (10, 15, 30)

#: Leiden resolution. The single most consequential unreported number in the field.
RESOLUTIONS = (0.5, 1.0)

#: Whether filtering happens before or after normalisation. Rarely stated, not neutral.
FILTER_ORDERS = ("filter_then_normalise", "normalise_then_filter")


@dataclass(frozen=True)
class Config:
    """One point in the grid."""

    normalisation: str
    hvg_method: str
    n_hvg: int
    n_neighbours: int
    resolution: float
    filter_order: str
    seed: int = 0

    @property
    def key(self) -> str:
        """Short stable identifier, used as the column name in the stability table."""
        return (
            f"{self.normalisation}|{self.hvg_method}{self.n_hvg}"
            f"|k{self.n_neighbours}|r{self.resolution}|{self.filter_order[:3]}"
        )

    def axes(self) -> dict:
        """The grid axes only, without the seed. Used to vary one axis at a time."""
        return {
            "normalisation": self.normalisation,
            "hvg_method": self.hvg_method,
            "n_hvg": self.n_hvg,
            "n_neighbours": self.n_neighbours,
            "resolution": self.resolution,
            "filter_order": self.filter_order,
        }

    def as_dict(self) -> dict:
        return {
            "key": self.key,
            "normalisation": self.normalisation,
            "hvg_method": self.hvg_method,
            "n_hvg": self.n_hvg,
            "n_neighbours": self.n_neighbours,
            "resolution": self.resolution,
            "filter_order": self.filter_order,
            "seed": self.seed,
        }


#: The configuration everything else is compared against. Not "the right answer" -- just
#: the one a tutorial would produce, which is what most published pipelines are.
REFERENCE = Config(
    normalisation="cp10k_log1p",
    hvg_method="seurat",
    n_hvg=2000,
    n_neighbours=15,
    resolution=1.0,
    filter_order="filter_then_normalise",
)


def full_grid(seed: int = 0) -> list[Config]:
    """Every combination. 144 configurations."""
    return [
        Config(
            normalisation=normalisation,
            hvg_method=hvg_method,
            n_hvg=n_hvg,
            n_neighbours=n_neighbours,
            resolution=resolution,
            filter_order=filter_order,
            seed=seed,
        )
        for normalisation, hvg_method, n_hvg, n_neighbours, resolution, filter_order in itertools.product(
            NORMALISATIONS, HVG_METHODS, HVG_COUNTS, NEIGHBOURS, RESOLUTIONS, FILTER_ORDERS
        )
    ]


def one_at_a_time(seed: int = 0) -> list[Config]:
    """The reference, plus one configuration per single deviation from it.

    Cheaper than the full grid and answers a different, sharper question: which *single*
    choice moves the answer most. A full grid tells you the total spread; this tells you
    where it comes from.
    """
    configs = [REFERENCE]
    axes: dict[str, tuple] = {
        "normalisation": NORMALISATIONS,
        "hvg_method": HVG_METHODS,
        "n_hvg": HVG_COUNTS,
        "n_neighbours": NEIGHBOURS,
        "resolution": RESOLUTIONS,
        "filter_order": FILTER_ORDERS,
    }
    for axis, values in axes.items():
        for value in values:
            if getattr(REFERENCE, axis) == value:
                continue
            configs.append(Config(**{**REFERENCE.axes(), axis: value, "seed": seed}))
    return configs


def iter_grid(mode: str = "one_at_a_time", seed: int = 0) -> Iterator[Config]:
    """Grid by name. ``full`` is 144 configurations; ``one_at_a_time`` is 9."""
    if mode == "full":
        yield from full_grid(seed)
    elif mode == "one_at_a_time":
        yield from one_at_a_time(seed)
    else:
        raise ValueError(f"unknown grid mode {mode!r}; expected 'full' or 'one_at_a_time'")
