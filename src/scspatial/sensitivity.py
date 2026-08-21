"""Agreement between two runs of the same analysis under different defensible choices.

Global agreement is not enough. The adjusted Rand index is dominated by abundant cell
types, so a rare population -- usually the interesting one -- can be entirely reassigned
while ARI stays near 1.0. Both numbers are computed here and both are reported.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from math import comb


def adjusted_rand_index(labels_a: Sequence[str], labels_b: Sequence[str]) -> float:
    """Adjusted Rand index between two labellings of the same cells.

    Implemented directly rather than imported so that the stability layer has no
    dependency on the analysis stack it is measuring.
    """
    if len(labels_a) != len(labels_b):
        raise ValueError(f"length mismatch: {len(labels_a)} vs {len(labels_b)} labels")
    n = len(labels_a)
    if n < 2:
        raise ValueError("adjusted Rand index needs at least two cells")

    contingency = Counter(zip(labels_a, labels_b, strict=True))
    counts_a = Counter(labels_a)
    counts_b = Counter(labels_b)

    index = sum(comb(count, 2) for count in contingency.values())
    sum_a = sum(comb(count, 2) for count in counts_a.values())
    sum_b = sum(comb(count, 2) for count in counts_b.values())
    total = comb(n, 2)

    expected = sum_a * sum_b / total
    maximum = (sum_a + sum_b) / 2
    if maximum == expected:
        return 1.0
    return (index - expected) / (maximum - expected)


def label_churn(labels_a: Sequence[str], labels_b: Sequence[str]) -> float:
    """Fraction of cells whose label differs between two runs.

    Unlike ARI this is not invariant to relabelling, which is the point: it answers the
    question an analyst actually asks -- "would my cell-type call for this cell change?"
    """
    if len(labels_a) != len(labels_b):
        raise ValueError(f"length mismatch: {len(labels_a)} vs {len(labels_b)} labels")
    if not labels_a:
        raise ValueError("cannot measure churn over zero cells")
    changed = sum(1 for a, b in zip(labels_a, labels_b, strict=True) if a != b)
    return changed / len(labels_a)


def rare_population_stability(
    labels_a: Sequence[str],
    labels_b: Sequence[str],
    *,
    rare_threshold: float = 0.05,
) -> dict[str, float]:
    """Per-population retention, restricted to populations below ``rare_threshold``.

    Returns, for each rare population in the first labelling, the fraction of its cells
    that remain together in the second. This is the number a global agreement score hides.
    """
    if len(labels_a) != len(labels_b):
        raise ValueError(f"length mismatch: {len(labels_a)} vs {len(labels_b)} labels")
    if not labels_a:
        raise ValueError("cannot measure stability over zero cells")

    counts = Counter(labels_a)
    total = len(labels_a)
    retention: dict[str, float] = {}
    for population, size in counts.items():
        if size / total >= rare_threshold:
            continue
        partners = Counter(
            b for a, b in zip(labels_a, labels_b, strict=True) if a == population
        )
        retention[population] = max(partners.values()) / size
    return retention
