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
    """Fraction of cells whose raw label differs between two runs.

    **Almost always the wrong number to report.** Clustering algorithms number their
    clusters arbitrarily, so two runs that partition the cells identically can still
    disagree on every single raw label. Use :func:`matched_label_churn`, which resolves the
    renumbering first; this function is kept because it is what the matched version is
    measured against, and the gap between the two is itself informative.
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


def align(
    cells_a: Sequence[str],
    labels_a: Sequence[str],
    cells_b: Sequence[str],
    labels_b: Sequence[str],
) -> tuple[list[str], list[str], int]:
    """Restrict two labellings to the cells they share, in a common order.

    Configurations that filter differently keep different cells. Comparing their label
    vectors by position would score unrelated cells against each other and return a
    number that looks like instability but is really a bookkeeping error. Returns the two
    aligned label lists and the count of shared cells.
    """
    if len(cells_a) != len(labels_a) or len(cells_b) != len(labels_b):
        raise ValueError("cells and labels must be the same length")
    lookup_b = dict(zip(cells_b, labels_b, strict=True))
    shared = [cell for cell in cells_a if cell in lookup_b]
    if not shared:
        raise ValueError("the two runs share no cells")
    lookup_a = dict(zip(cells_a, labels_a, strict=True))
    return [lookup_a[c] for c in shared], [lookup_b[c] for c in shared], len(shared)


def match_labels(labels_a: Sequence[str], labels_b: Sequence[str]) -> dict[str, str]:
    """Best one-to-one mapping from ``labels_b``'s clusters onto ``labels_a``'s.

    Cluster numbering is arbitrary, so before two labellings can be compared cell by cell
    the renumbering has to be undone. This solves the assignment problem on the contingency
    table, which is the standard way -- greedy matching on the largest overlap can assign
    two of B's clusters to the same one of A's and then double-count the agreement.
    """
    from scipy.optimize import linear_sum_assignment

    categories_a = sorted(set(labels_a))
    categories_b = sorted(set(labels_b))
    index_a = {name: i for i, name in enumerate(categories_a)}
    index_b = {name: i for i, name in enumerate(categories_b)}

    overlap = [[0] * len(categories_b) for _ in categories_a]
    for a, b in zip(labels_a, labels_b, strict=True):
        overlap[index_a[a]][index_b[b]] += 1

    # linear_sum_assignment minimises, so negate to maximise overlap.
    costs = [[-count for count in row] for row in overlap]
    rows, columns = linear_sum_assignment(costs)
    return {
        categories_b[column]: categories_a[row]
        for row, column in zip(rows, columns, strict=True)
    }


def matched_label_churn(labels_a: Sequence[str], labels_b: Sequence[str]) -> float:
    """Fraction of cells that change cluster once arbitrary renumbering is undone.

    This is the number an analyst actually cares about: would this cell's call change? Any
    cluster in B with no counterpart in A -- which happens whenever the two runs find
    different numbers of clusters -- counts as changed for every cell in it, because there
    is no honest way to call those cells unchanged.
    """
    if len(labels_a) != len(labels_b):
        raise ValueError(f"length mismatch: {len(labels_a)} vs {len(labels_b)} labels")
    if not labels_a:
        raise ValueError("cannot measure churn over zero cells")
    mapping = match_labels(labels_a, labels_b)
    changed = sum(1 for a, b in zip(labels_a, labels_b, strict=True) if mapping.get(b) != a)
    return changed / len(labels_a)
