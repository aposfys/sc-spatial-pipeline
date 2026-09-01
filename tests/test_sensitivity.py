"""Stability invariants, including the one that a global agreement score hides."""

from __future__ import annotations

import pytest

from scspatial import sensitivity
from scspatial.sensitivity import adjusted_rand_index, label_churn, rare_population_stability


def test_identical_labellings_agree_perfectly() -> None:
    labels = ["A", "A", "B", "B", "C"]
    assert adjusted_rand_index(labels, labels) == pytest.approx(1.0)
    assert label_churn(labels, labels) == 0.0


def test_ari_ignores_relabelling_but_churn_does_not() -> None:
    """Renaming every cluster is not a biological change; reassigning cells is."""
    original = ["A", "A", "B", "B"]
    renamed = ["X", "X", "Y", "Y"]
    assert adjusted_rand_index(original, renamed) == pytest.approx(1.0)
    assert label_churn(original, renamed) == 1.0


def test_random_relabelling_scores_near_zero() -> None:
    a = ["A", "A", "A", "A", "B", "B", "B", "B"]
    b = ["A", "B", "A", "B", "A", "B", "A", "B"]
    assert adjusted_rand_index(a, b) == pytest.approx(0.0, abs=0.2)


def test_a_rare_population_can_dissolve_while_ari_stays_high() -> None:
    """The central claim of this repo, asserted as a test rather than a README sentence."""
    abundant = ["A"] * 100
    rare_before = ["R"] * 5  # 5 / 105 sits just under the 5% rarity threshold
    rare_after = ["A", "A", "A", "B", "B"]

    before = abundant + rare_before
    after = abundant + rare_after

    assert adjusted_rand_index(before, after) < 1.0
    retention = rare_population_stability(before, after)
    assert retention["R"] == pytest.approx(0.6)
    assert "A" not in retention


def test_abundant_populations_are_not_reported_as_rare() -> None:
    labels = ["A"] * 50 + ["B"] * 50
    assert rare_population_stability(labels, labels) == {}


def test_length_mismatch_is_fatal_in_every_metric() -> None:
    for metric in (adjusted_rand_index, label_churn, rare_population_stability):
        with pytest.raises(ValueError):
            metric(["A", "B"], ["A"])


def test_empty_input_is_an_error_not_perfect_stability() -> None:
    with pytest.raises(ValueError):
        label_churn([], [])


def test_matched_churn_is_zero_for_a_pure_relabelling():
    """The bug this metric exists to avoid: renumbering is not instability."""
    a = ["1", "1", "2", "2", "3", "3"]
    b = ["7", "7", "9", "9", "4", "4"]  # identical partition, different numbers
    assert sensitivity.label_churn(a, b) == 1.0  # raw churn is fooled
    assert sensitivity.matched_label_churn(a, b) == 0.0  # matched churn is not
    assert sensitivity.adjusted_rand_index(a, b) == pytest.approx(1.0)


def test_matched_churn_counts_genuine_moves():
    a = ["1", "1", "1", "2", "2", "2"]
    b = ["1", "1", "2", "2", "2", "2"]
    # One cell genuinely moved from the first cluster to the second.
    assert sensitivity.matched_label_churn(a, b) == pytest.approx(1 / 6)


def test_match_labels_is_one_to_one():
    """Greedy matching can map two of B's clusters onto one of A's and double-count."""
    a = ["1"] * 10 + ["2"] * 10
    b = ["x"] * 9 + ["y"] + ["y"] * 9 + ["x"]
    mapping = sensitivity.match_labels(a, b)
    assert sorted(mapping) == ["x", "y"]
    assert len(set(mapping.values())) == 2


def test_align_restricts_to_shared_cells():
    labels_a, labels_b, shared = sensitivity.align(
        ["c1", "c2", "c3"],
        ["a", "b", "c"],
        ["c2", "c3", "c4"],
        ["B", "C", "D"],
    )
    assert shared == 2
    assert labels_a == ["b", "c"]
    assert labels_b == ["B", "C"]


def test_align_refuses_when_nothing_is_shared():
    with pytest.raises(ValueError, match="share no cells"):
        sensitivity.align(["c1"], ["a"], ["c9"], ["z"])


def test_align_refuses_mismatched_inputs():
    with pytest.raises(ValueError, match="same length"):
        sensitivity.align(["c1", "c2"], ["a"], ["c1"], ["a"])
