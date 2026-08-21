"""Stability invariants, including the one that a global agreement score hides."""

from __future__ import annotations

import pytest

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
