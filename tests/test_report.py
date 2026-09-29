"""The stability table, measured on hand-built runs rather than on a real grid."""

from __future__ import annotations

import sys

import pytest

from scspatial import pipeline, report
from scspatial.configs import REFERENCE, one_at_a_time
from scspatial.pipeline import RunResult


def _run(key: str, labels: list[str], pairs: list[tuple[str, str]]) -> RunResult:
    cells = [f"spot{i}" for i in range(len(labels))]
    return RunResult(
        key=key,
        cells=cells,
        labels=labels,
        n_cells=len(labels),
        n_clusters=len(set(labels)),
        enriched_pairs=pairs,
        seconds=0.0,
    )


def test_compare_scores_a_renumbered_run_as_identical() -> None:
    """The failure the conclusion column used to have: renumbering read as disagreement."""
    reference = _run("ref", ["0", "0", "1", "1", "2", "2"], [("0", "1")])
    renumbered = _run("other", ["5", "5", "3", "3", "4", "4"], [("3", "5")])
    config = REFERENCE.as_dict()
    row = report.compare(reference, renumbered, config, config)
    assert row.ari == pytest.approx(1.0)
    assert row.churn == 0.0
    assert row.worst_rare_retention == 1.0
    assert row.conclusion_jaccard == 1.0
    assert row.raw_conclusion_jaccard == 0.0


def test_render_reports_the_true_median() -> None:
    configs = [c.as_dict() for c in one_at_a_time()]
    rows = [
        {
            "changed_axis": report.changed_axis(configs[0], config),
            "n_clusters": 3,
            "ari": ari,
            "churn": 0.0,
            "worst_rare_retention": 1.0,
            "conclusion_jaccard": overlap,
        }
        for config, ari, overlap in zip(
            configs,
            [1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
            [1.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            strict=True,
        )
    ]
    findings = {"dataset": "toy", "reference": {"n_cells": 6}, "comparisons": rows}
    text = report.render(findings)
    # Eight deviations, so the median is the mean of the fourth and fifth values.
    assert "8 deviations that changed the clustering" in text
    assert "median ARI is **0.450**" in text
    assert "median conclusion overlap is **0.35**" in text
    assert "upper bound" not in text


def test_render_leaves_a_deviation_that_changed_nothing_out_of_the_median() -> None:
    configs = [c.as_dict() for c in one_at_a_time()]
    rows = [
        {
            "changed_axis": report.changed_axis(configs[0], config),
            "n_clusters": 3,
            "ari": ari,
            "churn": 0.0 if ari == 1.0 else 0.1,
            "worst_rare_retention": 1.0,
            "conclusion_jaccard": 1.0 if ari == 1.0 else 0.5,
        }
        for config, ari in zip(
            configs, [1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0], strict=True
        )
    ]
    text = report.render({"dataset": "toy", "reference": {"n_cells": 6}, "comparisons": rows})
    assert "7 deviations that changed the clustering" in text
    assert "median ARI is **0.400**" in text
    assert "no spot changed cluster: filter_order = normalise_then_filter" in text


def test_build_runs_stores_shared_barcodes_once() -> None:
    first = _run("a", ["0", "1"], [])
    second = _run("b", ["1", "0"], [("0", "1")])
    runs = report.build_runs([first, second])
    assert runs["cells"] == ["spot0", "spot1"]
    assert "cells" not in runs["runs"]["b"]
    assert runs["runs"]["b"]["enriched_pairs"] == [["0", "1"]]


def test_counts_check_accepts_counts_and_rejects_log_normalised_values() -> None:
    import numpy as np
    from scipy import sparse

    assert pipeline.looks_like_counts(sparse.csr_matrix(np.array([[0, 3], [1, 0]])))
    assert not pipeline.looks_like_counts(sparse.csr_matrix(np.log1p([[0, 3], [1, 0]])))


def test_missing_stack_exits_with_an_install_hint(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "scanpy", None)
    with pytest.raises(SystemExit, match=r"pip install -e"):
        pipeline.require_stack()
