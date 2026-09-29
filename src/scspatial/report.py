"""The stability table.

Everything is measured against the reference configuration, the scanpy tutorial path from raw
counts. The question is not which configuration is best. It is how far the answer moves when
one choice changes that nobody reports.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path

from scspatial.sensitivity import (
    adjusted_rand_index,
    align,
    conclusion_overlap,
    label_churn,
    matched_label_churn,
    rare_population_stability,
)


@dataclass
class Comparison:
    """One configuration measured against the reference."""

    key: str
    changed_axis: str
    n_clusters: int
    shared_cells: int
    ari: float
    #: Raw churn, which counts arbitrary cluster renumbering as change. Reported only so
    #: the gap against the matched figure is visible.
    raw_churn: float
    #: Churn after the renumbering is resolved. This is the honest one.
    churn: float
    #: Worst per-population retention among rare populations. The number ARI hides.
    worst_rare_retention: float
    n_rare_populations: int
    #: Jaccard of the enriched cluster-pair sets, after this run's clusters are matched
    #: one-to-one onto the reference's.
    conclusion_jaccard: float
    #: The same Jaccard on raw cluster numbers, without matching. Kept, like raw churn,
    #: so the size of the numbering artefact stays visible.
    raw_conclusion_jaccard: float
    seconds: float


def changed_axis(reference_config: dict, other: dict) -> str:
    """Which single axis differs between two configurations."""
    differing = [
        name
        for name in (
            "normalisation",
            "hvg_method",
            "n_hvg",
            "n_neighbours",
            "resolution",
            "filter_order",
        )
        if reference_config[name] != other[name]
    ]
    if not differing:
        return "(reference)"
    if len(differing) == 1:
        name = differing[0]
        return f"{name} = {other[name]}"
    return f"{len(differing)} axes"


def _raw_jaccard(pairs_a, pairs_b) -> float:
    a = {tuple(sorted(pair)) for pair in pairs_a}
    b = {tuple(sorted(pair)) for pair in pairs_b}
    return len(a & b) / len(a | b) if a | b else 1.0


def compare(reference, other, reference_config: dict, other_config: dict) -> Comparison:
    """Measure one run against the reference."""
    labels_a, labels_b, shared = align(
        reference.cells, reference.labels, other.cells, other.labels
    )
    retention = rare_population_stability(labels_a, labels_b)

    return Comparison(
        key=other.key,
        changed_axis=changed_axis(reference_config, other_config),
        n_clusters=other.n_clusters,
        shared_cells=shared,
        ari=adjusted_rand_index(labels_a, labels_b),
        raw_churn=label_churn(labels_a, labels_b),
        churn=matched_label_churn(labels_a, labels_b),
        worst_rare_retention=min(retention.values()) if retention else 1.0,
        n_rare_populations=len(retention),
        conclusion_jaccard=conclusion_overlap(
            labels_a, labels_b, reference.enriched_pairs, other.enriched_pairs
        ),
        raw_conclusion_jaccard=_raw_jaccard(reference.enriched_pairs, other.enriched_pairs),
        seconds=other.seconds,
    )


def render(findings: dict) -> str:
    """Render the findings as Markdown."""
    rows = findings["comparisons"]
    reference = findings["reference"]
    lines: list[str] = []

    lines.append("# Results\n")
    lines.append(
        f"{findings['dataset']}, {reference['n_cells']:,} spots. "
        f"{len(rows)} configurations, each differing from the reference in one choice.\n"
    )
    lines.append(
        "The reference is a standard scanpy path from the raw counts (normalise to "
        "10,000 counts per spot, log1p, 2,000 Seurat HVGs, 15 neighbours, Leiden at "
        "resolution 1.0). It is not a claim about what is correct.\n"
    )

    lines.append(
        "| Changed choice | Clusters | ARI | Cell churn | Worst rare retention | Conclusion overlap |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: |")
    for row in rows:
        lines.append(
            f"| {row['changed_axis']} | {row['n_clusters']} | {row['ari']:.3f} "
            f"| {row['churn']:.1%} | {row['worst_rare_retention']:.1%} "
            f"| {row['conclusion_jaccard']:.2f} |"
        )
    lines.append("")

    non_reference = [row for row in rows if row["changed_axis"] != "(reference)"]
    # A deviation that relabels no spot measures nothing, so it is left out of the medians.
    changed = [row for row in non_reference if row["ari"] < 1.0 or row["churn"] > 0.0]
    if non_reference:
        worst_ari = min(non_reference, key=lambda row: row["ari"])
        worst_churn = max(non_reference, key=lambda row: row["churn"])
        worst_rare = min(non_reference, key=lambda row: row["worst_rare_retention"])
        worst_conclusion = min(non_reference, key=lambda row: row["conclusion_jaccard"])

        lines.append("## What moved\n")
        lines.append(
            f"- **Lowest global agreement:** `{worst_ari['changed_axis']}` at ARI "
            f"{worst_ari['ari']:.3f}."
        )
        lines.append(
            f"- **Most cells relabelled:** `{worst_churn['changed_axis']}` moved "
            f"{worst_churn['churn']:.1%} of spots to a different cluster."
        )
        worst = worst_rare["worst_rare_retention"]
        tied = [row for row in non_reference if row["worst_rare_retention"] == worst]
        if len(tied) > 1:
            lines.append(
                f"- **Worst rare-population retention:** {len(tied)} of "
                f"{len(non_reference)} deviations leave at least one rare reference cluster "
                f"with {worst:.1%} of its spots in its matched counterpart."
            )
        else:
            lines.append(
                f"- **Worst rare-population retention:** `{worst_rare['changed_axis']}`. "
                f"One rare reference cluster kept {worst:.1%} of its spots in its matched "
                f"counterpart."
            )
        lines.append(
            f"- **Least stable spatial conclusion:** `{worst_conclusion['changed_axis']}` "
            f"shares only {worst_conclusion['conclusion_jaccard']:.0%} of its "
            f"neighbourhood-enrichment pairs with the reference.\n"
        )

    if changed:
        unchanged = [row["changed_axis"] for row in non_reference if row not in changed]
        median_ari = statistics.median(row["ari"] for row in changed)
        median_overlap = statistics.median(row["conclusion_jaccard"] for row in changed)
        lines.append(
            f"Across the {len(changed)} deviations that changed the clustering the median "
            f"ARI is **{median_ari:.3f}** and the median conclusion overlap is "
            f"**{median_overlap:.2f}**."
            + (
                f" Left out because no spot changed cluster: {', '.join(unchanged)}."
                if unchanged
                else ""
            )
            + "\n"
        )
        raw = [row.get("raw_conclusion_jaccard") for row in changed]
        if all(value is not None for value in raw):
            lines.append(
                "Compared on raw cluster numbers, without matching, the same conclusion "
                f"overlap would read {min(raw):.2f} to {max(raw):.2f}.\n"
            )

    lines.append("## Reading these columns\n")
    lines.append(
        "- **ARI** is the adjusted Rand index against the reference labels. It is "
        "dominated by the abundant clusters.\n"
        "- **Cell churn** is the fraction of spots whose cluster changes once the "
        "arbitrary cluster numbering is resolved by optimal one-to-one assignment. Raw "
        "churn, which skips that step, runs far higher and is mostly an artefact of "
        "numbering.\n"
        "- **Worst rare retention** is, for the rare reference clusters (under 5% of "
        "spots), the lowest fraction of a cluster's spots that land in its matched "
        "counterpart. A rare cluster absorbed into an abundant one scores 0.\n"
        "- **Conclusion overlap** is the Jaccard of the cluster pairs called "
        "significantly co-located by neighbourhood enrichment (z >= 2, self-pairs "
        "excluded), after each run's clusters are matched onto the reference's by the "
        "same assignment. A pair involving a cluster with no counterpart cannot match.\n"
    )
    versions = findings.get("versions")
    if versions:
        lines.append(
            "Run with "
            + ", ".join(f"{name} {value}" for name, value in versions.items())
            + ".\n"
        )
    return "\n".join(lines)


def write(findings_path: Path, out_path: Path) -> Path:
    findings = json.loads(findings_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render(findings))
    return out_path


def build_runs(results) -> dict:
    """Every run's spot labels and enriched pairs, keyed by configuration.

    Spot barcodes are stored once. A run that kept a different set of spots carries its
    own ``cells`` list.
    """
    cells = results[0].cells
    runs: dict[str, dict] = {}
    for result in results:
        run: dict = {
            "labels": result.labels,
            "enriched_pairs": [list(pair) for pair in result.enriched_pairs],
        }
        if result.cells != cells:
            run["cells"] = result.cells
        runs[result.key] = run
    return {"cells": cells, "runs": runs}


def build_findings(results, configs, dataset: str, versions: dict | None = None) -> dict:
    """Assemble the findings document from a list of runs."""
    reference_result = results[0]
    reference_config = configs[0]
    comparisons = [
        asdict(compare(reference_result, result, reference_config, config))
        for result, config in zip(results, configs, strict=True)
    ]
    return {
        "dataset": dataset,
        "reference": {
            "key": reference_result.key,
            "n_cells": reference_result.n_cells,
            "n_clusters": reference_result.n_clusters,
            "config": reference_config,
        },
        "comparisons": comparisons,
        "versions": versions or {},
    }
