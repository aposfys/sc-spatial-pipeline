"""The stability table.

Everything is measured against the reference configuration -- the one a tutorial produces,
which is what most published pipelines are. The question is not which configuration is best;
it is how far the answer moves when you change a choice nobody reports.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from scspatial.sensitivity import (
    adjusted_rand_index,
    align,
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
    #: Jaccard of the enriched cluster-pair sets. Cluster labels are not comparable across
    #: runs, so this is an upper bound on agreement rather than an exact match.
    conclusion_jaccard: float
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


def compare(reference, other, reference_config: dict, other_config: dict) -> Comparison:
    """Measure one run against the reference."""
    labels_a, labels_b, shared = align(
        reference.cells, reference.labels, other.cells, other.labels
    )
    retention = rare_population_stability(labels_a, labels_b)
    reference_pairs = {tuple(sorted(pair)) for pair in reference.enriched_pairs}
    other_pairs = {tuple(sorted(pair)) for pair in other.enriched_pairs}
    union = reference_pairs | other_pairs

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
        conclusion_jaccard=(len(reference_pairs & other_pairs) / len(union)) if union else 1.0,
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
        "The reference is the configuration a tutorial produces. It is not a claim about "
        "what is correct -- it is what most published pipelines are.\n"
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
            f"{worst_churn['churn']:.1%} of cells to a different cluster."
        )
        lines.append(
            f"- **Worst rare-population retention:** `{worst_rare['changed_axis']}` — the "
            f"least stable rare population kept only "
            f"{worst_rare['worst_rare_retention']:.1%} of its cells together."
        )
        lines.append(
            f"- **Least stable spatial conclusion:** `{worst_conclusion['changed_axis']}` "
            f"shares only {worst_conclusion['conclusion_jaccard']:.0%} of its "
            f"neighbourhood-enrichment pairs with the reference.\n"
        )

        median_ari = sorted(row["ari"] for row in non_reference)[len(non_reference) // 2]
        lines.append(
            f"Median ARI across single-choice deviations is **{median_ari:.3f}**. "
            "Global agreement that high is the number usually reported; the rare-population "
            "and conclusion columns are where the instability actually lives.\n"
        )

    lines.append("## Reading these columns\n")
    lines.append(
        "- **ARI** is dominated by abundant populations. It is the reassuring number.\n"
        "- **Cell churn** answers the question an analyst asks -- whether *this cell's* "
        "call would change -- and is computed after resolving arbitrary cluster "
        "renumbering by optimal assignment. Raw churn, which does not, runs far higher "
        "and is almost entirely an artefact of numbering.\n"
        "- **Worst rare retention** is the number ARI hides. A rare population can be "
        "entirely reassigned while ARI barely moves.\n"
        "- **Conclusion overlap** is the Jaccard of the significantly co-located cluster "
        "pairs. Cluster identities are not comparable across runs, so this is an upper "
        "bound on agreement, not an exact match -- it can only overstate stability.\n"
    )
    return "\n".join(lines)


def write(findings_path: Path, out_path: Path) -> Path:
    findings = json.loads(findings_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render(findings))
    return out_path


def build_findings(results, configs, dataset: str) -> dict:
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
    }
