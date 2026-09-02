# sc-spatial-pipeline
How much of a spatial transcriptomics result is the pipeline rather than the tissue?

[![CI](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A spatial analysis has several defensible choices at every stage, and the
published figure shows one path through that tree. This measures the width of the
tree: same dataset (squidpy Visium H&E, 2,688 spots), same question, one
preprocessing choice changed at a time.

```
make install
make data       # squidpy Visium H&E, ~330 MB, cached
make grid       # 9 configurations, ~5 s each
make analysis   # stability table -> results/RESULTS.md
make test       # 13 tests, no analysis stack needed
```

### Two metrics, opposite verdicts

| Changed choice | Clusters | ARI | Cell churn | Worst rare retention | **Conclusion overlap** |
| --- | ---: | ---: | ---: | ---: | ---: |
| (reference) | 18 | 1.000 | 0.0% | 100.0% | 1.00 |
| normalisation = cpm_log1p | 18 | 0.778 | 14.0% | 53.8% | **0.00** |
| normalisation = none | 19 | 0.810 | 15.4% | 53.8% | **0.00** |
| hvg_method = cell_ranger | 17 | 0.865 | 7.6% | 77.2% | 0.06 |
| n_hvg = 4000 | 17 | 0.833 | 10.5% | 72.3% | 0.06 |
| n_neighbours = 10 | 20 | 0.848 | 12.2% | 71.1% | 0.19 |
| n_neighbours = 30 | 15 | 0.883 | 9.5% | 50.0% | 0.07 |
| resolution = 0.5 | 11 | 0.738 | 25.1% | 53.8% | **0.00** |
| filter_order = normalise_then_filter | 18 | 1.000 | 0.0% | 100.0% | 1.00 |

**Adjusted Rand index says the analysis is stable. The spatial conclusion says it
is not.** Median ARI across single-choice deviations is 0.848 — the number that
would get reported, and it reads as reassuring. Meanwhile the set of cluster pairs
called significantly co-located by neighbourhood enrichment shares **0 to 19%** of
its members with the reference. Three of the seven deviations share *none*.

Rare populations are where the churn concentrates: the least stable retains as
little as **50%** of its cells while ARI stays above 0.88. Global agreement is
dominated by abundant cell types, and the interesting population is rarely one of
them.

One honest negative: **filter order made no difference at all** (ARI 1.000, zero
churn). At these thresholds no spot is filtered differently, so the order cannot
matter. Worth stating rather than quietly dropping the axis.

### Prior work

**The pattern measured here — a global metric staying reassuring while the downstream
conclusion moves — is published for bulk transcriptomics, and this is its spatial analogue
rather than a new phenomenon.**

- Paton et al., *Nucleic Acids Research* 2023 (FLOP) — end-to-end analysis of how pipeline
  choices propagate into functional enrichment. They report "effects not noticeable at the
  gene-level" that appear in gene-set space, and find filtering has the largest impact on
  agreement between pipelines. That is the same shape as ARI 0.848 alongside 0–19%
  conclusion overlap, one level up the stack.
- Chen et al., *iMeta* 2025 — 14 spatial clustering methods over ~600 datasets, explicitly
  investigating how preprocessing pipelines influence clustering outcomes.
- Heumos et al., *Nature Reviews Genetics* 2023 — the best-practice synthesis these
  defensible-choice ranges are drawn from.

What is specific here is the readout: neighbourhood enrichment as the *conclusion* metric
rather than cluster agreement, and rare-population retention tracked separately, which is
where the churn concentrates. One dataset, single-choice deviations only, so it measures the
width of the tree near the reference path and not the whole space.

### More

- [Analysis](ANALYSIS.md) — what was done and why, including two metrics that are easy to get wrong
- [Results](results/RESULTS.md) — full results
- [Design](docs/DESIGN.md) — what is varied and what is held fixed, the layout, and the traps this avoids
