# sc-spatial-pipeline
How much of a spatial transcriptomics result is the pipeline rather than the tissue?

[![CI](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A spatial analysis has several defensible choices at every stage, and the published figure shows one path through that tree. This measures the width of the tree: same dataset (squidpy Visium H&E, 2,688 spots), same question, one preprocessing choice changed at a time.

### The result

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

**Adjusted Rand index says the analysis is stable. The spatial conclusion says it is not.** Median ARI across single-choice deviations is 0.848 — the number that would get reported, and it reads as reassuring. Meanwhile the set of cluster pairs called significantly co-located by neighbourhood enrichment shares **0 to 19%** of its members with the reference. Three of the seven deviations share *none*.

Rare populations are where the churn concentrates: the least stable one retains as little as **50%** of its cells while ARI stays above 0.88. That is precisely the failure the design predicted — global agreement is dominated by abundant cell types, and the interesting population is rarely one of them.

One honest negative: **filter order made no difference at all** on this dataset (ARI 1.000, zero churn). At these thresholds no spot is filtered differently, so the order cannot matter. Worth stating rather than quietly dropping the axis.

### Running it

```
make install
make data       # squidpy Visium H&E, ~330 MB, cached
make grid       # 9 configurations, ~5 s each
make analysis   # stability table -> results/RESULTS.md
make test
```

The same steps are available as `scspatial fetch` / `grid` / `sensitivity`. `--mode full` runs all 144 combinations instead of the 9 single-axis deviations.

### A note on the metrics

Two of them are easy to get wrong, and both mistakes flatter the analysis:

- **Cell churn is computed after resolving cluster renumbering** by optimal assignment. Leiden numbers its clusters arbitrarily, so two runs that partition the cells identically can disagree on every raw label — raw churn on these same runs reads 49–89% and is almost entirely an artefact. A test asserts that a pure relabelling scores zero.
- **Runs are aligned by cell barcode, not by position.** Configurations that filter differently keep different cells, and comparing two label vectors by index would score unrelated cells against each other.

**Segmentation is not in the grid**, because spot-based Visium has no segmentation step — the spots are a fixed assay grid. On imaging-based data it would be the first axis, since it is the stage whose errors most often get reported as discoveries.

### Layout

```
src/scspatial/
  configs.py      the grid of defensible preprocessing choices
  pipeline.py     one configuration, end to end, via scverse
  sensitivity.py  ARI, matched churn, rare-population stability, barcode alignment
  report.py       the stability table
  cli.py          fetch / grid / sensitivity
```

13 tests, none of which need the analysis stack installed.

### More

- [Analysis: what was done, and why it was done that way](ANALYSIS.md)
- [Full results](results/RESULTS.md)
- [What is varied and what is held fixed, and the traps the pipeline is built to avoid](docs/DESIGN.md)
