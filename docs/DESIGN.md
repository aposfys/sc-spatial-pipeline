# sc-spatial-pipeline — design notes

A spatial transcriptomics analysis is a six-stage pipeline — raw processing and coordinate
alignment, QC and normalisation, segmentation, cell-type annotation, spatial statistics,
interpretation — and each stage has several defensible choices. The published figure shows
one path through that tree.

| | |
| --- | --- |
| **Stack** | scverse: `scanpy`, `squidpy`, `SpatialData`, `AnnData` |
| **Varied** | segmentation method · normalisation · HVG selection · neighbourhood size |
| **Fixed** | dataset, random seeds, annotation reference |
| **Readout** | Adjusted Rand Index across configurations, per-cell label churn, and whether the neighbourhood-enrichment conclusion flips |

## Why this framing, and not another benchmark

A 2026 harmonised benchmark of single-cell and spatial foundation models found that no model
dominates and that **rankings shift with preprocessing, tokenisation and metric choice**.
That result is usually read as a statement about models. It is equally a statement about
pipelines — and the pipeline half is what a working analyst actually controls.

So this is a sensitivity analysis, not a leaderboard.

## Traps this pipeline is built to avoid

- **Segmentation errors propagate as biology.** A cell wrongly merged with its neighbour
  becomes a chimeric expression profile, and chimeric profiles cluster together into what
  looks like a novel intermediate cell state. Segmentation is varied first because it is
  the stage whose errors are most often reported as discoveries.
- **Normalising before filtering changes the filter.** Order of operations inside stage 2 is
  rarely stated in methods sections and is not neutral; both orders are run.
- **The annotation reference is a prior, not ground truth.** Cell-type labels transferred
  from a reference atlas inherit that atlas's granularity and its biases. Held fixed here on
  purpose, so its contribution is not confounded with the choices under test.
- **ARI near 1.0 can still hide a flipped conclusion.** Global agreement is dominated by
  abundant cell types; a rare population can be entirely reassigned while ARI barely moves.
  Rare-population stability is reported separately, and it is the number that matters.

## Two metrics that are easy to get wrong

Both mistakes flatter the analysis:

- **Cell churn is computed after resolving cluster renumbering** by optimal
  assignment. Leiden numbers its clusters arbitrarily, so two runs that partition
  the cells identically can disagree on every raw label — raw churn on these same
  runs reads 49–89% and is almost entirely an artefact. A test asserts that a pure
  relabelling scores zero.
- **Runs are aligned by cell barcode, not by position.** Configurations that
  filter differently keep different cells, and comparing two label vectors by
  index would score unrelated cells against each other.

**Segmentation is not in the grid**, because spot-based Visium has no
segmentation step — the spots are a fixed assay grid. On imaging-based data it
would be the first axis, since it is the stage whose errors most often get
reported as discoveries.

## Layout

```
src/scspatial/
  configs.py      the grid of defensible preprocessing choices
  pipeline.py     one configuration, end to end, via scverse
  sensitivity.py  ARI, matched churn, rare-population stability, barcode alignment
  report.py       the stability table
  cli.py          fetch / grid / sensitivity
```

13 tests, none of which need the analysis stack installed. The same steps are
available as `scspatial fetch` / `grid` / `sensitivity`; `--mode full` runs all
144 combinations instead of the 9 single-axis deviations.
