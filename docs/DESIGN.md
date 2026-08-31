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
