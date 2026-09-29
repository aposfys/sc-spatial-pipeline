# sc-spatial-pipeline design notes

A spatial transcriptomics analysis runs through several stages (QC and normalisation,
feature selection, graph construction, clustering, spatial statistics, interpretation), and
each stage has several defensible choices. The published figure shows one path through that
tree. This project measures how far the answer moves when one of those choices changes.

| | |
| --- | --- |
| **Stack** | scverse (`scanpy`, `squidpy`, `AnnData`), Leiden via `igraph` |
| **Data** | squidpy Visium H&E mouse brain, 2,688 spots, raw counts from `.raw` |
| **Varied** | normalisation target · HVG method · HVG count · neighbourhood size · Leiden resolution · filter order |
| **Fixed** | dataset, seed, QC thresholds (`min_genes=200`, `min_cells=3`), 50 PCs, enrichment threshold z >= 2 |
| **Readout** | ARI, matched per-spot churn, one-to-one rare-cluster retention, and matched overlap of the neighbourhood-enrichment conclusion |

There is no segmentation step and no annotation reference. Visium spots are a fixed assay
grid, and clusters are compared with clusters, not with transferred cell-type labels.

## Traps this pipeline is built to avoid

- **Normalising data that is already normalised.** squidpy's copy of this dataset has a
  log-normalised `X` and raw counts in `.raw`. Every configuration starts from `.raw`, and
  non-integer input is refused.
- **Comparing cluster numbers across runs.** Leiden numbering is arbitrary. Churn, rare
  retention and conclusion overlap all resolve it first by optimal one-to-one assignment.
  Unmatched raw churn and raw conclusion overlap are kept in `findings.json` so the size of
  the artefact stays visible.
- **Comparing runs by position.** Runs are aligned by spot barcode, so configurations that
  keep different spots are never scored index by index.
- **ARI near 1.0 can hide a lost rare population.** ARI is dominated by abundant clusters.
  Rare retention counts a rare cluster's spots that land in its own matched counterpart, so
  a cluster absorbed into another scores 0.
- **Self-pairs in neighbourhood enrichment.** Every cluster is enriched with itself, so
  self-pairs are excluded from the conclusion.
- **Filter order is a no-op here.** Both orders are run because the order is rarely
  reported, but `min_genes` and `min_cells` count non-zero entries, which normalisation
  cannot change. The result (no spot changes cluster) is expected, not evidence that order
  never matters. A value-based filter such as a total-count or mitochondrial cut would be
  needed to test it.

## Prior work

The pattern measured here, a global agreement score staying moderate while a downstream
conclusion moves more, is published for bulk transcriptomics. This is a spatial analogue,
not a new phenomenon.

- Paton V et al. Assessing the impact of transcriptomics data analysis pipelines on
  downstream functional enrichment results. *Nucleic Acids Res* 2024, 52(14), 8100-8111.
  [doi:10.1093/nar/gkae552](https://doi.org/10.1093/nar/gkae552). FLOP runs end-to-end
  pipelines and reports effects in gene-set space that are not visible at the gene level,
  with filtering having the largest impact on agreement between pipelines.
- Chen R et al. A comprehensive benchmarking for spatially resolved transcriptomics
  clustering methods across variable technologies, organs, and replicates. *iMeta* 2025,
  4(6), e70084. [doi:10.1002/imt2.70084](https://doi.org/10.1002/imt2.70084). 14 spatial
  clustering methods over about 600 datasets, including how preprocessing influences
  clustering.
- Heumos L et al. Best practices for single-cell analysis across modalities. *Nat Rev
  Genet* 2023, 24(8), 550-572.
  [doi:10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w). The
  best-practice synthesis the defensible-choice ranges are drawn from.

What is specific here is the readout. The conclusion metric is neighbourhood enrichment on
matched clusters rather than cluster agreement, and rare-cluster retention is tracked
separately.

## Layout

```
src/scspatial/
  configs.py      the grid of defensible preprocessing choices
  pipeline.py     raw counts in, one configuration end to end via scverse
  sensitivity.py  ARI, cluster matching, matched churn, rare retention, conclusion overlap
  report.py       the stability table, findings.json and runs.json
  cli.py          fetch / grid / sensitivity
results/
  RESULTS.md      generated table
  findings.json   every comparison, plus package versions
  runs.json       per-configuration spot labels and enriched pairs
```

25 tests, none of which need the analysis stack installed. `--mode full` runs all 144
combinations instead of the reference plus 8 single-axis deviations.
