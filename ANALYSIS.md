# Analysis

What was built, why it was built that way, and the corrections that changed the result.

## The question

A spatial transcriptomics analysis has several defensible choices at every stage, and the
published figure shows one path through that tree. This measures the width of the tree
near one path. Same dataset, same question, one choice changed at a time.

## Design decisions, and the reasoning

**The reference is a standard scanpy path, not a claim about what is correct.** Raw counts,
normalised to 10,000 counts per spot, log1p, 2,000 Seurat highly variable genes, 15
neighbours, Leiden at resolution 1.0, then squidpy neighbourhood enrichment.

**One axis at a time, not the full grid, for the headline.** The full 144-combination grid
gives the total spread. Eight single-axis deviations show where the spread comes from.
`--mode full` runs the grid.

**Inclusion criterion for a choice is "would pass review unremarked".** Every normalisation
arm is total-count scaling plus log1p, differing only in the target sum. `median_log1p` is
scanpy's default and the scanpy spatial tutorial's choice.

**Segmentation is absent, because the data has none.** Spot-based Visium has a fixed assay
grid, and no analyst choice changes which transcripts land in which spot. On imaging-based
data segmentation would be the first axis.

**Each configuration works on a copy.** An in-place scanpy operation leaking into the next
configuration would make the analysis silently measure nothing.

## Corrections

**Every configuration starts from raw counts.** squidpy ships `visium_hne_adata` already
through the scanpy tutorial. Its `X` is log-normalised to the median total count (27,237)
and the integer counts sit in `.raw`. An earlier version normalised `X` again, so every
configuration ran on a double transform. The grid now starts from `.raw`, and a check
refuses non-integer input. The old `normalisation = none` arm was in effect the shipped,
median-normalised data. On raw counts that arm crashes Seurat HVG selection, so it is now
`median_log1p`, which is what the shipped data had.

**Clusters are matched before any comparison between runs.** Leiden numbers its clusters
arbitrarily, so two runs that partition the spots identically can disagree on every raw
label. Churn, rare retention and conclusion overlap all map each run's clusters onto the
reference by optimal one-to-one assignment on the contingency table first. Optimal rather
than greedy, because greedy can map two clusters onto one and double-count agreement. Tests
assert that a pure renumbering scores zero churn, full rare retention and a conclusion
overlap of 1.

An earlier version matched clusters for churn but not for conclusion overlap, which
compared raw cluster numbers of the enriched pairs. It reported 0 to 19% overlap and called
the score an upper bound on agreement. It was not a bound. A renumbered but identical run
scored 0. On the current runs the unmatched comparison reads 0.03 to 0.09, against 0.19 to
0.84 matched, so most of the old gap was numbering.

**Rare retention is one-to-one.** It used to take each rare cluster's best partner, which
scores a rare cluster absorbed whole into an abundant one as 1.0. It now counts the spots
that land in the rare cluster's matched counterpart, so absorption scores 0.

**Runs are aligned by spot barcode, not by position**, so configurations that keep
different spots are never compared index by index.

## What was measured

Visium H&E, 2,688 spots, nine configurations. The reference finds 21 clusters, 12 of them
rare (under 5% of spots), and 16 significantly co-located cluster pairs (z >= 2, self-pairs
excluded).

| Changed choice | ARI | Cell churn | Worst rare retention | Conclusion overlap | Unmatched overlap |
| --- | ---: | ---: | ---: | ---: | ---: |
| normalisation = cpm_log1p | 0.730 | 21.9% | 0.0% | 0.24 | 0.04 |
| normalisation = median_log1p | 0.848 | 11.5% | 0.0% | 0.65 | 0.08 |
| hvg_method = cell_ranger | 0.827 | 12.1% | 0.0% | 0.60 | 0.03 |
| n_hvg = 4000 | 0.817 | 13.1% | 0.0% | 0.52 | 0.03 |
| n_neighbours = 10 | 0.830 | 13.3% | 0.0% | 0.84 | 0.09 |
| n_neighbours = 30 | 0.913 | 7.7% | 0.0% | 0.53 | 0.04 |
| resolution = 0.5 | 0.805 | 19.5% | 0.0% | 0.19 | 0.06 |
| filter_order = normalise_then_filter | 1.000 | 0.0% | 100.0% | 1.00 | 1.00 |

Across the seven deviations that changed the clustering, median ARI is 0.827 and median
conclusion overlap is 0.53. The target sum of the normalisation moves the most spots
(21.9% at CPM against CP10k) and resolution 0.5 moves the spatial conclusion the most, with
3 enriched pairs against the reference's 16. Raw churn, without matching, reads 69 to 99%
on the same runs.

Every one of the seven leaves at least one rare cluster with none of its spots in its
matched counterpart, while ARI stays between 0.73 and 0.91. The per-run labels in
`results/runs.json` show that in each case at least one rare cluster moved as a block into
a cluster matched to a different reference cluster.

Filter order changed nothing, and it cannot. `min_genes` and `min_cells` count non-zero
entries, and normalisation and log1p keep zeros at zero. The axis is kept so the full grid
matches its description, and it is left out of the medians.

## Caveats

- One dataset and single-choice deviations only. This is the spread near one path, not
  across the whole space.
- The worst-rare column reports a minimum over 12 rare clusters, so a single absorbed
  cluster sets it to 0. The per-cluster values can be recomputed from `runs.json`.
- Conclusion overlap depends on the z >= 2 threshold and on the matching. A pair that
  involves a cluster with no counterpart in the reference cannot match.
- Results are deterministic for a fixed seed and package set. Two runs of the grid gave
  byte-identical labels. The versions used are recorded at the end of
  `results/RESULTS.md`.
