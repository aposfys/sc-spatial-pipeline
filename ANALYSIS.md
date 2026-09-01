# Analysis

What was built, why it was built that way, and two metric corrections that would each have
flattered the result.

## The question

A spatial transcriptomics analysis has several defensible choices at every stage, and the
published figure shows one path through that tree. This measures the width of the tree:
same dataset, same question, one choice changed at a time.

## Design decisions, and the reasoning

**The reference is a tutorial configuration, not a claim about what is correct.** Everything
is measured against it because that is what most published pipelines are, not because it is
right.

**One axis at a time, not the full grid, for the headline.** The full 144-combination grid
tells you the total spread. Nine single-axis deviations tell you *where the spread comes
from*, which is the more actionable question. Both are available; `--mode full` runs the
grid.

**Inclusion criterion for a choice: "would pass review unremarked".** Not merely plausible.
`normalisation = none` is in the grid because analysts do skip normalisation on
size-factor-corrected data, not because it is advisable.

**Segmentation is absent, and that is the data, not an oversight.** Spot-based Visium has no
segmentation step — the spots are a fixed assay grid, and no analyst choice changes which
transcripts land in which spot. The design notes are right that segmentation is the stage
whose errors most often become discoveries; it is simply not a choice this modality offers.
On imaging-based data it would be the first axis.

**Each configuration works on a copy.** An in-place scanpy operation leaking into the next
configuration is the single easiest way to make a sensitivity analysis silently measure
nothing.

## Two corrections that mattered

Both were caught after the first run produced numbers, and both had been making the analysis
look worse or better than it is.

**Runs are aligned by cell barcode, not by position.** Configurations that filter
differently keep different cells. Comparing two label vectors by index scores unrelated
cells against each other and returns a number that looks like instability but is a
bookkeeping error.

**Churn is computed after resolving cluster renumbering.** Leiden numbers its clusters
arbitrarily, so two runs that partition the cells *identically* can disagree on every raw
label. Raw churn on these runs reads 49–89%; matched churn, after optimal assignment,
reads 7.6–25.1%. The first number is almost entirely an artefact of numbering. A test
asserts that a pure relabelling scores exactly zero.

Optimal assignment rather than greedy matching: greedy can map two of B's clusters onto one
of A's and double-count the agreement.

## What was measured

Visium H&E, 2,688 spots, nine configurations.

| Changed choice | ARI | Cell churn | Worst rare retention | Conclusion overlap |
| --- | ---: | ---: | ---: | ---: |
| normalisation = cpm_log1p | 0.778 | 14.0% | 53.8% | **0.00** |
| normalisation = none | 0.810 | 15.4% | 53.8% | **0.00** |
| hvg_method = cell_ranger | 0.865 | 7.6% | 77.2% | 0.06 |
| n_hvg = 4000 | 0.833 | 10.5% | 72.3% | 0.06 |
| n_neighbours = 10 | 0.848 | 12.2% | 71.1% | 0.19 |
| n_neighbours = 30 | 0.883 | 9.5% | 50.0% | 0.07 |
| resolution = 0.5 | 0.738 | 25.1% | 53.8% | **0.00** |
| filter_order = normalise_then_filter | 1.000 | 0.0% | 100.0% | 1.00 |

**ARI says the analysis is stable. The spatial conclusion says it is not.** Median ARI
across deviations is 0.848 — the number that would get reported, and it reads as reassuring.
The set of cluster pairs called significantly co-located shares **0–19%** of its members
with the reference, and three of seven deviations share none.

Rare populations are where the churn concentrates: the least stable retains **50%** of its
cells while ARI stays above 0.88 — exactly the failure the design predicted, since global
agreement is dominated by abundant types and the interesting population rarely is one.

**One honest negative:** filter order changed nothing (ARI 1.000, zero churn). At these
thresholds no spot is filtered differently, so the order cannot matter. Reported rather than
quietly dropped.

## Caveats on the conclusion metric

Conclusion overlap is the Jaccard of enriched cluster-pair sets. Cluster identities are not
comparable across runs, so this is an **upper bound** on agreement — it can only overstate
stability. The instability is therefore at least as large as reported.

Self-pairs are excluded: every cluster is spatially enriched with itself, and including them
would put a guaranteed set into every configuration.

## What would change the conclusion

An imaging-based dataset, where segmentation enters the grid. The prediction from the design
notes is that it dominates everything here; that is untested.
