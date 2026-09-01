# Results

visium_hne, 2,688 spots. 9 configurations, each differing from the reference in one choice.

The reference is the configuration a tutorial produces. It is not a claim about what is correct -- it is what most published pipelines are.

| Changed choice | Clusters | ARI | Cell churn | Worst rare retention | Conclusion overlap |
| --- | ---: | ---: | ---: | ---: | ---: |
| (reference) | 18 | 1.000 | 0.0% | 100.0% | 1.00 |
| normalisation = cpm_log1p | 18 | 0.778 | 14.0% | 53.8% | 0.00 |
| normalisation = none | 19 | 0.810 | 15.4% | 53.8% | 0.00 |
| hvg_method = cell_ranger | 17 | 0.865 | 7.6% | 77.2% | 0.06 |
| n_hvg = 4000 | 17 | 0.833 | 10.5% | 72.3% | 0.06 |
| n_neighbours = 10 | 20 | 0.848 | 12.2% | 71.1% | 0.19 |
| n_neighbours = 30 | 15 | 0.883 | 9.5% | 50.0% | 0.07 |
| resolution = 0.5 | 11 | 0.738 | 25.1% | 53.8% | 0.00 |
| filter_order = normalise_then_filter | 18 | 1.000 | 0.0% | 100.0% | 1.00 |

## What moved

- **Lowest global agreement:** `resolution = 0.5` at ARI 0.738.
- **Most cells relabelled:** `resolution = 0.5` moved 25.1% of cells to a different cluster.
- **Worst rare-population retention:** `n_neighbours = 30` — the least stable rare population kept only 50.0% of its cells together.
- **Least stable spatial conclusion:** `normalisation = cpm_log1p` shares only 0% of its neighbourhood-enrichment pairs with the reference.

Median ARI across single-choice deviations is **0.848**. Global agreement that high is the number usually reported; the rare-population and conclusion columns are where the instability actually lives.

## Reading these columns

- **ARI** is dominated by abundant populations. It is the reassuring number.
- **Cell churn** answers the question an analyst asks -- whether *this cell's* call would change -- and is computed after resolving arbitrary cluster renumbering by optimal assignment. Raw churn, which does not, runs far higher and is almost entirely an artefact of numbering.
- **Worst rare retention** is the number ARI hides. A rare population can be entirely reassigned while ARI barely moves.
- **Conclusion overlap** is the Jaccard of the significantly co-located cluster pairs. Cluster identities are not comparable across runs, so this is an upper bound on agreement, not an exact match -- it can only overstate stability.
