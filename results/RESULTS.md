# Results

visium_hne, 2,688 spots. 9 configurations, each differing from the reference in one choice.

The reference is a standard scanpy path from the raw counts (normalise to 10,000 counts per spot, log1p, 2,000 Seurat HVGs, 15 neighbours, Leiden at resolution 1.0). It is not a claim about what is correct.

| Changed choice | Clusters | ARI | Cell churn | Worst rare retention | Conclusion overlap |
| --- | ---: | ---: | ---: | ---: | ---: |
| (reference) | 21 | 1.000 | 0.0% | 100.0% | 1.00 |
| normalisation = cpm_log1p | 18 | 0.730 | 21.9% | 0.0% | 0.24 |
| normalisation = median_log1p | 19 | 0.848 | 11.5% | 0.0% | 0.65 |
| hvg_method = cell_ranger | 21 | 0.827 | 12.1% | 0.0% | 0.60 |
| n_hvg = 4000 | 21 | 0.817 | 13.1% | 0.0% | 0.52 |
| n_neighbours = 10 | 22 | 0.830 | 13.3% | 0.0% | 0.84 |
| n_neighbours = 30 | 17 | 0.913 | 7.7% | 0.0% | 0.53 |
| resolution = 0.5 | 12 | 0.805 | 19.5% | 0.0% | 0.19 |
| filter_order = normalise_then_filter | 21 | 1.000 | 0.0% | 100.0% | 1.00 |

## What moved

- **Lowest global agreement:** `normalisation = cpm_log1p` at ARI 0.730.
- **Most cells relabelled:** `normalisation = cpm_log1p` moved 21.9% of spots to a different cluster.
- **Worst rare-population retention:** 7 of 8 deviations leave at least one rare reference cluster with 0.0% of its spots in its matched counterpart.
- **Least stable spatial conclusion:** `resolution = 0.5` shares only 19% of its neighbourhood-enrichment pairs with the reference.

Across the 7 deviations that changed the clustering the median ARI is **0.827** and the median conclusion overlap is **0.53**. `filter_order = normalise_then_filter` is left out because no spot changed cluster.

Compared on raw cluster numbers, without matching, the same conclusion overlap would read 0.03 to 0.09.

## Reading these columns

- **ARI** is the adjusted Rand index against the reference labels. It is dominated by the abundant clusters.
- **Cell churn** is the fraction of spots whose cluster changes once the arbitrary cluster numbering is resolved by optimal one-to-one assignment. Raw churn, which skips that step, runs far higher and is mostly an artefact of numbering.
- **Worst rare retention** is, for the rare reference clusters (under 5% of spots), the lowest fraction of a cluster's spots that land in its matched counterpart. A rare cluster absorbed into an abundant one scores 0.
- **Conclusion overlap** is the Jaccard of the cluster pairs called significantly co-located by neighbourhood enrichment (z >= 2, self-pairs excluded), after each run's clusters are matched onto the reference's by the same assignment. A pair involving a cluster with no counterpart cannot match.

Run with python 3.13.9, scanpy 1.12.4, squidpy 1.8.3, anndata 0.13.4, igraph 1.0.0, numpy 2.4.6, scipy 1.18.1.
