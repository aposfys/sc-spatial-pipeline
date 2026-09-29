# sc-spatial-pipeline
How much of a spatial transcriptomics result is the pipeline rather than the tissue?

[![CI](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

One dataset (squidpy Visium H&E mouse brain, 2,688 spots, raw counts), a standard scanpy
and squidpy path as the reference, and eight single-choice deviations from it. Each run is
compared with the reference on cluster agreement, rare-cluster retention and the
neighbourhood-enrichment conclusion, with cluster numbering resolved before every
comparison.

```
make install    # package, scanpy, squidpy, igraph and dev tools
make data       # squidpy Visium H&E, ~330 MB, cached in data/
make analysis   # 9 configurations, a few seconds each -> results/
make test       # 25 tests, no analysis stack needed
```

### Result

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

Across the seven deviations that changed the clustering, median ARI is **0.827** and the
median conclusion overlap is **0.53**. A single choice keeps between 19% and 84% of the
cluster pairs the reference calls significantly co-located, and each of the seven leaves at
least one rare cluster (under 5% of spots) with none of its spots in its matched
counterpart. Filter order changes nothing, because the count-based filters are blind to
normalisation.

An earlier version of this repo reported 0 to 19% conclusion overlap. That score compared
raw cluster numbers between runs. On the current runs the same unmatched comparison reads
0.03 to 0.09, so most of the old gap was numbering, not biology.

### More

- [Analysis](ANALYSIS.md) covers what was done, the metric corrections and the caveats
- [Results](results/RESULTS.md) is the generated table, with `findings.json` and per-run labels in `runs.json`
- [Design](docs/DESIGN.md) covers what is varied and held fixed, the layout and prior work
