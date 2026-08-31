# sc-spatial-pipeline
How much of a spatial transcriptomics result is the pipeline rather than the tissue?

[![CI](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/sc-spatial-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** The stability machinery and its tests exist; no dataset has been processed.

A spatial analysis is a six-stage pipeline and each stage has several defensible choices. The published figure shows one path through that tree; this repo measures the width of the tree. Same dataset, same biological question, a grid of reasonable preprocessing choices — and a direct answer to how many cells change cell type, and whether the spatial conclusion survives.

This is a sensitivity analysis, not a leaderboard. The output is not "method X wins"; it is "the conclusion is stable across the choices that were never reported."

### Running it
```
make install && make data && make grid && make analysis && make test
```
`scanpy`, `squidpy` and `spatialdata` are optional extras imported lazily; the stability machinery is pure Python and is what CI exercises.

### Layout
```
src/scspatial/
  sensitivity.py  agreement between runs: ARI, label churn, rare-population stability
```
Planned: `configs.py` (the grid of defensible choices), `pipeline.py` (one configuration end to end), `report.py`.

### Design notes
[What is varied and what is held fixed, and the traps the pipeline is built to avoid](docs/DESIGN.md)
