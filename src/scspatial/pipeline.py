"""One configuration, end to end, through scverse.

The pipeline itself is deliberately ordinary. It is a tutorial pipeline, because a
sensitivity analysis of an exotic pipeline would measure the exotic pipeline. What is not
ordinary is that every choice it makes is named in a :class:`~scspatial.configs.Config` and
nothing is hard-coded, so two runs differ in exactly the ways the grid says they do.

`scanpy` and `squidpy` are optional extras, imported inside the functions that need them,
so the stability layer and its tests run without the analysis stack.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

from scspatial.configs import Config


@dataclass
class RunResult:
    """One configuration's labels and the spatial conclusion drawn from them."""

    key: str
    #: Cell barcodes, parallel to ``labels``. Configurations that filter differently keep
    #: different cells, so any comparison between two runs has to align on these rather
    #: than on position. Comparing two label vectors by index would silently score
    #: different cells against each other.
    cells: list[str]
    labels: list[str]
    n_cells: int
    n_clusters: int
    #: Pairs of cluster labels called significantly co-located by neighbourhood enrichment.
    enriched_pairs: list[tuple[str, str]]
    seconds: float


INSTALL_HINT = "the analysis stack is not installed. Run: pip install -e '.[sc,spatial]'"


def require_stack() -> None:
    """Exit with an install hint, not a traceback, when scanpy or squidpy is missing."""
    try:
        import scanpy  # noqa: F401
        import squidpy  # noqa: F401
    except ImportError as error:
        raise SystemExit(
            f"{error.name or 'scanpy/squidpy'} not found. {INSTALL_HINT}"
        ) from None


def fetch(dataset: str = "visium_hne", cache_dir: Path | None = None):
    """Load a public spatial dataset, downloading it into ``cache_dir`` once.

    Visium H&E from squidpy. Spot-based, which is why segmentation is not in the grid (see
    the note in :mod:`scspatial.configs`). Without ``cache_dir`` squidpy uses its own cache.
    """
    import squidpy as sq

    warnings.filterwarnings("ignore")
    if dataset != "visium_hne":
        raise ValueError(f"unknown dataset {dataset!r}; only 'visium_hne' is wired up")
    if cache_dir is None:
        return sq.datasets.visium_hne_adata()
    cache_dir.mkdir(parents=True, exist_ok=True)
    return sq.datasets.visium_hne_adata(path=cache_dir / "visium_hne_adata.h5ad")


def looks_like_counts(matrix) -> bool:
    """Whether a matrix holds raw counts: non-negative whole numbers."""
    import numpy as np

    values = matrix.data if hasattr(matrix, "data") and hasattr(matrix, "nnz") else matrix
    values = np.asarray(values)
    if values.size == 0:
        return True
    return bool((values >= 0).all() and np.array_equal(values, np.round(values)))


def raw_counts(adata):
    """The dataset as raw counts, which is where every configuration has to start.

    squidpy ships ``visium_hne_adata`` already through the scanpy tutorial, so ``.X`` is
    normalised and log-transformed and the counts sit in ``.raw``. Normalising ``.X``
    again would put every configuration on a double transform. Results from the earlier
    tutorial run (clusters, neighbours, PCA) are dropped so nothing downstream can pick
    them up by accident.
    """
    import anndata as ad

    if adata.raw is not None:
        counts = ad.AnnData(
            X=adata.raw.X.copy(),
            obs=adata.obs[[]].copy(),
            var=adata.raw.var[[]].copy(),
            obsm={"spatial": adata.obsm["spatial"].copy()},
            uns={"spatial": adata.uns["spatial"]} if "spatial" in adata.uns else {},
        )
    else:
        counts = adata.copy()
    if not looks_like_counts(counts.X):
        raise ValueError("expected raw counts in .raw or .X, found non-integer values")
    return counts


def _normalise(adata, how: str) -> None:
    import scanpy as sc

    targets = {"cpm_log1p": 1e6, "cp10k_log1p": 1e4, "median_log1p": None}
    if how not in targets:
        raise ValueError(f"unknown normalisation {how!r}; expected one of {sorted(targets)}")
    sc.pp.normalize_total(adata, target_sum=targets[how])
    sc.pp.log1p(adata)


def _filter(adata) -> None:
    import scanpy as sc

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)


def run_one(adata, config: Config) -> RunResult:
    """Run one configuration and return its labels plus its spatial conclusion.

    ``adata`` must hold raw counts (see :func:`raw_counts`). It is copied, so a
    configuration cannot contaminate the next one through an in-place scanpy operation,
    which is the easiest way to make a sensitivity analysis silently measure nothing.
    """
    import time

    import scanpy as sc
    import squidpy as sq

    warnings.filterwarnings("ignore")
    started = time.time()
    working = adata.copy()

    # Both orders are run because the order is rarely stated in methods sections. With
    # these filters it cannot matter. min_genes and min_cells count non-zero entries,
    # and normalisation and log1p leave zeros at zero.
    if config.filter_order == "filter_then_normalise":
        _filter(working)
        _normalise(working, config.normalisation)
    else:
        _normalise(working, config.normalisation)
        _filter(working)

    sc.pp.highly_variable_genes(working, flavor=config.hvg_method, n_top_genes=config.n_hvg)
    working = working[:, working.var["highly_variable"]].copy()

    sc.pp.scale(working, max_value=10)
    sc.tl.pca(working, n_comps=50, svd_solver="arpack", random_state=config.seed)
    sc.pp.neighbors(working, n_neighbors=config.n_neighbours, random_state=config.seed)
    sc.tl.leiden(
        working,
        resolution=config.resolution,
        random_state=config.seed,
        key_added="cluster",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    # The spatial conclusion under test: which cluster pairs are called co-located.
    sq.gr.spatial_neighbors(working)
    sq.gr.nhood_enrichment(
        working, cluster_key="cluster", seed=config.seed, show_progress_bar=False
    )
    enriched = _significant_pairs(working)

    labels = [str(value) for value in working.obs["cluster"]]
    return RunResult(
        key=config.key,
        cells=[str(name) for name in working.obs_names],
        labels=labels,
        n_cells=len(labels),
        n_clusters=len(set(labels)),
        enriched_pairs=enriched,
        seconds=time.time() - started,
    )


def _significant_pairs(adata, *, z_threshold: float = 2.0) -> list[tuple[str, str]]:
    """Cluster pairs whose neighbourhood-enrichment z-score clears ``z_threshold``.

    Self-pairs are excluded: every cluster is spatially enriched with itself, so including
    them would put a guaranteed set of pairs into every configuration and make the
    conclusion look far more stable than it is.
    """
    import numpy as np

    scores = adata.uns["cluster_nhood_enrichment"]["zscore"]
    categories = list(adata.obs["cluster"].cat.categories)
    pairs: list[tuple[str, str]] = []
    for i in range(len(categories)):
        for j in range(i + 1, len(categories)):
            if np.isfinite(scores[i, j]) and scores[i, j] >= z_threshold:
                pairs.append((str(categories[i]), str(categories[j])))
    return pairs
