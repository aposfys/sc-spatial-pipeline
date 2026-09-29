"""Command line entry point: ``python -m scspatial.cli``."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

from scspatial import __version__

DATASETS = ("visium_hne",)
GRID_MODES = ("one_at_a_time", "full")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scspatial",
        description="How much of a spatial result is the pipeline rather than the tissue?",
    )
    parser.add_argument("--version", action="version", version=f"scspatial {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser("fetch", help="download the dataset")
    fetch.add_argument("--dataset", choices=DATASETS, default="visium_hne")

    grid = sub.add_parser("grid", help="run every configuration, storing labels")
    grid.add_argument("--dataset", choices=DATASETS, default="visium_hne")
    grid.add_argument("--mode", choices=GRID_MODES, default="one_at_a_time")
    grid.add_argument("--seed", type=int, default=0)

    sensitivity = sub.add_parser("sensitivity", help="stability across configurations")
    sensitivity.add_argument("--dataset", choices=DATASETS, default="visium_hne")

    return parser


def package_versions() -> dict[str, str]:
    """Versions of the packages that decide the clustering, for the findings file."""
    import platform
    from importlib.metadata import PackageNotFoundError, version

    versions = {"python": platform.python_version()}
    for name in ("scanpy", "squidpy", "anndata", "igraph", "numpy", "scipy"):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = "not installed"
    return versions


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command in ("fetch", "grid"):
        from scspatial.pipeline import require_stack

        require_stack()

    if args.command == "fetch":
        from scspatial.pipeline import fetch

        adata = fetch(args.dataset, cache_dir=args.data_dir)
        print(f"{args.dataset}: {adata.shape[0]:,} spots x {adata.shape[1]:,} genes")
        return 0

    if args.command == "grid":
        from scspatial.configs import iter_grid
        from scspatial.pipeline import fetch, raw_counts, run_one

        adata = raw_counts(fetch(args.dataset, cache_dir=args.data_dir))
        configs = list(iter_grid(args.mode, seed=args.seed))
        print(f"{len(configs)} configurations on {adata.shape[0]:,} spots")

        results = []
        for i, config in enumerate(configs, start=1):
            result = run_one(adata, config)
            results.append(result)
            print(
                f"  [{i}/{len(configs)}] {config.key}  "
                f"{result.n_clusters} clusters, {len(result.enriched_pairs)} pairs, "
                f"{result.seconds:.0f}s",
                flush=True,
            )

        args.results_dir.mkdir(parents=True, exist_ok=True)
        # Labels are pickled rather than written as JSON: they are the intermediate the
        # sensitivity step consumes, not a deliverable.
        with (args.results_dir / "grid.pkl").open("wb") as handle:
            pickle.dump(
                {
                    "results": results,
                    "configs": [c.as_dict() for c in configs],
                    "dataset": args.dataset,
                    "versions": package_versions(),
                },
                handle,
            )
        print(f"wrote {args.results_dir / 'grid.pkl'}")
        return 0

    if args.command == "sensitivity":
        from scspatial.report import build_findings, build_runs, write

        grid_path = args.results_dir / "grid.pkl"
        if not grid_path.exists():
            raise SystemExit(f"{grid_path} not found; run `scspatial grid` first")
        with grid_path.open("rb") as handle:
            stored = pickle.load(handle)

        findings = build_findings(
            stored["results"],
            stored["configs"],
            stored["dataset"],
            versions=stored.get("versions"),
        )
        (args.results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
        # Per-configuration labels and enriched pairs, so every column can be recomputed
        # from committed files without rerunning the grid.
        (args.results_dir / "runs.json").write_text(json.dumps(build_runs(stored["results"])))
        out = write(args.results_dir / "findings.json", args.results_dir / "RESULTS.md")
        print(f"wrote {out}")
        return 0

    raise SystemExit(f"unknown command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
