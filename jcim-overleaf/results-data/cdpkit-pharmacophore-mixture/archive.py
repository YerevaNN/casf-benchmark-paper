#!/usr/bin/env python3
"""Archive the compact CDPKit mixture-policy evidence used by the manuscript."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = (
    HERE.parents[3]
    / "casf-benchmark"
    / "data"
    / "results"
    / "cdpkit_plinder22_exploration_v1"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def primary_arms(frame: pd.DataFrame) -> pd.DataFrame:
    return frame[
        (
            frame["family"].str.startswith("individual")
            & (frame["construction_count"] == 1)
            & (frame["support_rule"] == "single_pose")
        )
        | (
            frame["family"].str.startswith("shared")
            & (frame["construction_count"] == 2)
            & (frame["support_rule"] == "all")
        )
    ].copy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE,
        help="Completed cdpkit_plinder22_exploration_v1 output directory.",
    )
    args = parser.parse_args()
    source = args.source.resolve()
    tables = source / "tables"

    required = {
        "coverage_summary": tables / "coverage_summary.csv",
        "paired_recipe_differences": tables / "paired_recipe_differences.csv",
        "diversity_summary": tables / "diversity_summary.csv",
        "query_manifest": source / "queries" / "query_manifest.csv",
        "parent_manifest": source / "queries" / "parent_manifest.csv",
        "validity_counts": source / "validity" / "method_counts.csv",
        "protocol": source / "protocol.json",
        "run_status": source / "run_status.json",
        "query_examples": source / "figures" / "query_family_examples.pdf",
        "generator_grid": source / "figures" / "generator_separation_heatmap.pdf",
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing source artifacts: {missing}")

    status = json.loads(required["run_status"].read_text())
    if status.get("status") != "complete" or status.get("failed_screening_shards") != 0:
        raise ValueError("The source analysis is not a complete zero-failure run")

    coverage = primary_arms(pd.read_csv(required["coverage_summary"], low_memory=False))
    expected_cells = 4 * 7 * 3 * 7 * 4
    if len(coverage) != expected_cells:
        raise ValueError(f"Expected {expected_cells} primary coverage rows, found {len(coverage)}")
    coverage = coverage.sort_values(
        ["family", "condition", "radius_angstrom", "recipe", "budget"]
    )
    coverage.to_csv(HERE / "coverage-primary.csv", index=False)

    paired = primary_arms(pd.read_csv(required["paired_recipe_differences"], low_memory=False))
    paired = paired[
        (paired["condition"] == "4/4")
        & (paired["radius_angstrom"] == 0.5)
        & (paired["budget"] == 20)
        & (
            (
                (paired["left"] == "qwen50_loqi50")
                & paired["right"].isin(["qwen", "loqi", "qwen75_loqi25"])
            )
            | (
                (paired["right"] == "qwen50_loqi50")
                & paired["left"].isin(["qwen", "loqi", "qwen75_loqi25"])
            )
        )
    ].sort_values(["family", "left", "right"])
    if len(paired) != 12:
        raise ValueError(f"Expected 12 K=20 paired contrasts, found {len(paired)}")
    paired.to_csv(HERE / "paired-k20.csv", index=False)

    paired_ligand = primary_arms(
        pd.read_csv(required["paired_recipe_differences"], low_memory=False)
    )
    paired_ligand = paired_ligand[
        paired_ligand["family"].isin(["individual_ligand", "shared_ligand"])
        & (paired_ligand["condition"] == "6/6")
        & (paired_ligand["radius_angstrom"] == 0.5)
    ].sort_values(["family", "budget", "left", "right"])
    if len(paired_ligand) != 168:
        raise ValueError(
            f"Expected 168 stringent ligand-only paired contrasts, found {len(paired_ligand)}"
        )
    paired_ligand.to_csv(HERE / "paired-ligand-only.csv", index=False)

    diversity = pd.read_csv(required["diversity_summary"], low_memory=False)
    diversity = diversity.sort_values(["cluster_radius", "recipe", "budget"])
    diversity.to_csv(HERE / "diversity-summary.csv", index=False)

    queries = pd.read_csv(required["query_manifest"], low_memory=False)
    parents = pd.read_csv(required["parent_manifest"], low_memory=False)
    scope = []
    for family in sorted(queries["family"].unique()):
        family_queries = queries[queries["family"] == family]
        family_parents = parents[parents["family"] == family]
        scope.append(
            {
                "family": family,
                "queries": len(family_queries),
                "parents": len(family_parents),
                "molecules_with_queries": family_queries["smiles_id"].nunique(),
            }
        )
    pd.DataFrame(scope).to_csv(HERE / "query-scope.csv", index=False)
    shutil.copy2(required["validity_counts"], HERE / "validity-method-counts.csv")
    shutil.copy2(required["protocol"], HERE / "protocol.json")
    shutil.copy2(required["run_status"], HERE / "run-status.json")
    shutil.copy2(required["query_examples"], HERE / "query-family-examples.pdf")
    shutil.copy2(required["generator_grid"], HERE / "generator-separation-grid.pdf")

    archived = [
        "coverage-primary.csv",
        "paired-k20.csv",
        "paired-ligand-only.csv",
        "diversity-summary.csv",
        "query-scope.csv",
        "validity-method-counts.csv",
        "protocol.json",
        "run-status.json",
        "query-family-examples.pdf",
        "generator-separation-grid.pdf",
    ]
    provenance = {
        "source_directory": str(source),
        "source_files": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in required.items()
        },
        "source_run_status": status,
        "archived_files": {
            name: sha256(HERE / name)
            for name in archived
        },
        "selection": (
            "Primary individual arms use construction_count=1/single_pose; primary shared arms "
            "use construction_count=2/all. The main allocation decision uses the two ligand-only "
            "families at the prespecified stringent 6/6-feature, 0.5-angstrom condition; the two "
            "protein-informed families are secondary validation. The complete grid and unequal "
            "allocations remain archived."
        ),
    }
    (HERE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
