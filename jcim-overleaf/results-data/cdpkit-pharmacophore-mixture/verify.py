#!/usr/bin/env python3
"""Validate the archived allocation evidence and its rendered paper assets."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIGAND_FAMILIES = ["individual_ligand", "shared_ligand"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paired_difference(
    frame: pd.DataFrame,
    family: str,
    budget: int,
    left: str,
    right: str,
) -> tuple[float, float, float]:
    rows = frame[
        (frame["family"] == family)
        & (frame["budget"] == budget)
        & (
            ((frame["left"] == left) & (frame["right"] == right))
            | ((frame["left"] == right) & (frame["right"] == left))
        )
    ]
    assert len(rows) == 1, (family, budget, left, right)
    row = rows.iloc[0]
    if row["left"] == left:
        return row["difference"], row["ci_low"], row["ci_high"]
    return -row["difference"], -row["ci_high"], -row["ci_low"]


def main() -> None:
    provenance = json.loads((HERE / "provenance.json").read_text())
    for name, expected in provenance["archived_files"].items():
        assert digest(HERE / name) == expected, name
    for record in provenance["source_files"].values():
        path = Path(record["path"])
        if path.is_file():
            assert digest(path) == record["sha256"], path

    rendered = json.loads((HERE / "rendered-files.json").read_text())
    for name, expected in rendered["inputs"].items():
        assert digest(ROOT / name) == expected, name
    for name, expected in rendered["outputs"].items():
        assert digest(ROOT / name) == expected, name

    coverage = pd.read_csv(HERE / "coverage-primary.csv")
    paired = pd.read_csv(HERE / "paired-k20.csv")
    paired_ligand = pd.read_csv(HERE / "paired-ligand-only.csv")
    diversity = pd.read_csv(HERE / "diversity-summary.csv")
    assert (
        len(coverage) == 2352
        and len(paired) == 12
        and len(paired_ligand) == 168
        and len(diversity) == 56
    )

    reference = coverage[
        (coverage["condition"] == "6/6") & np.isclose(coverage["radius_angstrom"], 0.5)
    ]
    ligand_reference = reference[reference["family"].isin(LIGAND_FAMILIES)]
    ligand_means = ligand_reference.pivot_table(
        index=["family", "budget"], columns="recipe", values="coverage"
    )
    for family in LIGAND_FAMILIES:
        for budget in [20, 100, 500, 1000]:
            row = ligand_means.loc[(family, budget)]
            assert row["qwen50_loqi50"] > max(row["qwen"], row["loqi"])

    expected = {
        ("individual_ligand", 20, "loqi", "qwen"): (0.0996559, 0.0050695, 0.1950935),
        ("shared_ligand", 20, "loqi", "qwen"): (0.1459853, 0.0243874, 0.2764553),
        ("individual_ligand", 1000, "qwen", "loqi"): (0.0555556, 0.0079365, 0.1158730),
        ("shared_ligand", 1000, "qwen", "loqi"): (0.0284879, -0.0506266, 0.1353571),
    }
    for key, values in expected.items():
        np.testing.assert_allclose(
            paired_difference(paired_ligand, *key), values, atol=5e-7
        )

    for family, budget in [
        ("individual_ligand", 100),
        ("individual_ligand", 500),
    ]:
        pure = max(
            ["qwen", "loqi"],
            key=lambda recipe: ligand_means.loc[(family, budget), recipe],
        )
        assert paired_difference(
            paired_ligand, family, budget, "qwen50_loqi50", pure
        )[1] > 0

    for family in LIGAND_FAMILIES:
        rdkit = max(
            ["rdkit_raw", "rdkit_minimized"],
            key=lambda recipe: ligand_means.loc[(family, 1000), recipe],
        )
        difference = paired_difference(
            paired_ligand, family, 1000, "qwen50_loqi50", rdkit
        )
        assert difference[1] > 0

    diversity_reference = diversity[
        np.isclose(diversity["cluster_radius"], 0.5)
    ].pivot(index="budget", columns="recipe", values="effective_clusters")
    for budget in [20, 100, 500, 1000]:
        row = diversity_reference.loc[budget]
        assert (
            row["qwen"]
            > row["qwen75_loqi25"]
            > row["qwen50_loqi50"]
            > row["loqi"]
        )
    np.testing.assert_allclose(
        diversity_reference.loc[
            1000,
            ["qwen", "qwen75_loqi25", "qwen50_loqi50", "loqi"],
        ],
        [233.185991, 185.090962, 131.757643, 66.719306],
        atol=5e-7,
    )

    main_table = (ROOT / "tables" / "pharmacophore-mixture.tex").read_text()
    results = (ROOT / "results.tex").read_text()
    discussion = (ROOT / "discussion.tex").read_text()
    conclusions = (ROOT / "conclusions.tex").read_text()
    for token in [
        "+9.97",
        "+14.60",
        "+3.97",
        "+3.38",
        "+25.04",
        "+24.87",
        "+2.53",
        "+3.23",
        "+0.26",
        "-0.01",
        "233.19",
        "185.09",
        "131.76",
        "66.72",
    ]:
        assert token in main_table, token
    assert "2,352" in results and "57,459" in results
    assert "allocation decision uses the two ligand-only" in results
    assert "75:25 allocation is shown under the same coverage" in results
    assert "Generation time was not used" in discussion
    assert "recommend generating 50\\%" in conclusions
    print(
        "Verified 2,352 primary coverage cells, 168 stringent ligand-only contrasts, "
        "all archived/rendered hashes, and the coverage--diversity 50:50 claims."
    )


if __name__ == "__main__":
    main()
