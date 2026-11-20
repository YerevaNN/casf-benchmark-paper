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
FAMILIES = [
    "individual_ligand",
    "individual_interaction",
    "shared_ligand",
    "shared_protein",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q50_difference(row: pd.Series) -> tuple[str, float, float, float]:
    if row["left"] == "qwen50_loqi50":
        return row["right"], row["difference"], row["ci_low"], row["ci_high"]
    if row["right"] == "qwen50_loqi50":
        return row["left"], -row["difference"], -row["ci_high"], -row["ci_low"]
    raise AssertionError("Archived contrast does not include 50:50")


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
    diversity = pd.read_csv(HERE / "diversity-summary.csv")
    assert len(coverage) == 2352 and len(paired) == 12 and len(diversity) == 56

    reference = coverage[
        (coverage["condition"] == "4/4") & np.isclose(coverage["radius_angstrom"], 0.5)
    ]
    k20 = reference[reference["budget"] == 20].pivot(
        index="recipe", columns="family", values="coverage"
    )
    for family in FAMILIES:
        assert k20[family].idxmax() == "qwen50_loqi50"
    q50_minus_qwen = 100 * (k20.loc["qwen50_loqi50"] - k20.loc["qwen"])
    np.testing.assert_allclose(
        [q50_minus_qwen.min(), q50_minus_qwen.max()],
        [2.43177732531387, 3.19396186182783],
        atol=1e-10,
    )

    contrast = {}
    for _, row in paired.iterrows():
        comparator, difference, low, high = q50_difference(row)
        contrast[(row["family"], comparator)] = (difference, low, high)
    for family in FAMILIES:
        assert contrast[(family, "loqi")][1] > 0
        assert contrast[(family, "qwen75_loqi25")][1] < 0 < contrast[(family, "qwen75_loqi25")][2]
    assert contrast[("individual_ligand", "qwen")][1] > 0
    assert contrast[("individual_interaction", "qwen")][1] > 0
    assert contrast[("shared_ligand", "qwen")][1] < 0 < contrast[("shared_ligand", "qwen")][2]
    assert contrast[("shared_protein", "qwen")][1] < 0 < contrast[("shared_protein", "qwen")][2]

    k1000 = reference[reference["budget"] == 1000].pivot(
        index="recipe", columns="family", values="coverage"
    )
    q50 = k1000.loc["qwen50_loqi50"]
    np.testing.assert_allclose([q50.min(), q50.max()], [0.983899, 0.995910], atol=5e-7)
    assert 100 * (k1000.max(axis=0) - q50).max() < 0.29

    endpoint = diversity[
        (diversity["budget"] == 1000) & np.isclose(diversity["cluster_radius"], 0.5)
    ].set_index("recipe")
    np.testing.assert_allclose(
        endpoint.loc[
            ["qwen", "qwen75_loqi25", "qwen50_loqi50", "qwen25_loqi75", "loqi"],
            "effective_clusters",
        ],
        [233.185991, 185.090962, 131.757643, 92.973430, 66.719306],
        atol=5e-7,
    )

    main_table = (ROOT / "tables" / "pharmacophore-mixture.tex").read_text()
    results = (ROOT / "results.tex").read_text()
    conclusions = (ROOT / "conclusions.tex").read_text()
    for token in ["82.26", "80.26", "87.54", "80.84", "0.28", "131.8"]:
        assert token in main_table, token
    assert "2,352" in results and "57,459" in results
    assert "recommend generating 50\\%" in conclusions
    print(
        "Verified 2,352 primary coverage cells, 12 paired contrasts, 56 diversity rows, "
        "all archived/rendered hashes, and the bounded 50:50 manuscript claims."
    )


if __name__ == "__main__":
    main()
