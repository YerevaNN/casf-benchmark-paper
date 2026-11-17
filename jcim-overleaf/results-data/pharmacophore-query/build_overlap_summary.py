"""Build the descriptive ATC48 exact-training-identity endpoint summary."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
METHODS = {
    "qwen_1p7b_v2_coordtoken_step41054_atc48": "Qwen 1.7B FSQ",
    "loqi_atc48": "LoQI",
    "qwen_loqi_equal_mixture": "Qwen/LoQI 50:50",
}
REPLICATES = 2000
BASE_SEED = 20261116


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def bootstrap_interval(values: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    sampled = rng.choice(values, size=(REPLICATES, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(sampled, [0.025, 0.975]).tolist())


def main() -> None:
    curves = pd.read_csv(HERE / "atc48-per-molecule-budget-curves.csv")
    overlap = pd.read_csv(HERE / "atc48-exact-identity-overlap.csv")
    overlap["exact_identity_overlap"] = overlap["exact_identity_overlap"].map(parse_bool)
    endpoint = curves[
        (curves["analysis"] == "primary")
        & np.isclose(curves["radius"], 0.5)
        & (curves["budget"] == 1000)
        & curves["method"].isin(METHODS)
    ].copy()
    endpoint = endpoint.merge(
        overlap[["expected_smiles", "exact_identity_overlap"]],
        left_on="smiles",
        right_on="expected_smiles",
        how="left",
        validate="many_to_one",
    )
    if endpoint["exact_identity_overlap"].isna().any():
        raise ValueError("An eligible molecule is missing from the exact-identity audit")

    rows = []
    for (method, exact_overlap), group in endpoint.groupby(
        ["method", "exact_identity_overlap"], sort=True
    ):
        values = group["molecule_coverage"].to_numpy(float)
        training_group = "Exact identity" if exact_overlap else "No exact identity"
        seed = BASE_SEED + sum(map(ord, f"training_overlap|{training_group}|{method}"))
        low, high = bootstrap_interval(values, seed)
        rows.append(
            {
                "method": method,
                "method_label": METHODS[method],
                "training_group": training_group,
                "coverage": float(values.mean()),
                "ci_low": low,
                "ci_high": high,
                "molecule_count": len(values),
                "bootstrap_replicates": REPLICATES,
                "bootstrap_seed": seed,
            }
        )
    summary = pd.DataFrame(rows)
    if set(summary["molecule_count"]) != {17, 27}:
        raise ValueError("Eligible exact-identity strata do not reproduce 17/27 molecules")
    summary.to_csv(HERE / "atc48-training-overlap-endpoints.csv", index=False)


if __name__ == "__main__":
    main()
