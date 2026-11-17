"""Validate the archived pharmacophore-query summaries and rendered figure."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUDGETS = [10, 20, 50, 100, 200, 500, 1000]
RADII = [0.5, 0.75, 1.0]
PILOT_ENDPOINTS = {
    "qwen_1p7b_v2_coordtoken_step41054": 0.7122,
    "loqi_druglike": 0.7363,
    "rdkit_random_raw_druglike": 0.5211,
    "rdkit_random_minimized_22_druglike": 0.5293,
    "qwen_loqi_equal_mixture": 0.7812,
}
EXTENSION_ENDPOINTS = {
    "qwen_1p7b_v2_coordtoken_step41054_atc48": 0.7031,
    "loqi_atc48": 0.6830,
    "rdkit_etkdg_raw_atc48": 0.4666,
    "rdkit_etkdg_mmff94s_atc48": 0.5729,
    "qwen_loqi_equal_mixture": 0.7169,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primary_endpoint(frame: pd.DataFrame, method: str) -> float:
    rows = frame[
        (frame["analysis"] == "primary")
        & np.isclose(frame["radius"], 0.5)
        & (frame["budget"] == 1000)
        & (frame["method"] == method)
    ]
    if len(rows) != 1:
        raise AssertionError(f"Expected one endpoint row for {method}; found {len(rows)}")
    return float(rows.iloc[0]["coverage"])


def validate_curves(path: Path, expected_molecules: int) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"analysis", "method", "radius", "budget", "coverage", "molecule_count"}
    if missing := required.difference(frame.columns):
        raise AssertionError(f"{path.name} lacks {sorted(missing)}")
    if set(frame["budget"]) != set(BUDGETS):
        raise AssertionError(f"{path.name} does not contain the frozen budget grid")
    primary = frame[frame["analysis"] == "primary"]
    if set(np.round(primary["radius"].unique(), 2)) != set(RADII):
        raise AssertionError(f"{path.name} does not contain the frozen radii")
    if set(primary["molecule_count"]) != {expected_molecules}:
        raise AssertionError(f"{path.name} has an unexpected molecule denominator")
    ordered = primary.sort_values(["method", "radius", "budget"])
    decreases = ordered.groupby(["method", "radius"])["coverage"].diff() < -1e-12
    if decreases.any():
        raise AssertionError(f"{path.name} has a coverage curve that decreases with budget")
    radius_order = primary.sort_values(["method", "budget", "radius"])
    radius_decreases = radius_order.groupby(["method", "budget"])["coverage"].diff() < -1e-12
    if radius_decreases.any():
        raise AssertionError(f"{path.name} has a coverage curve that decreases with radius")
    return frame


def validate_per_molecule() -> None:
    pilot = pd.read_csv(HERE / "pilot-per-molecule-budget-curves.csv")
    extension = pd.read_csv(HERE / "atc48-per-molecule-budget-curves.csv")
    for name, frame, expected in (
        ("pilot", pilot, 22),
        ("ATC48", extension, 44),
    ):
        primary = frame[frame["analysis"] == "primary"]
        molecule_column = "smiles" if "smiles" in primary.columns else "molecule"
        if primary[molecule_column].nunique() != expected:
            raise AssertionError(f"{name} per-molecule table has an unexpected denominator")
        if primary["molecule_coverage"].lt(0).any() or primary["molecule_coverage"].gt(1).any():
            raise AssertionError(f"{name} per-molecule coverage is outside [0, 1]")


def validate_manifests() -> None:
    pilot_match = json.loads((HERE / "pilot-match-manifest.json").read_text())
    extension_match = json.loads((HERE / "atc48-match-manifest.json").read_text())
    for name, manifest in (("pilot", pilot_match), ("ATC48", extension_match)):
        if manifest.get("status") != "complete":
            raise AssertionError(f"{name} match manifest is not complete")
        if manifest.get("matching_algorithm") != "exhaustive_typed_feature_correspondence_then_proper_kabsch":
            raise AssertionError(f"{name} did not use the complete positional matcher")
        if manifest.get("radius_membership_monotonicity_violations") != 0:
            raise AssertionError(f"{name} has radius-monotonicity violations")
        if not manifest.get("all_candidates_exported") or not manifest.get(
            "all_searches_complete"
        ):
            raise AssertionError(f"{name} has incomplete candidate matching")

    expected_extension_methods = {
        "qwen_1p7b_v2_coordtoken_step41054_atc48",
        "loqi_atc48",
        "rdkit_etkdg_raw_atc48",
        "rdkit_etkdg_mmff94s_atc48",
    }
    if set(extension_match.get("methods", [])) != expected_extension_methods:
        raise AssertionError("ATC48 match manifest does not contain the four frozen methods")

    for prefix in ("pilot", "atc48"):
        summary = json.loads((HERE / f"{prefix}-summary-manifest.json").read_text())
        validation = summary.get("validation", {})
        if summary.get("status") != "complete" or not validation.get(
            "all_searches_resolved"
        ):
            raise AssertionError(f"{prefix} summary is not complete")
        if validation.get("budget_monotonicity_violations") != 0 or validation.get(
            "radius_monotonicity_violations"
        ) != 0:
            raise AssertionError(f"{prefix} summary is not monotone")

        pb = json.loads((HERE / f"{prefix}-posebusters-manifest.json").read_text())
        if (
            pb.get("status") != "complete"
            or not pb.get("failed_candidates_consume_budget")
            or pb.get("post_gate_replenishment")
        ):
            raise AssertionError(f"{prefix} PoseBusters semantics changed")

        rmsd = json.loads((HERE / f"{prefix}-rmsd-manifest.json").read_text())
        if (
            rmsd.get("status") != "complete"
            or not np.isclose(rmsd.get("pharmacophore_radius"), 0.5)
            or not np.isclose(rmsd.get("whole_molecule_rmsd_threshold"), 0.75)
        ):
            raise AssertionError(f"{prefix} RMSD agreement semantics changed")

    overlap = pd.read_csv(HERE / "atc48-exact-identity-overlap.csv")
    extension = pd.read_csv(HERE / "atc48-per-molecule-budget-curves.csv")
    eligible_smiles = set(extension.loc[extension["analysis"] == "primary", "smiles"])
    eligible = overlap[overlap["expected_smiles"].isin(eligible_smiles)]
    exact = eligible["exact_identity_overlap"].astype(str).str.lower().isin({"true", "1", "yes"})
    if (len(eligible), int(exact.sum()), int((~exact).sum())) != (44, 17, 27):
        raise AssertionError("ATC48 exact-identity strata do not reproduce 44 = 17 + 27")

    inventory = pd.read_csv(HERE / "atc48-pool-inventory.csv")
    counts = inventory.groupby("method")["stored_entries"].agg(["min", "max", "sum"])
    expected_counts = {
        "qwen_1p7b_v2_coordtoken_step41054_atc48": (992, 1000, 47943),
        "loqi_atc48": (1000, 1000, 48000),
        "rdkit_etkdg_raw_atc48": (1000, 1000, 48000),
        "rdkit_etkdg_mmff94s_atc48": (1000, 1000, 48000),
    }
    for method, expected in expected_counts.items():
        actual = tuple(map(int, counts.loc[method, ["min", "max", "sum"]]))
        if actual != expected:
            raise AssertionError(f"Unexpected ATC48 pool counts for {method}: {actual}")

    loqi_run = json.loads((HERE / "atc48-loqi-run-manifest.json").read_text())
    expected_run_values = {
        "status": "complete",
        "n_molecules": 48,
        "conformers_per_molecule": 1000,
        "total_conformers": 48000,
        "unique_coordinate_hashes": 48000,
        "duplicate_coordinate_hashes": 0,
        "merged_shard_order_validated": True,
        "shard_file_count": 232,
        "ignored_unsharded_part_count": 5,
    }
    for field, expected in expected_run_values.items():
        if loqi_run.get(field) != expected:
            raise AssertionError(
                f"Unexpected LoQI generation manifest value for {field}: "
                f"{loqi_run.get(field)!r}"
            )
    layout = loqi_run.get("shard_layout", {})
    if len(layout.get("4", [])) != 38 or len(layout.get("8", [])) != 10:
        raise AssertionError("LoQI shard layout does not reproduce 38x4 plus 10x8")


def validate_rendered_files() -> None:
    record = json.loads((HERE / "rendered-files.json").read_text())
    for relative, expected in {**record["inputs"], **record["outputs"]}.items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != expected:
            raise AssertionError(f"Rendered-file hash mismatch: {relative}")


def validate_provenance() -> None:
    record = json.loads((HERE / "provenance.json").read_text())
    for name, expected in record["copied_results"].items():
        path = HERE / name
        if not path.is_file() or sha256(path) != expected["sha256"]:
            raise AssertionError(f"Archived-result hash mismatch: {name}")
    for name, expected in record["derived_archive_files"].items():
        path = HERE / name
        if not path.is_file() or sha256(path) != expected["sha256"]:
            raise AssertionError(f"Derived-archive hash mismatch: {name}")
    loqi_run = json.loads((HERE / "atc48-loqi-run-manifest.json").read_text())
    if (
        loqi_run["generation_output_sha256"]
        != record["external_artifacts"]["loqi_atc48"]["sha256"]
    ):
        raise AssertionError("LoQI run manifest and external-pool provenance disagree")


def main() -> None:
    pilot = validate_curves(HERE / "pilot-aggregate-budget-curves.csv", 22)
    extension = validate_curves(HERE / "atc48-aggregate-budget-curves.csv", 44)
    for method, expected in PILOT_ENDPOINTS.items():
        actual = primary_endpoint(pilot, method)
        if not np.isclose(actual, expected, atol=5e-5):
            raise AssertionError(f"Pilot endpoint changed for {method}: {actual:.8f}")
    for method, expected in EXTENSION_ENDPOINTS.items():
        actual = primary_endpoint(extension, method)
        if not np.isclose(actual, expected, atol=5e-5):
            raise AssertionError(f"ATC48 endpoint changed for {method}: {actual:.8f}")
    validate_per_molecule()
    validate_manifests()
    validate_rendered_files()
    validate_provenance()
    print("pharmacophore-query archive verified")


if __name__ == "__main__":
    main()
