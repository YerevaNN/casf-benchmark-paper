"""PLINDER panel metrics for the dashboard.

Full-pool greedy clustering, common-hydrogen MMFF94s on the PoseBusters-passing
subset, and inclusive COV-R/COV-P at the headline cutoffs. RMSD matrices are
reused; coverage and K curves never recompute them.
"""

from __future__ import annotations

import math
import pickle
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from casf_benchmark.analysis.druglike_covmat import (
    compute_rmsd_matrices,
    coverage_fraction,
    finite_axis_min,
    load_rmsd_matrices,
    save_rmsd_matrices,
)
from casf_benchmark.analysis.metrics import greedy_cluster_metrics, mol_no_h, threshold_tag

HEADLINE_THRESHOLDS = (0.5, 0.75, 1.0)
RECOVERY_STOP = 1.0
RECOVERY_STEP = 0.05
MMFF_VARIANT = "MMFF94s"
MMFF_NONBONDED_THRESHOLD = 100.0

# Per-molecule column -> summary column. Summary values are equal-weight means.
SUMMARY_MEAN_COLUMNS = (
    ("n_generated", "mean_n_generated"),
    ("pb_pass_rate", "mean_pb_pass_rate"),
    ("greedy_clusters_0p5", "mean_clusters_0p5"),
    ("greedy_clusters_1p0", "mean_clusters_1p0"),
    ("greedy_clusters_2p0", "mean_clusters_2p0"),
    ("greedy_clusters_3p0", "mean_clusters_3p0"),
    ("clusters_per_100_1p0", "clusters_per_100_1p0"),
    ("cluster_entropy_1p0", "cluster_entropy_1p0"),
    ("largest_cluster_fraction_0p5", "largest_cluster_fraction_0p5"),
    ("largest_cluster_fraction_1p0", "largest_cluster_fraction_1p0"),
    ("effective_clusters_1p0", "effective_clusters_1p0"),
    ("singleton_fraction_1p0", "singleton_fraction_1p0"),
    ("energy_count", "energy_count"),
    ("energy_min", "energy_min"),
    ("energy_max", "energy_max"),
    ("energy_median", "energy_median"),
    ("energy_std", "energy_std"),
    ("cov_r_0p5", "cov_r_0p5_mean"),
    ("cov_p_0p5", "cov_p_0p5_mean"),
    ("cov_r_0p75", "cov_r_0p75_mean"),
    ("cov_p_0p75", "cov_p_0p75_mean"),
    ("cov_r_1p0", "cov_r_1p0_mean"),
    ("cov_p_1p0", "cov_p_1p0_mean"),
    ("mat_r", "mat_r_mean"),
    ("mat_p", "mat_p_mean"),
)

METHOD_FILES = ("summary.csv", "per_molecule.csv", "recovery.csv", "k_efficiency.csv")
DEFAULT_REFERENCE_DIR = Path(
    "/mnt/weka/mbedrosian/data/plinder/cleaned/split_by_molecule_dedup_0.05_anchor_triplet"
)

# Columns produced from the reference×generated RMSD matrix. A new reference
# set replaces these; clustering, PoseBusters, and energy stay with the pool.
_REFERENCE_COLUMNS = frozenset({"mat_r", "mat_p", "num_true_confs", "num_gen_confs"})
_DESCRIPTOR_COLUMNS = frozenset(
    {"label", "generator", "display_label", "model_size", "tokenizer", "recipe", "step"}
)
PLINDER_TABLES = (
    "plinder_summary",
    "plinder_per_molecule",
    "plinder_recovery_curve",
    "plinder_k_efficiency",
)


def recovery_thresholds(
    stop: float = RECOVERY_STOP,
    step: float = RECOVERY_STEP,
) -> tuple[float, ...]:
    """Cutoffs from 0 through `stop` at `step`, always including the headlines."""
    count = int(round(float(stop) / float(step)))
    grid = [round(index * float(step), 2) for index in range(count + 1)]
    for headline in HEADLINE_THRESHOLDS:
        if headline not in grid:
            grid.append(float(headline))
    return tuple(sorted(set(grid)))


def _finite(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array[np.isfinite(array)]


def _matching(values: np.ndarray) -> float:
    finite = _finite(values)
    if finite.size == 0:
        return float("nan")
    return float(np.mean(finite))


def headline_coverage(matrix: np.ndarray, *, inclusive: bool = True) -> dict[str, float | int]:
    """COV at 0.5, 0.75, and 1.0 Å plus threshold-free MAT for one molecule."""
    min_true = finite_axis_min(matrix, axis=1)
    min_generated = finite_axis_min(matrix, axis=0)
    finite_true = _finite(min_true)
    finite_generated = _finite(min_generated)
    row: dict[str, float | int] = {
        "mat_r": _matching(finite_true),
        "mat_p": _matching(finite_generated),
        "num_true_confs": int(matrix.shape[0]),
        "num_gen_confs": int(matrix.shape[1]),
    }
    for threshold in HEADLINE_THRESHOLDS:
        tag = threshold_tag(threshold)
        row[f"cov_r_{tag}"] = coverage_fraction(finite_true, threshold, inclusive=inclusive)
        row[f"cov_p_{tag}"] = coverage_fraction(finite_generated, threshold, inclusive=inclusive)
    return row


def recovery_curve(
    matrices: dict[str, np.ndarray],
    thresholds: tuple[float, ...] | list[float],
    *,
    inclusive: bool = True,
) -> list[dict[str, float | int]]:
    """Equal-weight mean COV-R and COV-P at each cutoff. Minima are computed once."""
    prepared: list[tuple[np.ndarray, np.ndarray]] = []
    for matrix in matrices.values():
        prepared.append((_finite(finite_axis_min(matrix, axis=1)), _finite(finite_axis_min(matrix, axis=0))))
    rows: list[dict[str, float | int]] = []
    for threshold in thresholds:
        cov_r = [coverage_fraction(true, threshold, inclusive=inclusive) for true, _generated in prepared]
        cov_p = [coverage_fraction(generated, threshold, inclusive=inclusive) for _true, generated in prepared]
        rows.append(
            {
                "threshold": float(threshold),
                "n_molecules": len(prepared),
                "cov_r_mean": float(np.nanmean(cov_r)) if cov_r else float("nan"),
                "cov_p_mean": float(np.nanmean(cov_p)) if cov_p else float("nan"),
            }
        )
    return rows


def common_hydrogen_mmff_energy(mol: Chem.Mol) -> float:
    """MMFF94s after stripping and rebuilding hydrogens, with no minimization.

    Matches the CASF energy rescore: heavy-atom coordinates stay put, explicit
    hydrogens are rebuilt with ``AddHs(addCoords=True)``, and molecules without
    MMFF parameters (bortezomib) return NaN.
    """
    try:
        heavy = Chem.RemoveHs(mol)
    except Exception:
        return math.nan
    if heavy.GetNumConformers() == 0:
        return math.nan
    coordinates = np.asarray(heavy.GetConformer().GetPositions(), dtype=float)
    if coordinates.size == 0 or not np.isfinite(coordinates).all():
        return math.nan
    try:
        prepared = Chem.AddHs(heavy, addCoords=True)
    except Exception:
        return math.nan
    try:
        if not AllChem.MMFFHasAllMoleculeParams(prepared):
            return math.nan
        properties = AllChem.MMFFGetMoleculeProperties(prepared, mmffVariant=MMFF_VARIANT)
        if properties is None:
            return math.nan
        forcefield = AllChem.MMFFGetMoleculeForceField(
            prepared,
            properties,
            confId=0,
            nonBondedThresh=MMFF_NONBONDED_THRESHOLD,
            ignoreInterfragInteractions=True,
        )
        if forcefield is None:
            return math.nan
        return float(forcefield.CalcEnergy())
    except Exception:
        return math.nan


def energy_summary(mols: list[Chem.Mol]) -> dict[str, float]:
    values = np.asarray([common_hydrogen_mmff_energy(mol) for mol in mols], dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {
            "energy_count": 0.0,
            "energy_min": math.nan,
            "energy_max": math.nan,
            "energy_median": math.nan,
            "energy_std": math.nan,
        }
    return {
        "energy_count": float(finite.size),
        "energy_min": float(np.min(finite)),
        "energy_max": float(np.max(finite)),
        "energy_median": float(np.median(finite)),
        "energy_std": float(np.std(finite)),
    }


def _empty_coverage() -> dict[str, float]:
    row = {"mat_r": math.nan, "mat_p": math.nan, "num_true_confs": 0, "num_gen_confs": 0}
    for threshold in HEADLINE_THRESHOLDS:
        tag = threshold_tag(threshold)
        row[f"cov_r_{tag}"] = math.nan
        row[f"cov_p_{tag}"] = math.nan
    return row


def reference_smiles(mol: Chem.Mol) -> str:
    """Isomeric heavy-atom SMILES used as the generation-pickle key."""
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def load_reference_sdf_dir(path: Path) -> dict[str, dict]:
    """Load one cleaned SDF per molecule into the ground-truth dict.

    Each file is a deduplicated set of experimental conformers for one CCD
    code. Records carry ``molecule_name``. The dict key is the heavy-atom
    SMILES so it joins the generation pickles. Molecules are ordered by name.
    """
    sdf_paths = sorted(Path(path).glob("*.sdf"))
    if not sdf_paths:
        raise FileNotFoundError(f"No SDF files in {path}")
    parsed: list[tuple[str, str, dict]] = []
    for sdf_path in sdf_paths:
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
        mols = [mol for mol in supplier if mol is not None and mol.GetNumConformers() > 0]
        if not mols:
            raise ValueError(f"{sdf_path} has no readable conformers")
        name = mols[0].GetProp("molecule_name") if mols[0].HasProp("molecule_name") else sdf_path.stem
        ccd = mols[0].GetProp("ccd") if mols[0].HasProp("ccd") else sdf_path.stem
        smiles = reference_smiles(mols[0])
        parsed.append((name, smiles, {"confs": mols, "names": [name], "ccd": ccd}))
    parsed.sort(key=lambda item: item[0])
    ground_truth: dict[str, dict] = {}
    for name, smiles, row in parsed:
        previous = ground_truth.get(smiles)
        if previous is not None:
            raise ValueError(f"{name} and {previous['names'][0]} share SMILES {smiles}")
        ground_truth[smiles] = row
    return ground_truth


def pool_frame_for_reference(frame: pd.DataFrame, ground_truth: dict) -> pd.DataFrame | None:
    """Keep clustering, PoseBusters, and energy for the current reference molecules.

    Returns None when any reference SMILES is missing, so the caller rescores
    the pools. Reference-derived columns are dropped and filled again from the
    new RMSD matrices.
    """
    if frame.empty or "smiles" not in frame.columns:
        return None
    if "n_generated" not in frame.columns or "greedy_clusters_1p0" not in frame.columns:
        return None
    needed = [str(smiles) for smiles in ground_truth]
    have = set(frame["smiles"].astype(str))
    if not set(needed).issubset(have):
        return None
    subset = frame[frame["smiles"].astype(str).isin(needed)].drop_duplicates("smiles").copy()
    order = {smiles: index for index, smiles in enumerate(needed)}
    subset["_order"] = subset["smiles"].astype(str).map(order)
    subset = subset.sort_values("_order").drop(columns="_order")
    drop = [
        column
        for column in subset.columns
        if column in _REFERENCE_COLUMNS or column in _DESCRIPTOR_COLUMNS or column.startswith("cov_")
    ]
    return subset.drop(columns=drop).reset_index(drop=True)


def heavy_atom_conformers(mols: list) -> list[Chem.Mol]:
    """Drop explicit hydrogens so every pool matches the heavy-atom reference.

    FlowR writes explicit hydrogens. PoseBusters copies a generated conformer
    onto the SMILES reference, and that copy fails when the atom counts differ.
    """
    heavy: list[Chem.Mol] = []
    for mol in mols:
        if mol is None or mol.GetNumConformers() == 0:
            continue
        stripped = mol_no_h(mol)
        if stripped is None or stripped.GetNumConformers() == 0:
            continue
        heavy.append(stripped)
    return heavy


def score_conformer_pool(task: tuple[str, list, str]) -> dict[str, object]:
    """Cluster the full pool, then PoseBusters-filter it for the energy stats."""
    smiles, mols, name = task
    RDLogger.DisableLog("rdApp.*")
    generated = heavy_atom_conformers(mols)
    n_generated = len(generated)
    clusters = greedy_cluster_metrics(generated)
    if n_generated == 0:
        return {
            "smiles": smiles,
            "name": name,
            "n_generated": 0,
            "n_pb_pass": 0,
            "pb_pass_rate": math.nan,
            **clusters,
            **energy_summary([]),
        }

    from casf_benchmark.generation.conformer_sets import posebusters_geometry_passes

    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        reference = Chem.Mol(generated[0])
        reference.RemoveAllConformers()
    result = posebusters_geometry_passes(generated, reference, max_workers=1, energy_num_threads=1)
    n_pass = int(sum(result.passes))
    passing = [mol for mol, ok in zip(generated, result.passes) if ok]
    return {
        "smiles": smiles,
        "name": name,
        "n_generated": n_generated,
        "n_pb_pass": n_pass,
        "pb_pass_rate": n_pass / n_generated,
        **clusters,
        **energy_summary(passing),
    }


def _column_mean(frame: pd.DataFrame, column: str) -> float:
    if column not in frame.columns or frame.empty:
        return math.nan
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().all():
        return math.nan
    return float(values.mean())


def summarize_molecules(per_molecule: pd.DataFrame) -> dict[str, float | int]:
    row: dict[str, float | int] = {"n_molecules": int(len(per_molecule))}
    for source, destination in SUMMARY_MEAN_COLUMNS:
        row[destination] = _column_mean(per_molecule, source)
    return row


def cached_rmsd_matrices(
    *,
    label: str,
    gen_pickle: Path,
    ground_truth: dict,
    cache_dir: Path,
    workers: int,
) -> dict[str, np.ndarray]:
    cache_path = cache_dir / f"{label}.npz"
    if cache_path.is_file():
        print(f"{label}: loading cached RMSD matrices from {cache_path}", flush=True)
        return load_rmsd_matrices(cache_path)
    print(f"{label}: computing RMSD matrices from {gen_pickle}", flush=True)
    with gen_pickle.open("rb") as handle:
        generated = pickle.load(handle)
    matrices, missing, all_nan = compute_rmsd_matrices(ground_truth, generated, workers)
    if missing:
        print(f"WARNING: {label} missing generated confs for {len(missing)} molecule(s)", flush=True)
    if all_nan:
        print(f"WARNING: {label} all-NaN RMSD matrices for {len(all_nan)} molecule(s)", flush=True)
    save_rmsd_matrices(cache_path, matrices)
    print(f"{label}: cached RMSD matrices -> {cache_path}", flush=True)
    return matrices


def method_dir(output_dir: Path, label: str) -> Path:
    return output_dir / label


def method_outputs_complete(path: Path) -> bool:
    return all((path / name).is_file() for name in METHOD_FILES)


def write_method_outputs(
    path: Path,
    *,
    summary: pd.DataFrame,
    per_molecule: pd.DataFrame,
    recovery: pd.DataFrame,
    k_efficiency: pd.DataFrame,
) -> None:
    path.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path / "summary.csv", index=False)
    per_molecule.to_csv(path / "per_molecule.csv", index=False)
    recovery.to_csv(path / "recovery.csv", index=False)
    k_efficiency.to_csv(path / "k_efficiency.csv", index=False)


def _read_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    return frame


def load_panel_frames(output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Concatenate every completed method directory under `output_dir`."""
    summaries: list[pd.DataFrame] = []
    molecules: list[pd.DataFrame] = []
    curves: list[pd.DataFrame] = []
    curves_k: list[pd.DataFrame] = []
    if output_dir.is_dir():
        for path in sorted(child for child in output_dir.iterdir() if child.is_dir()):
            if not method_outputs_complete(path):
                continue
            summaries.append(_read_csv(path / "summary.csv"))
            molecules.append(_read_csv(path / "per_molecule.csv"))
            curves.append(_read_csv(path / "recovery.csv"))
            curves_k.append(_read_csv(path / "k_efficiency.csv"))
    return (
        pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame(),
        pd.concat(molecules, ignore_index=True) if molecules else pd.DataFrame(),
        pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(),
        pd.concat(curves_k, ignore_index=True) if curves_k else pd.DataFrame(),
    )


def publish_plinder_tables(
    extended_db: Path,
    summary: pd.DataFrame,
    per_molecule: pd.DataFrame,
    recovery: pd.DataFrame,
    k_efficiency: pd.DataFrame,
) -> None:
    extended_db.parent.mkdir(parents=True, exist_ok=True)
    frames = {
        "plinder_summary": summary,
        "plinder_per_molecule": per_molecule,
        "plinder_recovery_curve": recovery,
        "plinder_k_efficiency": k_efficiency,
    }
    with sqlite3.connect(str(extended_db)) as connection:
        for table, frame in frames.items():
            frame.to_sql(table, connection, if_exists="replace", index=False)
        connection.commit()
