#!/usr/bin/env python3
"""Evaluate generated conformer pools against the "druglike" test set.

The druglike set is a small hand-picked list of named drugs/candidates with an
ensemble of experimentally observed conformers per molecule. There is no single
crystal pose, so ensemble recovery is handled by eval_druglike_covmat.py while
this script scores each generated pool on validity/diversity/energy:

  - PoseBusters validity via `posebusters_geometry_passes`. The reference molecule
    is the SMILES topology with no coordinates, so the checks fall back to the
    first generated conformer for the identity/geometry reference -- the same
    fallback the CASF generation pipeline uses when no external reference exists.
  - Diversity (torsion std, pairwise RMSD, greedy clustering) via
    `diversity_metrics`. This dominates the runtime, so molecules are evaluated in
    parallel processes, mirroring how the CASF analyzer parallelizes across ligands.
  - Energy (MMFF94s) via `energy_stats`, over the PoseBusters-passing subset only.

Input: one molgen3D `generation_results.pickle` (dict SMILES -> list[Chem.Mol]).
Output: a per-molecule CSV and a per-checkpoint summary CSV under --out-dir.
"""

from __future__ import annotations

import argparse
import json
import pickle
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

from casf_benchmark.analysis.metrics import diversity_metrics, energy_stats, sample_mols
from casf_benchmark.generation.conformer_sets import posebusters_geometry_passes

DIVERSITY_SAMPLE_CAP = 200


def _init_worker() -> None:
    RDLogger.DisableLog("rdApp.*")


def load_druglike_metadata(path: Path) -> dict[str, dict]:
    with open(path, "rb") as f:
        data = pickle.load(f)
    out = {}
    for smi, row in data.items():
        out[smi] = {
            "name": (row.get("names") or [smi])[0],
            "category": ";".join(row.get("categories") or []),
            "num_heavy_atoms": row.get("num_heavy_atoms"),
            "rotatable_bonds": (row.get("properties") or {}).get("rotb"),
            "mw": (row.get("properties") or {}).get("mw"),
        }
    return out


def _eval_one_molecule(
    label: str,
    smi: str,
    mols: list,
    info: dict,
    pb_workers: int,
    energy_num_threads: int,
) -> dict:
    n_generated = len(mols)
    if n_generated == 0:
        return {"label": label, "smiles": smi, **info, "n_generated": 0, "n_pb_pass": 0, "pb_pass_rate": float("nan")}

    reference_mol = Chem.MolFromSmiles(smi)
    if reference_mol is None:
        reference_mol = Chem.Mol(mols[0])
        reference_mol.RemoveAllConformers()

    result = posebusters_geometry_passes(
        mols, reference_mol, max_workers=pb_workers, energy_num_threads=energy_num_threads
    )
    n_pass = sum(result.passes)
    fail_counts: Counter[str] = Counter()
    for row in result.check_rows:
        for check, ok in row.items():
            if not ok:
                fail_counts[check] += 1

    pb_passing_mols = [mol for mol, ok in zip(mols, result.passes) if ok]
    # Greedy clustering is O(n x clusters); high-rotatable-bond molecules can produce
    # hundreds of clusters out of ~1000 conformers, so subsample the pool first, as
    # analysis.metrics does with PAIRWISE_SAMPLE_CAP.
    div = diversity_metrics(sample_mols(mols, cap=DIVERSITY_SAMPLE_CAP))
    energy = energy_stats(pb_passing_mols)

    return {
        "label": label,
        "smiles": smi,
        **info,
        "n_generated": n_generated,
        "n_pb_pass": n_pass,
        "pb_pass_rate": n_pass / n_generated,
        "pb_check_fail_counts_json": json.dumps(dict(fail_counts)),
        **{f"div_{k}": v for k, v in div.items()},
        **{f"energy_{k}": v for k, v in energy.items()},
    }


def eval_one_checkpoint(
    label: str,
    pickle_path: Path,
    meta: dict[str, dict],
    molecule_workers: int,
    pb_workers_per_molecule: int,
    energy_num_threads: int,
) -> tuple[list[dict], dict]:
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)

    items = []
    for smi, mols in data.items():
        info = meta.get(smi, {"name": smi, "category": "", "num_heavy_atoms": None, "rotatable_bonds": None, "mw": None})
        items.append((smi, mols, info))

    per_molecule_rows = []
    with ProcessPoolExecutor(max_workers=molecule_workers, initializer=_init_worker) as executor:
        futures = {
            executor.submit(_eval_one_molecule, label, smi, mols, info, pb_workers_per_molecule, energy_num_threads): info["name"]
            for smi, mols, info in items
        }
        for future in as_completed(futures):
            row = future.result()
            per_molecule_rows.append(row)
            name = futures[future]
            n_generated = row["n_generated"]
            if n_generated:
                print(
                    f"[{label}] {name:25s} n={n_generated:4d} "
                    f"pb_pass={row['n_pb_pass']:4d} ({row['pb_pass_rate']:.1%})",
                    flush=True,
                )
            else:
                print(f"[{label}] {name:25s} n=0", flush=True)

    pb_check_fail_counter: Counter[str] = Counter()
    total_confs = 0
    total_pb_pass = 0
    for row in per_molecule_rows:
        total_confs += row["n_generated"]
        total_pb_pass += row.get("n_pb_pass", 0)
        fail_json = row.get("pb_check_fail_counts_json")
        if fail_json:
            for check, count in json.loads(fail_json).items():
                pb_check_fail_counter[check] += count

    summary = {
        "label": label,
        "n_molecules": len(per_molecule_rows),
        "total_confs": total_confs,
        "total_pb_pass": total_pb_pass,
        "overall_pb_pass_rate": total_pb_pass / total_confs if total_confs else float("nan"),
        "mean_per_molecule_pb_pass_rate": pd.Series([r["pb_pass_rate"] for r in per_molecule_rows]).mean(),
        "pb_check_fail_counts_json": json.dumps(dict(pb_check_fail_counter)),
    }
    return per_molecule_rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--druglike-pickle",
        type=Path,
        required=True,
        help="Druglike test-set pickle (SMILES -> metadata) defining the molecules to score",
    )
    parser.add_argument("--label", required=True, help="Checkpoint label for this run")
    parser.add_argument("--gen-pickle", type=Path, required=True, help="Path to generation_results.pickle")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--molecule-workers", type=int, default=16, help="Parallelism across molecules")
    parser.add_argument("--pb-workers-per-molecule", type=int, default=1)
    parser.add_argument("--energy-num-threads", type=int, default=1)
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    meta = load_druglike_metadata(args.druglike_pickle)
    per_molecule_rows, summary = eval_one_checkpoint(
        args.label,
        args.gen_pickle,
        meta,
        args.molecule_workers,
        args.pb_workers_per_molecule,
        args.energy_num_threads,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(per_molecule_rows).to_csv(args.out_dir / f"{args.label}_per_molecule.csv", index=False)
    pd.DataFrame([summary]).to_csv(args.out_dir / f"{args.label}_summary.csv", index=False)
    print(f"Wrote {args.out_dir / f'{args.label}_per_molecule.csv'}")
    print(f"Wrote {args.out_dir / f'{args.label}_summary.csv'}")


if __name__ == "__main__":
    main()
