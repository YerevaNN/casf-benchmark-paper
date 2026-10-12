#!/usr/bin/env python3
"""Compute GEOM-style COV/MAT metrics for a generated druglike conformer pool.

Both inputs use the established SMILES-keyed pickle contract. Ground-truth
values are metadata dictionaries containing ``confs``; generated values are
lists of RDKit molecules. Outputs match the files consumed by
``build_druglike_covmat.py``.
"""

from __future__ import annotations

import argparse
import pickle
import time
from pathlib import Path

import pandas as pd
from rdkit import RDLogger

from casf_benchmark.analysis.druglike_covmat import (
    aggregate,
    compute_rmsd_matrices,
    per_molecule_metrics,
)


def write_report(
    path: Path,
    *,
    gen_pickle: Path,
    gt_pickle: Path,
    ground_truth: dict,
    generated: dict,
    missing: list[str],
    all_nan: list[str],
    metrics: dict[str, float],
    threshold: float,
    dmax: float,
    elapsed: float,
) -> None:
    successful = len(ground_truth) - len(set(missing) | set(all_nan))
    success_rate = successful / len(ground_truth) if ground_truth else float("nan")
    text = f"""\
================================================================================
COVMAT EVALUATION RESULTS
================================================================================

EVALUATION SUMMARY
----------------------------------------
Processed file: {gen_pickle}
Ground truth file: {gt_pickle}
Total molecules generated: {len(generated)}
Total conformers generated: {sum(len(value) for value in generated.values())}
Total molecules in ground truth: {len(ground_truth)}
Total conformers in ground truth: {sum(len(value.get("confs", [])) for value in ground_truth.values())}
Missing molecules (no conformers): {len(missing)}
All-NaN RMSD keys: {len(all_nan)}

EXECUTION RUNTIME
----------------------------------------
Total processing: {elapsed / 60:.4f} min
CovMat processing: {elapsed / 60:.4f} min
PoseBusters processing: 0 min

COVERAGE AND RECALL METRICS
----------------------------------------
Threshold: {threshold:g}
Coverage-Recall (COV-R):
  Mean:   {metrics["cov_r_mean"]:.4f}
  Median: {metrics["cov_r_median"]:.4f}
Coverage-Precision (COV-P):
  Mean:   {metrics["cov_p_mean"]:.4f}
  Median: {metrics["cov_p_median"]:.4f}
Matching-Recall (MAT-R):
  Mean:   {metrics["mat_r_mean"]:.4f}
  Median: {metrics["mat_r_median"]:.4f}
Matching-Precision (MAT-P):
  Mean:   {metrics["mat_p_mean"]:.4f}
  Median: {metrics["mat_p_median"]:.4f}

Censored matching (restricted mean, dmax={dmax:g} A)
Molecule success rate: {success_rate:.4f}
Censored Matching-Recall (cMAT-R):
  Mean:   {metrics["cmat_r_mean"]:.4f}
  Median: {metrics["cmat_r_median"]:.4f}
Censored Matching-Precision (cMAT-P):
  Mean:   {metrics["cmat_p_mean"]:.4f}
  Median: {metrics["cmat_p_median"]:.4f}
"""
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth-pickle", type=Path, required=True)
    parser.add_argument("--gen-pickle", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--num-workers", type=int, default=80)
    parser.add_argument("--threshold", type=float, default=0.75)
    parser.add_argument("--dmax", type=float, default=3.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    RDLogger.DisableLog("rdApp.*")
    with args.ground_truth_pickle.open("rb") as handle:
        ground_truth = pickle.load(handle)
    with args.gen_pickle.open("rb") as handle:
        generated = pickle.load(handle)

    started = time.perf_counter()
    matrices, missing, all_nan = compute_rmsd_matrices(
        ground_truth, generated, workers=args.num_workers
    )
    rows = []
    for smiles, row in ground_truth.items():
        if smiles not in matrices:
            continue
        metrics = per_molecule_metrics(matrices[smiles], args.threshold, args.dmax)
        rows.append({"geom_smiles": smiles, **metrics, "sub_smiles": ""})
    summary = aggregate(rows)
    elapsed = time.perf_counter() - started

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out_dir / "rmsd_matrix.csv", index=False)
    with (args.out_dir / "rmsd_matrix.pickle").open("wb") as handle:
        pickle.dump(
            {"threshold": args.threshold, "dmax": args.dmax, "summary": summary},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    write_report(
        args.out_dir / "covmat_results.txt",
        gen_pickle=args.gen_pickle,
        gt_pickle=args.ground_truth_pickle,
        ground_truth=ground_truth,
        generated=generated,
        missing=missing,
        all_nan=all_nan,
        metrics=summary,
        threshold=args.threshold,
        dmax=args.dmax,
        elapsed=elapsed,
    )
    print(f"Wrote COV/MAT outputs to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
