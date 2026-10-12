#!/usr/bin/env python3
"""Generate the CASF raw RDKit / torsion baselines on the 23-molecule druglike set.

Mirrors the `rdkit_random_raw_fixed` and `torsion_raw_fixed` pipelines in
casf_benchmark.generation.conformer_sets, minus the CASF-only parts (dynamic /
ChEMBL-count subsets and the PoseBusters filter -- the druglike evaluators score
PB validity themselves, so pools are written unfiltered like every other
druglike generator):

  - rdkit_random_raw: ETKDGv3 with random coordinates, no clash filter, no
    minimization.
  - torsion_raw: rotatable torsions of a seed structure perturbed by up to
    +/-120 deg, kept if they pass the DG-bounds clash check at 0.7.

The druglike set has no single reference pose, so the torsion seed (and the
RDKit embedding topology) is the molecule's first crystal conformer from the
druglike pickle, with hydrogens added -- as CASF seeds from the ChEMBL3D
reference structure.

Writes `{output-root}/{method}_druglike/generation_results.pickle`
(SMILES -> list of single-conformer, heavy-atom Chem.Mol), the contract the
druglike evaluators read.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import tempfile
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger

from casf_benchmark.chembl3d.loader import prepare_torsion_ref_mol
from casf_benchmark.generation.conformer_sets import (
    accumulate_pre_clash_pool_rdkit,
    accumulate_pre_clash_pool_torsion,
    build_arg_parser as build_casf_arg_parser,
    clone_args,
    get_rotatable_torsions,
)

DEFAULT_DRUGLIKE_PICKLE = Path("/mnt/weka/vtarasov/druglike_smi.pickle")
METHODS = ("rdkit_random_raw", "torsion_raw")


def _atomic_pickle_dump(value: object, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=destination.parent, prefix=f".{destination.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, destination)


def load_seeds(path: Path) -> list[tuple[str, Chem.Mol]]:
    with path.open("rb") as handle:
        data = pickle.load(handle)
    seeds = []
    for smiles, row in data.items():
        seed = prepare_torsion_ref_mol(row["confs"][0])
        if seed is None:
            raise ValueError(f"No usable crystal conformer for {smiles!r}")
        seeds.append((smiles, seed))
    return seeds


def generate(
    method: str, index: int, smiles: str, seed_mol: Chem.Mol, casf_args: argparse.Namespace
) -> tuple[str, int, str, list[Chem.Mol], dict]:
    RDLogger.DisableLog("rdApp.*")
    mol_id = f"druglike_{index:03d}"
    if method == "rdkit_random_raw":
        base_mol = Chem.Mol(seed_mol)
        base_mol.RemoveAllConformers()
        pool, stats, status = accumulate_pre_clash_pool_rdkit(
            base_mol=base_mol,
            mol_id=mol_id,
            target=casf_args.fixed_set_size,
            seed=casf_args.seed,
            num_threads=1,
            args=clone_args(casf_args, pre_clash_cutoff=False, post_clash_cutoff=False),
        )
    elif not get_rotatable_torsions(seed_mol):
        return method, index, smiles, [], {"status": "no_rotatable_bonds"}
    else:
        pool, stats, status = accumulate_pre_clash_pool_torsion(
            ref_mol=seed_mol,
            mol_id=mol_id,
            target=casf_args.fixed_set_size,
            seed=casf_args.seed,
            args=clone_args(casf_args, pre_clash_cutoff=0.7, minimize_workers=1),
            seed_key="torsion_raw",
        )
    mols = [Chem.RemoveHs(mol) for mol, _energy, _ff in pool]
    return method, index, smiles, mols, {
        "status": status,
        "generated_candidates": stats.generated_candidates,
        "clash_rejected": stats.clash_rejected,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--druglike-pickle", type=Path, default=DEFAULT_DRUGLIKE_PICKLE)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--num-conformers", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()

    # Start from the CASF generator's own defaults (torsion delta, perturb
    # fraction, batch size) so both pipelines stay in lockstep with CASF.
    casf_args = build_casf_arg_parser().parse_args(
        ["--seed", str(args.seed), "--fixed_set_size", str(args.num_conformers)]
    )
    seeds = load_seeds(args.druglike_pickle)

    results: dict[str, dict[int, tuple[str, list[Chem.Mol], dict]]] = {m: {} for m in args.methods}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(generate, method, index, smiles, seed_mol, casf_args)
            for method in args.methods
            for index, (smiles, seed_mol) in enumerate(seeds)
        ]
        for future in futures:
            method, index, smiles, mols, info = future.result()
            results[method][index] = (smiles, mols, info)
            print(f"[{method}] {index:03d} n={len(mols)} {info}", flush=True)

    for method in args.methods:
        output_dir = args.output_root / f"{method}_druglike"
        rows = [results[method][index] for index in range(len(seeds))]
        _atomic_pickle_dump(
            OrderedDict((smiles, mols) for smiles, mols, _info in rows),
            output_dir / "generation_results.pickle",
        )
        summary = {
            "method": method,
            "seed": args.seed,
            "n_molecules": len(rows),
            "total_conformers": sum(len(mols) for _smiles, mols, _info in rows),
            "molecules": [
                {"molecule_index": index, "smiles": smiles, "n_conformers": len(mols), **info}
                for index, (smiles, mols, info) in enumerate(rows)
            ],
        }
        (output_dir / "generation_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        print(f"Wrote {output_dir / 'generation_results.pickle'}", flush=True)


if __name__ == "__main__":
    main()
