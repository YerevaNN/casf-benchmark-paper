#!/usr/bin/env python3
"""Generate reference-independent PLINDER torsion baselines.

Each molecule is built from its isomeric SMILES. Five ETKDGv3 embeddings,
seeded from the global seed, the canonical SMILES, and the seed index, are
torsion-perturbed to 200 accepted conformers each (1,000 raw conformers).
The minimized pool is the MMFF94s optimization of those same conformers.

Experimental coordinates are never used as a seed. Writes:

    {output-root}/torsion_independent_raw_druglike/generation_results.pickle
    {output-root}/torsion_independent_minimized_druglike/generation_results.pickle
    {output-root}/provenance/{molecule}.json

The older crystal-seeded pools from scripts/run_baseline_druglike.py are left
in place.
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

from rdkit import RDLogger

from casf_benchmark.analysis.plinder_panel import DEFAULT_REFERENCE_DIR, load_reference_sdf_dir
from casf_benchmark.generation.plinder_independent_torsion import (
    MINIMIZED_METHOD,
    N_ETKDG_STARTS,
    RAW_METHOD,
    CONFORMERS_PER_START,
    annotate_pickle_hashes,
    generate_from_reference_row,
    reference_topology_matches,
    write_molecule_provenance,
)

DEFAULT_OUTPUT_ROOT = Path("data/results/generation/plinder_independent_torsion")


def _atomic_pickle_dump(value: object, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=destination.parent, prefix=f".{destination.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, destination)


def _worker(payload: tuple[object, ...]) -> object:
    RDLogger.DisableLog("rdApp.*")
    smiles, molecule_id, global_seed, n_starts, per_start, max_embed_attempts = payload
    return generate_from_reference_row(
        str(smiles),
        None,
        global_seed=int(global_seed),
        molecule_id=str(molecule_id),
        n_starts=int(n_starts),
        per_start=int(per_start),
        max_embed_attempts=int(max_embed_attempts),
    )


def _summary(method: str, seed: int, rows: list[object]) -> dict[str, object]:
    return {
        "method": method,
        "seed": seed,
        "seed_source": "rdkit_etkdg_v3",
        "topology_source": "isomeric_smiles",
        "experimental_coordinates_used": False,
        "n_etkdg_starts": N_ETKDG_STARTS,
        "conformers_per_start": CONFORMERS_PER_START,
        "n_molecules": len(rows),
        "n_success": sum(1 for row in rows if row.ok),
        "total_conformers": sum(len(row.raw_mols if method == RAW_METHOD else row.minimized_mols) for row in rows),
        "molecules": [
            {
                "molecule_id": row.molecule_id,
                "smiles": row.canonical_isomeric_smiles,
                "status": row.status,
                "n_raw": len(row.raw_mols),
                "n_minimized": len(row.minimized_mols),
                "etkdg_seeds": row.etkdg_seeds,
            }
            for row in rows
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--starts", type=int, default=N_ETKDG_STARTS)
    parser.add_argument("--per-start", type=int, default=CONFORMERS_PER_START)
    parser.add_argument("--max-embed-attempts", type=int, default=8)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    ground_truth = load_reference_sdf_dir(args.reference_dir)
    jobs = []
    references = []
    for smiles, row in ground_truth.items():
        name = str((row.get("names") or [smiles])[0])
        jobs.append((smiles, name, args.seed, args.starts, args.per_start, args.max_embed_attempts))
        references.append(row["confs"][0])
    print(f"Generating {len(jobs)} molecules, {args.starts} x {args.per_start}", flush=True)

    results = []
    workers = max(1, min(int(args.workers), len(jobs) or 1))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for generated, reference in zip(executor.map(_worker, jobs), references, strict=True):
            if generated.canonical_isomeric_smiles:
                generated.reference_topology_matches = reference_topology_matches(
                    reference,
                    generated.canonical_isomeric_smiles,
                )
            results.append(generated)
            print(
                f"[{generated.molecule_id}] status={generated.status} "
                f"raw={len(generated.raw_mols)} minimized={len(generated.minimized_mols)}",
                flush=True,
            )

    provenance_dir = args.output_root / "provenance"
    for generated in results:
        write_molecule_provenance(generated, provenance_dir)

    raw_pickle = args.output_root / RAW_METHOD / "generation_results.pickle"
    minimized_pickle = args.output_root / MINIMIZED_METHOD / "generation_results.pickle"
    _atomic_pickle_dump(
        OrderedDict((row.canonical_isomeric_smiles, row.raw_mols) for row in results if row.ok),
        raw_pickle,
    )
    _atomic_pickle_dump(
        OrderedDict((row.canonical_isomeric_smiles, row.minimized_mols) for row in results if row.ok),
        minimized_pickle,
    )
    for method, path in ((RAW_METHOD, raw_pickle), (MINIMIZED_METHOD, minimized_pickle)):
        summary_path = path.parent / "generation_summary.json"
        summary_path.write_text(json.dumps(_summary(method, args.seed, results), indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {path}", flush=True)
    annotate_pickle_hashes(provenance_dir, raw_pickle, minimized_pickle)
    failed = [row.molecule_id for row in results if not row.ok]
    if failed:
        raise SystemExit(f"Generation failed for: {', '.join(failed)}")


if __name__ == "__main__":
    main()
