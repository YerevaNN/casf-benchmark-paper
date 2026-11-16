#!/usr/bin/env python3
"""Generate reference-independent PLINDER torsion baselines.

Each molecule is built from its isomeric SMILES. One ETKDGv3 embedding,
seeded from the global seed and the canonical SMILES, is torsion-perturbed
to 1,000 accepted conformers. The minimized pool is a fresh torsion sample
from that same embedding, minimized with MMFF94s and filtered again for
clashes, until 1,000 survivors remain.

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
    MINIMIZED_PROTOCOL,
    N_ETKDG_STARTS,
    RAW_METHOD,
    CONFORMERS_PER_START,
    annotate_pickle_hashes,
    generate_from_reference_row,
    mmff_status_counts,
    provenance_document,
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


def _merge_minimized_provenance(generated: object, provenance_dir: Path) -> None:
    """Record the fresh minimized pool without discarding the existing raw record."""
    provenance_dir.mkdir(parents=True, exist_ok=True)
    from casf_benchmark.generation.plinder_independent_torsion import _slug, _write_sdf, sha256_file

    slug = _slug(generated.molecule_id)
    path = provenance_dir / f"{slug}.json"
    document = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else provenance_document(generated)
    document["minimized_conformer_count"] = len(generated.minimized_mols)
    document["minimized_protocol"] = MINIMIZED_PROTOCOL
    document["minimized_sampling"] = dict(generated.minimized_sampling)
    document["mmff_status_counts"] = mmff_status_counts(generated.mmff_statuses)
    document["mmff_statuses"] = generated.mmff_statuses
    document["status"] = generated.status
    document["experimental_coordinates_used"] = False
    if generated.minimized_mols:
        minimized_sdf = provenance_dir / f"{slug}_minimized.sdf"
        _write_sdf(generated.minimized_mols, minimized_sdf)
        hashes = dict(document.get("output_hashes") or {})
        hashes["minimized_sdf"] = sha256_file(minimized_sdf)
        document["output_hashes"] = hashes
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def _worker(payload: tuple[object, ...]) -> object:
    RDLogger.DisableLog("rdApp.*")
    smiles, molecule_id, global_seed, n_starts, per_start, max_embed_attempts, pools = payload
    return generate_from_reference_row(
        str(smiles),
        None,
        global_seed=int(global_seed),
        molecule_id=str(molecule_id),
        n_starts=int(n_starts),
        per_start=int(per_start),
        max_embed_attempts=int(max_embed_attempts),
        pools=str(pools),
    )


def _summary(method: str, seed: int, n_starts: int, per_start: int, rows: list[object]) -> dict[str, object]:
    return {
        "method": method,
        "seed": seed,
        "seed_source": "rdkit_etkdg_v3",
        "topology_source": "isomeric_smiles",
        "experimental_coordinates_used": False,
        "n_etkdg_starts": n_starts,
        "conformers_per_start": per_start,
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
    parser.add_argument("--pools", choices=("both", "raw", "minimized"), default="both")
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    ground_truth = load_reference_sdf_dir(args.reference_dir)
    jobs = []
    references = []
    for smiles, row in ground_truth.items():
        name = str((row.get("names") or [smiles])[0])
        jobs.append((smiles, name, args.seed, args.starts, args.per_start, args.max_embed_attempts, args.pools))
        references.append(row["confs"][0])
    print(f"Generating {len(jobs)} molecules, {args.pools}, {args.starts} x {args.per_start}", flush=True)

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
        if args.pools == "minimized":
            _merge_minimized_provenance(generated, provenance_dir)
        else:
            write_molecule_provenance(generated, provenance_dir)

    raw_pickle = args.output_root / RAW_METHOD / "generation_results.pickle"
    minimized_pickle = args.output_root / MINIMIZED_METHOD / "generation_results.pickle"
    written = []
    if args.pools != "minimized":
        _atomic_pickle_dump(
            OrderedDict((row.canonical_isomeric_smiles, row.raw_mols) for row in results if row.ok),
            raw_pickle,
        )
        written.append((RAW_METHOD, raw_pickle))
    if args.pools != "raw":
        _atomic_pickle_dump(
            OrderedDict((row.canonical_isomeric_smiles, row.minimized_mols) for row in results if row.ok),
            minimized_pickle,
        )
        written.append((MINIMIZED_METHOD, minimized_pickle))
    for method, path in written:
        summary_path = path.parent / "generation_summary.json"
        summary_path.write_text(
            json.dumps(_summary(method, args.seed, args.starts, args.per_start, results), indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {path}", flush=True)
    annotate_pickle_hashes(provenance_dir, raw_pickle, minimized_pickle)
    failed = [row.molecule_id for row in results if not row.ok]
    if failed:
        raise SystemExit(f"Generation failed for: {', '.join(failed)}")


if __name__ == "__main__":
    main()
