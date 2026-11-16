#!/usr/bin/env python3
"""Generate raw and minimized conformer sets for CASF / ChEMBL3D intersection ligands."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import random
import re
import statistics
import time
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDistGeom, rdMolTransforms

from casf_benchmark.chembl3d.identity import parse_optional_int, validate_stereo_mapping_csv
from casf_benchmark.chembl3d.loader import resolve_torsion_generation
from casf_benchmark.paths import (
    DEFAULT_CASF16_DATA,
    DEFAULT_CHEMBL_DATASET_ROOT,
    DEFAULT_CHEMBL_MAP_CSV,
    DEFAULT_CORE_INTERSECTION_ROOT,
    DEFAULT_CORE_LIGAND_DIR,
    resolve_generation_dir,
)

DEFAULT_LIGAND_DIR = DEFAULT_CORE_LIGAND_DIR
DEFAULT_TOPOLOGY_ROOT = DEFAULT_CHEMBL_DATASET_ROOT / "topologies"
DEFAULT_ZARR_ROOT = DEFAULT_CHEMBL_DATASET_ROOT / "zarr_database"
RDKIT_MIN_PRE_POOL_SIZE = 2500
TORSION_MIN_REFILL_POOL_SIZE = 500

RAW_METHODS = (
    "rdkit_random_raw_fixed",
    "rdkit_random_raw_dynamic",
    "rdkit_random_raw_chembl_count",
    "torsion_raw_fixed",
    "torsion_raw_dynamic",
    "torsion_raw_chembl_count",
)
MINIMIZED_METHODS = (
    "rdkit_random_minimized_fixed",
    "rdkit_random_minimized_dynamic",
    "rdkit_random_minimized_chembl_count",
    "torsion_minimized_fixed",
    "torsion_minimized_dynamic",
    "torsion_minimized_chembl_count",
)
ALL_METHODS = (*RAW_METHODS, *MINIMIZED_METHODS)


def get_rotatable_torsions(mol: Chem.Mol) -> list[tuple[int, int, int, int]]:
    patt = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
    matches = mol.GetSubstructMatches(patt)

    torsions: list[tuple[int, int, int, int]] = []
    seen: set[tuple[int, int]] = set()

    for j, k in matches:
        bond_key = tuple(sorted((j, k)))
        if bond_key in seen:
            continue
        seen.add(bond_key)

        aj = mol.GetAtomWithIdx(j)
        ak = mol.GetAtomWithIdx(k)

        nj = [
            a.GetIdx()
            for a in aj.GetNeighbors()
            if a.GetIdx() != k and a.GetAtomicNum() > 1
        ]
        nk = [
            a.GetIdx()
            for a in ak.GetNeighbors()
            if a.GetIdx() != j and a.GetAtomicNum() > 1
        ]

        if not nj:
            nj = [a.GetIdx() for a in aj.GetNeighbors() if a.GetIdx() != k]
        if not nk:
            nk = [a.GetIdx() for a in ak.GetNeighbors() if a.GetIdx() != j]

        if nj and nk:
            torsions.append((nj[0], j, k, nk[0]))

    return torsions


def make_candidate_from_ref(ref_mol: Chem.Mol) -> Chem.Mol:
    cand = Chem.Mol(ref_mol)
    cand.RemoveAllConformers()
    conf = Chem.Conformer(ref_mol.GetConformer(0))
    conf.SetId(0)
    cand.AddConformer(conf, assignId=True)
    return cand


def perturb_torsions(
    mol: Chem.Mol,
    torsions: list[tuple[int, int, int, int]],
    max_delta_deg: float,
    perturb_fraction: float,
) -> None:
    conf = mol.GetConformer(0)
    for i, j, k, l in torsions:
        if random.random() > perturb_fraction:
            continue
        old_angle = rdMolTransforms.GetDihedralDeg(conf, i, j, k, l)
        delta = random.uniform(-max_delta_deg, max_delta_deg)
        rdMolTransforms.SetDihedralDeg(conf, i, j, k, l, old_angle + delta)


def get_forcefield(
    mol: Chem.Mol,
    conf_id: int = 0,
    ff_variant: str = "MMFF94s",
) -> tuple[object, str]:
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError(f"Missing MMFF parameters for {ff_variant}")
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ff_variant)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    return ff, ff_variant


def minimize(
    mol: Chem.Mol,
    max_iters: int = 500,
    ff_variant: str = "MMFF94s",
) -> tuple[int, float, str]:
    ff, ff_name = get_forcefield(mol, conf_id=0, ff_variant=ff_variant)
    status = ff.Minimize(maxIts=max_iters)
    energy = ff.CalcEnergy()
    return status, energy, ff_name


DG_BOUND_MATRIX_PARAMS = {
    "set15bounds": True,
    "scaleVDW": True,
    "doTriangleSmoothing": True,
    "useMacrocycle14config": False,
}


def _get_bond_pairs(mol: Chem.Mol) -> set[tuple[int, int]]:
    return {
        tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        for bond in mol.GetBonds()
    }


def _get_angle_pairs(mol: Chem.Mol) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()

    for atom in mol.GetAtoms():
        neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        for a in range(len(neighbors)):
            for b in range(a + 1, len(neighbors)):
                pairs.add(tuple(sorted((neighbors[a], neighbors[b]))))

    return pairs


def get_dg_bounds(
    mol: Chem.Mol,
    bound_matrix_params: dict[str, bool] = DG_BOUND_MATRIX_PARAMS,
) -> tuple[object, set[tuple[int, int]], set[tuple[int, int]]]:
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(mol, **bound_matrix_params)
    return bounds, _get_bond_pairs(mol), _get_angle_pairs(mol)


def has_clash(
    mol: Chem.Mol,
    bounds: object,
    bond_pairs: set[tuple[int, int]],
    angle_pairs: set[tuple[int, int]],
    cutoff: float,
    ignore_hydrogens: bool = True,
) -> bool:
    conf = mol.GetConformer(0)
    n = mol.GetNumAtoms()

    for i in range(n):
        if ignore_hydrogens and mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue

        for j in range(i + 1, n):
            if ignore_hydrogens and mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue

            pair = (i, j)
            if pair in bond_pairs or pair in angle_pairs:
                continue

            lower_bound = bounds[j, i]
            d = conf.GetAtomPosition(i).Distance(conf.GetAtomPosition(j))
            if d < cutoff * lower_bound:
                return True

    return False


@dataclass(frozen=True)
class InputMolecule:
    mol_id: str
    smiles: str
    source_input: str
    chembl3d_group: str = ""
    chembl3d_mol_id: str = ""
    chembl3d_conformer_count: int = 0
    chembl3d_sdf_record_index: int | None = None


@dataclass(frozen=True)
class FilterStats:
    generated_candidates: int = 0
    finite_rejected: int = 0
    clash_rejected: int = 0
    bond_rejected: int = 0
    stereo_rejected: int = 0
    rmsd_rejected: int = 0
    pre_clash_passed: int = 0
    generation_batches: int = 0


@dataclass(frozen=True)
class MethodResult:
    mol_id: str
    input_smiles: str
    source_input: str
    generation_method: str
    set_tier: str
    num_target_confs: int
    rotatable_bonds: int
    generated_candidates: int
    kept_confs: int
    finite_rejected: int
    clash_rejected: int
    bond_rejected: int
    stereo_rejected: int
    rmsd_rejected: int
    pre_clash_passed: int
    generation_batches: int
    minimization_applied: bool
    minimization_input_confs: int
    minimization_failed: int
    minimization_invalid: int
    minimization_error: int
    post_min_clash_rejected: int
    selected_confs: int
    waste_ratio: float
    minimization_failed_rate: float
    minimization_error_rate: float
    post_min_clash_rejected_rate: float
    pb_input_confs: int
    pb_pass_confs: int
    pb_fail_confs: int
    pb_pass_rate: float
    pb_check_fail_counts_json: str
    pb_check_fail_rates_json: str
    walltime_seconds: float
    status: str
    sdf_path: str


@dataclass(frozen=True)
class PoseBustersResult:
    passes: list[bool]
    check_rows: list[dict[str, bool]]


def stable_seed(base_seed: int, key: str) -> int:
    digest = hashlib.blake2b(key.encode("utf-8"), digest_size=8).digest()
    offset = int.from_bytes(digest, byteorder="little", signed=False)
    return (int(base_seed) + offset) % (2**31 - 1)


def safe_mol_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("._") or "molecule"


def parse_limit_molecules(value: str) -> int | None:
    if value.lower() == "all":
        return None
    parsed = int(value)
    if parsed < 0:
        raise ValueError("--limit_molecules must be non-negative or 'all'")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def float_or_false(value: str) -> float | bool:
    if value.strip().lower() in ("false", "off"):
        return False
    return nonnegative_float(value)


def clone_args(args: argparse.Namespace, **overrides: object) -> argparse.Namespace:
    cloned = argparse.Namespace(**vars(args))
    for key, value in overrides.items():
        setattr(cloned, key, value)
    return cloned


def safe_ratio(numerator: float, denominator: float) -> float:
    return (float(numerator) / float(denominator)) if denominator else math.nan


def json_dumps_sorted(value: dict[str, int | float]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def summarize_posebusters_checks(check_rows: list[dict[str, bool]]) -> tuple[str, str]:
    if not check_rows:
        return "{}", "{}"

    check_names = sorted({name for row in check_rows for name in row})
    fail_counts = {
        name: sum(1 for row in check_rows if row.get(name) is False)
        for name in check_names
    }
    total = len(check_rows)
    fail_rates = {name: safe_ratio(count, total) for name, count in fail_counts.items()}
    return json_dumps_sorted(fail_counts), json_dumps_sorted(fail_rates)


def dynamic_candidate_count(rotatable_bonds: int) -> int:
    return max(1, -20 + (22 * int(rotatable_bonds)))


def chembl_count_target(input_mol: InputMolecule) -> int:
    return max(1, int(input_mol.chembl3d_conformer_count))


def report_dynamic_target_cap(
    input_mol: InputMolecule,
    rotatable_bonds: int,
    dynamic_target: int,
    fixed_set_size: int,
) -> None:
    if dynamic_target <= fixed_set_size:
        return
    print(
        f"[dynamic-target-cap] {input_mol.mol_id}: requested dynamic target "
        f"{dynamic_target} (-20 + 22 * {rotatable_bonds} rotatable bonds) "
        f"exceeds fixed pool size {fixed_set_size}; dynamic subsets are capped at the fixed pool.",
        flush=True,
    )


def report_chembl_count_target_cap(
    input_mol: InputMolecule,
    chembl_target: int,
    fixed_set_size: int,
) -> None:
    if chembl_target <= fixed_set_size:
        return
    print(
        f"[chembl-count-target-cap] {input_mol.mol_id}: requested ChEMBL-count target "
        f"{chembl_target} exceeds fixed pool size {fixed_set_size}; ChEMBL-count subsets "
        "are capped at the fixed pool.",
        flush=True,
    )


def etkdg_params(seed: int, num_threads: int) -> AllChem.EmbedParameters:
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = int(num_threads)
    params.pruneRmsThresh = -1.0
    params.useRandomCoords = True
    params.enforceChirality = True
    return params


def has_finite_coordinates(mol: Chem.Mol, conf_id: int) -> bool:
    conf = mol.GetConformer(conf_id)
    for atom_idx in range(mol.GetNumAtoms()):
        point = conf.GetAtomPosition(atom_idx)
        if not (math.isfinite(point.x) and math.isfinite(point.y) and math.isfinite(point.z)):
            return False
    return True


def single_conf_mol(source: Chem.Mol, conf_id: int) -> Chem.Mol:
    single = Chem.Mol(source)
    single.RemoveAllConformers()
    conf = Chem.Conformer(source.GetConformer(conf_id))
    conf.SetId(0)
    single.AddConformer(conf, assignId=True)
    return single


def build_mol_from_candidates(
    template: Chem.Mol,
    candidates: list[tuple[Chem.Mol, float, str, float | None]],
    source_raw_conf_ids: list[int] | None = None,
) -> Chem.Mol:
    out = Chem.Mol(template)
    out.RemoveAllConformers()
    for new_conf_id, (cand, energy, ff_name, min_rmsd) in enumerate(candidates):
        conf = Chem.Conformer(cand.GetConformer(0))
        conf.SetId(new_conf_id)
        conf.SetDoubleProp("post_min_energy", float(energy))
        conf.SetProp("force_field_used", ff_name)
        if min_rmsd is not None:
            conf.SetDoubleProp("min_rmsd_to_selected", float(min_rmsd))
        if source_raw_conf_ids is not None:
            conf.SetIntProp("source_raw_conf_id", int(source_raw_conf_ids[new_conf_id]))
        out.AddConformer(conf, assignId=True)
    return out


def select_dynamic_indices(
    fixed_count: int,
    dynamic_target: int,
    seed: int,
    key: str,
) -> list[int]:
    if fixed_count <= 0:
        return []
    count = min(int(dynamic_target), int(fixed_count))
    rng = random.Random(stable_seed(seed, key))
    return sorted(rng.sample(range(fixed_count), count))


def subset_mol_by_indices(
    source: Chem.Mol,
    indices: list[int],
    template: Chem.Mol,
) -> Chem.Mol:
    out = Chem.Mol(template)
    out.RemoveAllConformers()
    for new_conf_id, source_idx in enumerate(indices):
        conf = Chem.Conformer(source.GetConformer(source_idx))
        conf.SetId(new_conf_id)
        conf.SetIntProp("fixed_pool_index", int(source_idx))
        out.AddConformer(conf, assignId=True)
    return out


def write_indices_sidecar(
    path: Path,
    mol_id: str,
    family: str,
    tier: str,
    indices: list[int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"mol_id\t{mol_id}\n")
        handle.write(f"family\t{family}\n")
        handle.write(f"tier\t{tier}\n")
        handle.write("fixed_pool_index\n")
        for index in indices:
            handle.write(f"{index}\n")


def indices_sidecar_path(sdf_path: Path) -> Path:
    return sdf_path.with_suffix(".indices.tsv")


def build_raw_mol_from_pool(
    pool: list[tuple[Chem.Mol, float, str]],
    template: Chem.Mol,
) -> Chem.Mol:
    candidates = [(mol, energy, ff_name, None) for mol, energy, ff_name in pool]
    return build_mol_from_candidates(template, candidates)


def pool_status(target: int, pool_len: int, embed_failed: bool = False) -> str:
    if embed_failed and pool_len == 0:
        return "embedding_failed"
    if pool_len < target:
        return "failed_to_fill_pool"
    return "ok"


def generate_rdkit_batch(
    base_mol: Chem.Mol,
    mol_id: str,
    batch_size: int,
    batch_index: int,
    seed: int,
    num_threads: int,
    args: argparse.Namespace,
) -> tuple[list[tuple[Chem.Mol, float, str]], int, int, int]:
    mol = Chem.Mol(base_mol)
    mol.RemoveAllConformers()
    params = etkdg_params(
        stable_seed(seed, f"{mol_id}:rdkit:batch{batch_index}"),
        num_threads,
    )
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=batch_size, params=params))
    if not conf_ids:
        return [], 0, 0, 0

    survivors: list[tuple[Chem.Mol, float, str]] = []
    finite_rejected = 0
    pre_clash_rejected = 0
    pre_clash_cutoff = args.pre_clash_cutoff
    dg_clash_data = get_dg_bounds(mol) if pre_clash_cutoff is not False else None

    for conf_id in conf_ids:
        if not has_finite_coordinates(mol, conf_id):
            finite_rejected += 1
            continue
        single = single_conf_mol(mol, conf_id)
        if pre_clash_cutoff is not False and dg_clash_data is not None and has_clash(
            single,
            *dg_clash_data,
            cutoff=float(pre_clash_cutoff),
        ):
            pre_clash_rejected += 1
            continue
        survivors.append((single, 0.0, "none"))

    return survivors, len(conf_ids), finite_rejected, pre_clash_rejected


def generate_torsion_batch(
    ref_mol: Chem.Mol,
    torsions: list[tuple[int, int, int, int]],
    mol_id: str,
    batch_size: int,
    batch_index: int,
    seed: int,
    args: argparse.Namespace,
    seed_key: str = "torsion",
    dg_clash_data: tuple[object, set[tuple[int, int]], set[tuple[int, int]]] | None = None,
) -> tuple[list[tuple[Chem.Mol, float, str]], int, int]:
    random.seed(stable_seed(seed, f"{mol_id}:{seed_key}:batch{batch_index}"))
    survivors: list[tuple[Chem.Mol, float, str]] = []
    pre_clash_rejected = 0
    pre_clash_cutoff = getattr(args, "pre_clash_cutoff", False)

    for _ in range(batch_size):
        cand = make_candidate_from_ref(ref_mol)
        perturb_torsions(
            cand,
            torsions,
            max_delta_deg=args.max_torsion_delta_deg,
            perturb_fraction=args.perturb_fraction,
        )
        if pre_clash_cutoff is not False and dg_clash_data is not None and has_clash(
            cand,
            *dg_clash_data,
            cutoff=float(pre_clash_cutoff),
        ):
            pre_clash_rejected += 1
            continue
        survivors.append((cand, 0.0, "none"))

    return survivors, batch_size, pre_clash_rejected


def _split_work(total: int, workers: int) -> list[int]:
    workers = max(1, min(int(workers), int(total)))
    base = int(total) // workers
    remainder = int(total) % workers
    return [base + (1 if idx < remainder else 0) for idx in range(workers)]


def _torsion_batch_worker(
    payload: tuple[object, ...],
) -> tuple[list[bytes], int, int]:
    (
        ref_mol_bytes,
        torsions,
        mol_id,
        batch_size,
        batch_index,
        seed,
        seed_key,
        max_torsion_delta_deg,
        perturb_fraction,
        pre_clash_cutoff,
        dg_clash_data,
    ) = payload
    ref_mol = Chem.Mol(ref_mol_bytes)
    worker_args = argparse.Namespace(
        max_torsion_delta_deg=max_torsion_delta_deg,
        perturb_fraction=perturb_fraction,
        pre_clash_cutoff=pre_clash_cutoff,
    )
    survivors, generated, pre_clash_rejected = generate_torsion_batch(
        ref_mol=ref_mol,
        torsions=torsions,
        mol_id=mol_id,
        batch_size=batch_size,
        batch_index=batch_index,
        seed=seed,
        args=worker_args,
        seed_key=seed_key,
        dg_clash_data=dg_clash_data,
    )
    survivor_bytes = [mol.ToBinary() for mol, _, _ in survivors]
    del survivors, ref_mol
    gc.collect()
    return survivor_bytes, generated, pre_clash_rejected


def generate_torsion_batch_parallel(
    ref_mol: Chem.Mol,
    torsions: list[tuple[int, int, int, int]],
    mol_id: str,
    batch_size: int,
    batch_index: int,
    seed: int,
    args: argparse.Namespace,
    seed_key: str = "torsion",
    dg_clash_data: tuple[object, set[tuple[int, int]], set[tuple[int, int]]] | None = None,
) -> tuple[list[tuple[Chem.Mol, float, str]], int, int, int]:
    workers = max(1, min(int(args.minimize_workers), int(batch_size)))
    if workers == 1:
        survivors, generated, pre_clash_rejected = generate_torsion_batch(
            ref_mol=ref_mol,
            torsions=torsions,
            mol_id=mol_id,
            batch_size=batch_size,
            batch_index=batch_index,
            seed=seed,
            args=args,
            seed_key=seed_key,
            dg_clash_data=dg_clash_data,
        )
        return survivors, generated, pre_clash_rejected, 1

    ref_mol_bytes = ref_mol.ToBinary()
    payloads = [
        (
            ref_mol_bytes,
            torsions,
            mol_id,
            chunk_size,
            batch_index + chunk_idx,
            seed,
            seed_key,
            float(args.max_torsion_delta_deg),
            float(args.perturb_fraction),
            getattr(args, "pre_clash_cutoff", False),
            dg_clash_data,
        )
        for chunk_idx, chunk_size in enumerate(_split_work(batch_size, workers))
    ]

    survivors: list[tuple[Chem.Mol, float, str]] = []
    total_generated = 0
    total_pre_clash_rejected = 0
    payload_count = len(payloads)
    chunksize = max(1, len(payloads) // max(1, workers * 4))
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_single_thread_worker) as executor:
        for survivor_bytes, generated, pre_clash_rejected in executor.map(
            _torsion_batch_worker,
            payloads,
            chunksize=chunksize,
        ):
            total_generated += generated
            total_pre_clash_rejected += pre_clash_rejected
            survivors.extend((Chem.Mol(mol_bytes), 0.0, "none") for mol_bytes in survivor_bytes)
            del survivor_bytes

    del payloads, ref_mol_bytes
    gc.collect()
    return survivors, total_generated, total_pre_clash_rejected, payload_count


def accumulate_pre_clash_pool_rdkit(
    base_mol: Chem.Mol,
    mol_id: str,
    target: int,
    seed: int,
    num_threads: int,
    args: argparse.Namespace,
) -> tuple[list[tuple[Chem.Mol, float, str]], FilterStats, str]:
    pool: list[tuple[Chem.Mol, float, str]] = []
    total_generated = 0
    finite_rejected = 0
    pre_clash_rejected = 0
    batch_index = 0
    embed_failed = False

    while len(pool) < target:
        remaining = target - len(pool)
        if args.pre_clash_cutoff is False:
            batch_size = remaining
        else:
            batch_size = int(args.generation_batch_size)
        survivors, generated, batch_finite, batch_clash = generate_rdkit_batch(
            base_mol=base_mol,
            mol_id=mol_id,
            batch_size=batch_size,
            batch_index=batch_index,
            seed=seed,
            num_threads=num_threads,
            args=args,
        )
        batch_index += 1
        total_generated += generated
        finite_rejected += batch_finite
        pre_clash_rejected += batch_clash

        if generated == 0:
            embed_failed = True
            break

        for survivor in survivors:
            pool.append(survivor)
            if len(pool) == target:
                break

    stats = FilterStats(
        generated_candidates=total_generated,
        finite_rejected=finite_rejected,
        clash_rejected=pre_clash_rejected,
        pre_clash_passed=len(pool),
        generation_batches=batch_index,
    )
    return pool, stats, pool_status(target, len(pool), embed_failed)


def accumulate_pre_clash_pool_torsion(
    ref_mol: Chem.Mol,
    mol_id: str,
    target: int,
    seed: int,
    args: argparse.Namespace,
    seed_key: str = "torsion_raw",
) -> tuple[list[tuple[Chem.Mol, float, str]], FilterStats, str]:
    torsions = get_rotatable_torsions(ref_mol)
    dg_clash_data = (
        get_dg_bounds(ref_mol)
        if getattr(args, "pre_clash_cutoff", False) is not False
        else None
    )
    pool: list[tuple[Chem.Mol, float, str]] = []
    total_generated = 0
    pre_clash_rejected = 0
    batch_index = 0

    while len(pool) < target:
        batch_size = int(args.generation_batch_size)
        survivors, generated, batch_clash, batches_used = generate_torsion_batch_parallel(
            ref_mol=ref_mol,
            torsions=torsions,
            mol_id=mol_id,
            batch_size=batch_size,
            batch_index=batch_index,
            seed=seed,
            args=args,
            seed_key=seed_key,
            dg_clash_data=dg_clash_data,
        )
        batch_index += batches_used
        total_generated += generated
        pre_clash_rejected += batch_clash

        for survivor in survivors:
            pool.append(survivor)
            if len(pool) == target:
                break
        del survivors
        gc.collect()

    stats = FilterStats(
        generated_candidates=total_generated,
        clash_rejected=pre_clash_rejected,
        pre_clash_passed=len(pool),
        generation_batches=batch_index,
    )
    status = "failed_to_fill_pool" if len(pool) < target else "ok"
    if not pool:
        status = "all_candidates_invalid"
    return pool, stats, status


def _init_single_thread_worker() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"


def _minimize_worker(payload: tuple[object, ...]) -> tuple[bytes | None, float, str, str]:
    (
        mol_bytes,
        max_iters,
        ff_variant,
        post_clash_cutoff,
        dg_clash_data,
    ) = payload
    single = Chem.Mol(mol_bytes)
    try:
        status, energy, ff_name = minimize(
            single,
            max_iters=max_iters,
            ff_variant=ff_variant,
        )
        if status not in (0, 1):
            return None, 0.0, ff_name, "min_failed"
        if post_clash_cutoff is not False and dg_clash_data is not None and has_clash(
            single,
            *dg_clash_data,
            cutoff=float(post_clash_cutoff),
        ):
            return None, 0.0, ff_name, "post_clash"
        return single.ToBinary(), energy, ff_name, "ok"
    except Exception:
        return None, 0.0, "failed", "error"


def minimize_candidates_parallel(
    candidates: list[tuple[Chem.Mol, float, str]],
    args: argparse.Namespace,
) -> tuple[list[tuple[Chem.Mol, float, str]], int, int, int]:
    if not candidates:
        return [], 0, 0, 0

    workers = max(1, min(int(args.minimize_workers), len(candidates)))
    post_clash_cutoff = args.post_clash_cutoff
    dg_clash_data = get_dg_bounds(candidates[0][0]) if post_clash_cutoff is not False else None
    payloads = [
        (
            cand[0].ToBinary(),
            int(args.max_minimize_iters),
            args.ff_variant,
            post_clash_cutoff,
            dg_clash_data,
        )
        for cand in candidates
    ]

    survivors: list[tuple[Chem.Mol, float, str]] = []
    min_failed = 0
    min_error = 0
    post_clash_rejected = 0

    chunksize = max(1, len(payloads) // max(1, workers * 4))
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_single_thread_worker) as executor:
        for mol_bytes, energy, ff_name, outcome in executor.map(
            _minimize_worker,
            payloads,
            chunksize=chunksize,
        ):
            if outcome == "ok" and mol_bytes is not None:
                survivors.append((Chem.Mol(mol_bytes), energy, ff_name))
            elif outcome == "post_clash":
                post_clash_rejected += 1
            elif outcome == "min_failed":
                min_failed += 1
            else:
                min_error += 1

    del payloads
    gc.collect()
    return survivors, min_failed, min_error, post_clash_rejected


def accumulate_torsion_minimized_pool(
    ref_mol: Chem.Mol,
    mol_id: str,
    target: int,
    seed: int,
    args: argparse.Namespace,
) -> tuple[list[tuple[Chem.Mol, float, str]], FilterStats, str, int, int, int, int]:
    torsions = get_rotatable_torsions(ref_mol)
    final_pool: list[tuple[Chem.Mol, float, str]] = []
    total_generated = 0
    pre_clash_rejected = 0
    minimization_input = 0
    min_failed = 0
    min_error = 0
    post_clash_rejected = 0
    generation_batches = 0

    args_pre = clone_args(args, pre_clash_cutoff=0.7)
    args_post = clone_args(args, post_clash_cutoff=0.7)
    dg_clash_data = get_dg_bounds(ref_mol)
    batch_index = 0

    while len(final_pool) < target:
        pre_pool: list[tuple[Chem.Mol, float, str]] = []
        pre_target = (
            int(args.torsion_min_pre_pool_size)
            if not final_pool
            else TORSION_MIN_REFILL_POOL_SIZE
        )

        while len(pre_pool) < pre_target:
            batch_size = int(args.generation_batch_size)
            survivors, generated, batch_clash, batches_used = generate_torsion_batch_parallel(
                ref_mol=ref_mol,
                torsions=torsions,
                mol_id=mol_id,
                batch_size=batch_size,
                batch_index=batch_index,
                seed=seed,
                args=args_pre,
                seed_key="torsion_min_pre",
                dg_clash_data=dg_clash_data,
            )
            batch_index += batches_used
            generation_batches += batches_used
            total_generated += generated
            pre_clash_rejected += batch_clash
            for survivor in survivors:
                pre_pool.append(survivor)
                if len(pre_pool) >= pre_target:
                    break
            del survivors
            gc.collect()

        minimization_input += len(pre_pool)
        minimized, batch_min_failed, batch_min_error, batch_post_rejected = minimize_candidates_parallel(
            pre_pool,
            args_post,
        )
        min_failed += batch_min_failed
        min_error += batch_min_error
        post_clash_rejected += batch_post_rejected

        for item in minimized:
            final_pool.append(item)
            if len(final_pool) == target:
                break
        del pre_pool, minimized
        gc.collect()

    stats = FilterStats(
        generated_candidates=total_generated,
        clash_rejected=pre_clash_rejected,
        pre_clash_passed=len(final_pool),
        generation_batches=generation_batches,
    )
    status = "failed_to_fill_pool" if len(final_pool) < target else "ok"
    return final_pool, stats, status, minimization_input, min_failed, min_error, post_clash_rejected


_POSEBUSTERS_RUNNERS: dict[tuple[int, int | None, int], object] = {}


_POSEBUSTERS_CONFORMER_VALIDITY_CONFIG = {
    "modules": [
        {
            "name": "Loading",
            "function": "loading",
            "chosen_binary_test_output": ["mol_pred_loaded", "mol_true_loaded"],
            "rename_outputs": {
                "mol_pred_loaded": "File loads",
                "mol_true_loaded": "Reference file loads",
            },
        },
        {
            "name": "Chemistry",
            "function": "rdkit_sanity",
            "chosen_binary_test_output": ["passes_rdkit_sanity_checks"],
            "rename_outputs": {"passes_rdkit_sanity_checks": "Sanitisation"},
        },
        {
            "name": "Chemistry",
            "function": "inchi_convertible",
            "chosen_binary_test_output": ["inchi_convertible"],
            "rename_outputs": {"inchi_convertible": "InChI convertible"},
        },
        {
            "name": "Chemistry",
            "function": "atoms_connected",
            "chosen_binary_test_output": ["all_atoms_connected"],
            "rename_outputs": {"all_atoms_connected": "All atoms connected"},
        },
        {
            "name": "Chemistry",
            "function": "check_radicals",
            "chosen_binary_test_output": ["no_radicals"],
            "rename_outputs": {"no_radicals": "No radicals"},
        },
        {
            "name": "Chemistry",
            "function": "identity",
            "parameters": {"inchi_options": "w"},
            "chosen_binary_test_output": [
                "formula",
                "connections",
                "stereo_tetrahedral",
                "stereo_dbond",
            ],
            "rename_outputs": {
                "formula": "Molecular formula",
                "connections": "Bonds",
                "stereo_tetrahedral": "Tetrahedral chirality",
                "stereo_dbond": "Double bond stereochemistry",
            },
        },
        {
            "name": "Geometry",
            "function": "distance_geometry",
            "parameters": {
                "bound_matrix_params": {
                    "set15bounds": True,
                    "scaleVDW": True,
                    "doTriangleSmoothing": True,
                    "useMacrocycle14config": False,
                },
                "threshold_bad_bond_length": 0.25,
                "threshold_bad_angle": 0.25,
                "threshold_clash": 0.3,
                "ignore_hydrogens": True,
                "sanitize": True,
            },
            "chosen_binary_test_output": [
                "bond_lengths_within_bounds",
                "bond_angles_within_bounds",
                "no_internal_clash",
            ],
            "rename_outputs": {
                "bond_lengths_within_bounds": "Bond lengths",
                "bond_angles_within_bounds": "Bond angles",
                "no_internal_clash": "Internal steric clash",
            },
        },
        {
            "name": "Ring flatness",
            "function": "flatness",
            "parameters": {
                "flat_systems": {
                    "aromatic_5_membered_rings_sp2": "[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1",
                    "aromatic_6_membered_rings_sp2": "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1",
                },
                "threshold_flatness": 0.25,
            },
            "chosen_binary_test_output": ["flatness_passes"],
            "rename_outputs": {
                "flatness_passes": "Planar aromatic rings",
                "num_systems_checked": "number_aromatic_rings_checked",
                "num_systems_passed": "number_aromatic_rings_pass",
                "max_distance": "aromatic_ring_maximum_distance_from_plane",
            },
        },
        {
            "name": "Double bond flatness",
            "function": "flatness",
            "parameters": {
                "flat_systems": {
                    "trigonal_planar_double_bonds": "[C;X3;^2](*)(*)=[C;X3;^2](*)(*)",
                },
                "threshold_flatness": 0.25,
            },
            "chosen_binary_test_output": ["flatness_passes"],
            "rename_outputs": {
                "flatness_passes": "Planar double bonds",
                "num_systems_checked": "number_double_bonds_checked",
                "num_systems_passed": "number_double_bonds_pass",
                "max_distance": "double_bond_maximum_distance_from_plane",
            },
        },
        {
            "name": "Energy ratio",
            "function": "energy_ratio",
            "parameters": {
                "threshold_energy_ratio": 100.0,
                "ensemble_number_conformations": 50,
                "inchi_strict": False,
            },
            "chosen_binary_test_output": ["energy_ratio_passes"],
            "rename_outputs": {"energy_ratio_passes": "Energy ratio"},
        },
    ],
    "loading": {
        "mol_pred": {"cleanup": False, "sanitize": False, "add_hs": False, "assign_stereo": False, "load_all": True},
        "mol_true": {"cleanup": False, "sanitize": False, "add_hs": False, "assign_stereo": False, "load_all": True},
        "mol_cond": {
            "cleanup": False,
            "sanitize": False,
            "add_hs": False,
            "assign_stereo": False,
            "proximityBonding": False,
        },
    },
}


def _posebusters_config(energy_num_threads: int) -> dict[str, object]:
    config = deepcopy(_POSEBUSTERS_CONFORMER_VALIDITY_CONFIG)
    for module in config["modules"]:
        if module.get("function") == "energy_ratio":
            module.setdefault("parameters", {})["num_threads"] = int(energy_num_threads)
    return config


def _posebusters_runner(
    max_workers: int,
    chunk_size: int | None,
    energy_num_threads: int,
) -> object:
    from posebusters import PoseBusters

    cache_key = (int(max_workers), chunk_size, int(energy_num_threads))
    runner = _POSEBUSTERS_RUNNERS.get(cache_key)
    if runner is None:
        runner = PoseBusters(
            config=_posebusters_config(energy_num_threads),
            max_workers=max_workers,
            chunk_size=chunk_size,
        )
        _POSEBUSTERS_RUNNERS[cache_key] = runner
    return runner


def _posebusters_reference_mol(reference_mol: Chem.Mol, fallback_mol: Chem.Mol) -> Chem.Mol:
    reference = Chem.Mol(reference_mol)
    if reference.GetNumConformers() == 0 and fallback_mol.GetNumConformers() > 0:
        conf = Chem.Conformer(fallback_mol.GetConformer(0))
        conf.SetId(0)
        reference.AddConformer(conf, assignId=True)
    return reference


def _posebusters_failure_row(exc: Exception) -> dict[str, bool]:
    return {
        "posebusters_runtime_error": False,
        f"posebusters_exception_{type(exc).__name__}": False,
    }


def _posebusters_result_from_frame(frame: object) -> PoseBustersResult:
    bool_df = frame.select_dtypes(include=["bool"])
    if bool_df.empty:
        raise RuntimeError("PoseBusters returned no boolean conformer-validity checks.")
    check_rows = [
        {str(column): bool(value) for column, value in row.items()}
        for row in bool_df.to_dict(orient="records")
    ]
    passes = [all(row.values()) for row in check_rows]
    return PoseBustersResult(passes, check_rows)


def _run_posebusters_records(
    records: list[Chem.Mol],
    reference_mol: Chem.Mol,
    max_workers: int,
    chunk_size: int | None,
    energy_num_threads: int,
) -> PoseBustersResult:
    runner = _posebusters_runner(
        max_workers=max_workers,
        chunk_size=chunk_size,
        energy_num_threads=energy_num_threads,
    )
    frame = runner.bust(
        records,
        _posebusters_reference_mol(reference_mol, records[0]),
        None,
        full_report=False,
    )
    return _posebusters_result_from_frame(frame)


def _posebusters_isolating_fallback(
    records: list[Chem.Mol],
    reference_mol: Chem.Mol,
    max_workers: int,
    energy_num_threads: int,
) -> PoseBustersResult:
    if not records:
        return PoseBustersResult([], [])

    workers = max(1, min(int(max_workers), len(records)))
    energy_threads = 1 if workers > 1 else max(1, int(energy_num_threads))
    try:
        result = _run_posebusters_records(
            records,
            reference_mol,
            max_workers=workers,
            chunk_size=max(1, math.ceil(len(records) / max(1, workers * 4))),
            energy_num_threads=energy_threads,
        )
    except ImportError:
        raise
    except Exception as exc:  # noqa: BLE001 - split to preserve good conformers.
        if len(records) == 1:
            return PoseBustersResult([False], [_posebusters_failure_row(exc)])
        midpoint = len(records) // 2
        left = _posebusters_isolating_fallback(
            records[:midpoint],
            reference_mol,
            max_workers=max_workers,
            energy_num_threads=energy_num_threads,
        )
        right = _posebusters_isolating_fallback(
            records[midpoint:],
            reference_mol,
            max_workers=max_workers,
            energy_num_threads=energy_num_threads,
        )
        return PoseBustersResult(left.passes + right.passes, left.check_rows + right.check_rows)

    if len(result.passes) != len(records) or len(result.check_rows) != len(records):
        raise RuntimeError("PoseBusters fallback returned wrong row count.")
    return result


def posebusters_geometry_passes(
    mols: list[Chem.Mol],
    reference_mol: Chem.Mol,
    max_workers: int,
    energy_num_threads: int,
) -> PoseBustersResult:
    if not mols:
        return PoseBustersResult([], [])

    records = []
    for mol in mols:
        record = Chem.Mol(mol)
        record.RemoveAllConformers()
        record.AddConformer(Chem.Conformer(mol.GetConformer(0)), assignId=True)
        records.append(record)

    workers = max(1, min(int(max_workers), len(records)))
    energy_threads = 1 if workers > 1 else max(1, int(energy_num_threads))
    chunk_size = max(1, math.ceil(len(records) / max(1, workers * 4)))

    try:
        return _run_posebusters_records(
            records,
            reference_mol,
            max_workers=workers,
            chunk_size=chunk_size,
            energy_num_threads=energy_threads,
        )
    except ImportError as exc:
        raise RuntimeError(
            "PoseBusters is required for CASF conformer generation; refusing to mark "
            "conformers as passing when the dependency is unavailable."
        ) from exc
    except Exception:
        return _posebusters_isolating_fallback(
            records,
            reference_mol,
            max_workers=workers,
            energy_num_threads=energy_num_threads,
        )


def filter_mol_by_passes(source: Chem.Mol, passes: list[bool]) -> Chem.Mol:
    out = Chem.Mol(source)
    out.RemoveAllConformers()
    new_id = 0
    for conf_idx, conf in enumerate(source.GetConformers()):
        if conf_idx >= len(passes) or not passes[conf_idx]:
            continue
        new_conf = Chem.Conformer(conf)
        new_conf.SetId(new_id)
        out.AddConformer(new_conf, assignId=True)
        new_id += 1
    return out


def set_record_props(
    record: Chem.Mol,
    input_mol: InputMolecule,
    method: str,
    conf_id: int,
    rotatable_bonds: int,
    set_tier: str,
    num_target_confs: int,
    minimization_applied: bool,
    conf: Chem.Conformer,
) -> None:
    for prop in list(record.GetPropNames()):
        record.ClearProp(prop)
    record.SetProp("_Name", f"{input_mol.mol_id}_{method}_conf{conf_id}")
    record.SetProp("mol_id", input_mol.mol_id)
    record.SetProp("input_smiles", input_mol.smiles)
    record.SetProp("source_input", input_mol.source_input)
    record.SetProp("generation_method", method)
    record.SetProp("set_tier", set_tier)
    record.SetProp("conf_id", str(conf_id))
    record.SetProp("num_rotatable_bonds", str(rotatable_bonds))
    record.SetProp("num_target_confs", str(num_target_confs))
    record.SetProp("minimization_applied", str(bool(minimization_applied)))
    if conf.HasProp("min_rmsd_to_selected"):
        record.SetDoubleProp("min_rmsd_to_selected", conf.GetDoubleProp("min_rmsd_to_selected"))
    if conf.HasProp("force_field_used"):
        record.SetProp("force_field_used", conf.GetProp("force_field_used"))
    if conf.HasProp("post_min_energy"):
        record.SetDoubleProp("post_min_energy", conf.GetDoubleProp("post_min_energy"))
    if conf.HasProp("fixed_pool_index"):
        record.SetIntProp("fixed_pool_index", conf.GetIntProp("fixed_pool_index"))


def write_molecule_sdf(
    mol: Chem.Mol,
    input_mol: InputMolecule,
    output_path: Path,
    method: str,
    set_tier: str,
    rotatable_bonds: int,
    num_target_confs: int,
    minimization_applied: bool,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if mol.GetNumConformers() == 0:
        output_path.write_text("", encoding="utf-8")
        return
    writer = Chem.SDWriter(str(output_path))
    try:
        for conf_id, conf in enumerate(mol.GetConformers()):
            record = Chem.Mol(mol)
            record.RemoveAllConformers()
            single_conf = Chem.Conformer(conf)
            single_conf.SetId(0)
            record.AddConformer(single_conf, assignId=True)
            set_record_props(
                record=record,
                input_mol=input_mol,
                method=method,
                conf_id=conf_id,
                set_tier=set_tier,
                rotatable_bonds=rotatable_bonds,
                num_target_confs=num_target_confs,
                minimization_applied=minimization_applied,
                conf=conf,
            )
            writer.write(record, confId=0)
    finally:
        writer.close()


def empty_result(
    input_mol: InputMolecule,
    method: str,
    status: str,
    sdf_path: Path,
    walltime_seconds: float = 0.0,
) -> MethodResult:
    if method.endswith("_fixed"):
        tier = "fixed"
    elif method.endswith("_dynamic"):
        tier = "dynamic"
    elif method.endswith("_chembl_count"):
        tier = "chembl_count"
    else:
        tier = "unknown"
    return MethodResult(
        mol_id=input_mol.mol_id,
        input_smiles=input_mol.smiles,
        source_input=input_mol.source_input,
        generation_method=method,
        set_tier=tier,
        num_target_confs=0,
        rotatable_bonds=0,
        generated_candidates=0,
        kept_confs=0,
        finite_rejected=0,
        clash_rejected=0,
        bond_rejected=0,
        stereo_rejected=0,
        rmsd_rejected=0,
        pre_clash_passed=0,
        generation_batches=0,
        minimization_applied=method in MINIMIZED_METHODS,
        minimization_input_confs=0,
        minimization_failed=0,
        minimization_invalid=0,
        minimization_error=0,
        post_min_clash_rejected=0,
        selected_confs=0,
        waste_ratio=math.nan,
        minimization_failed_rate=math.nan,
        minimization_error_rate=math.nan,
        post_min_clash_rejected_rate=math.nan,
        pb_input_confs=0,
        pb_pass_confs=0,
        pb_fail_confs=0,
        pb_pass_rate=math.nan,
        pb_check_fail_counts_json="{}",
        pb_check_fail_rates_json="{}",
        walltime_seconds=walltime_seconds,
        status=status,
        sdf_path=str(sdf_path),
    )


def make_method_result(
    input_mol: InputMolecule,
    method: str,
    set_tier: str,
    num_target_confs: int,
    rotatable_bonds: int,
    stats: FilterStats,
    pre_pb_confs: int,
    pb_input: int,
    pb_pass: int,
    pb_fail: int,
    kept_confs: int,
    status: str,
    sdf_path: Path,
    minimization_applied: bool,
    minimization_input_confs: int = 0,
    minimization_failed: int = 0,
    minimization_invalid: int = 0,
    minimization_error: int = 0,
    post_min_clash_rejected: int = 0,
    pb_check_fail_counts_json: str = "{}",
    pb_check_fail_rates_json: str = "{}",
    walltime_seconds: float = 0.0,
) -> MethodResult:
    pb_rate = (pb_pass / pb_input) if pb_input > 0 else math.nan
    minimization_failed_rate = safe_ratio(minimization_failed, minimization_input_confs)
    minimization_error_rate = safe_ratio(minimization_error, minimization_input_confs)
    post_min_clash_rejected_rate = safe_ratio(post_min_clash_rejected, minimization_input_confs)
    return MethodResult(
        mol_id=input_mol.mol_id,
        input_smiles=input_mol.smiles,
        source_input=input_mol.source_input,
        generation_method=method,
        set_tier=set_tier,
        num_target_confs=num_target_confs,
        rotatable_bonds=rotatable_bonds,
        generated_candidates=stats.generated_candidates,
        kept_confs=kept_confs,
        finite_rejected=stats.finite_rejected,
        clash_rejected=stats.clash_rejected,
        bond_rejected=stats.bond_rejected,
        stereo_rejected=stats.stereo_rejected,
        rmsd_rejected=stats.rmsd_rejected,
        pre_clash_passed=pre_pb_confs,
        generation_batches=stats.generation_batches,
        minimization_applied=minimization_applied,
        minimization_input_confs=minimization_input_confs,
        minimization_failed=minimization_failed,
        minimization_invalid=minimization_invalid,
        minimization_error=minimization_error,
        post_min_clash_rejected=post_min_clash_rejected,
        selected_confs=pre_pb_confs,
        waste_ratio=safe_ratio(stats.generated_candidates, pre_pb_confs),
        minimization_failed_rate=minimization_failed_rate,
        minimization_error_rate=minimization_error_rate,
        post_min_clash_rejected_rate=post_min_clash_rejected_rate,
        pb_input_confs=pb_input,
        pb_pass_confs=pb_pass,
        pb_fail_confs=pb_fail,
        pb_pass_rate=pb_rate,
        pb_check_fail_counts_json=pb_check_fail_counts_json,
        pb_check_fail_rates_json=pb_check_fail_rates_json,
        walltime_seconds=walltime_seconds,
        status=status,
        sdf_path=str(sdf_path),
    )


def finalize_pipeline_pair(
    input_mol: InputMolecule,
    fixed_method: str,
    dynamic_method: str,
    chembl_count_method: str,
    family: str,
    fixed_pre_pb: Chem.Mol,
    dynamic_indices: list[int],
    chembl_count_indices: list[int],
    template: Chem.Mol,
    reference_mol: Chem.Mol,
    paths: dict[str, Path],
    rotatable_bonds: int,
    fixed_target: int,
    dynamic_target: int,
    chembl_count_target: int,
    stats: FilterStats,
    pool_status_value: str,
    minimization_applied: bool,
    minimization_input_confs: int = 0,
    minimization_failed: int = 0,
    minimization_invalid: int = 0,
    minimization_error: int = 0,
    post_min_clash_rejected: int = 0,
    posebusters_workers: int = 1,
    posebusters_energy_threads: int = 1,
    pipeline_start_time: float = 0.0,
) -> list[MethodResult]:
    # The sampling frame is the full 1k pre-PoseBusters pool. Subset indices are
    # chosen before PB so failed conformers do not disappear before sampling.
    dynamic_pre_pb = subset_mol_by_indices(fixed_pre_pb, dynamic_indices, template)
    chembl_count_pre_pb = subset_mol_by_indices(fixed_pre_pb, chembl_count_indices, template)

    fixed_records = [single_conf_mol(fixed_pre_pb, conf.GetId()) for conf in fixed_pre_pb.GetConformers()]
    fixed_pb = posebusters_geometry_passes(
        fixed_records,
        reference_mol,
        max_workers=posebusters_workers,
        energy_num_threads=posebusters_energy_threads,
    )
    fixed_passes = fixed_pb.passes
    if len(fixed_passes) != fixed_pre_pb.GetNumConformers():
        raise RuntimeError(
            f"PoseBusters returned {len(fixed_passes)} result(s) for "
            f"{fixed_pre_pb.GetNumConformers()} conformer(s)."
        )
    dynamic_passes = [fixed_passes[source_idx] for source_idx in dynamic_indices]
    chembl_count_passes = [fixed_passes[source_idx] for source_idx in chembl_count_indices]
    dynamic_check_rows = [fixed_pb.check_rows[source_idx] for source_idx in dynamic_indices]
    chembl_count_check_rows = [fixed_pb.check_rows[source_idx] for source_idx in chembl_count_indices]

    fixed_out = filter_mol_by_passes(fixed_pre_pb, fixed_passes)
    dynamic_out = filter_mol_by_passes(dynamic_pre_pb, dynamic_passes)
    chembl_count_out = filter_mol_by_passes(chembl_count_pre_pb, chembl_count_passes)
    if dynamic_out.GetNumConformers() == 0:
        dynamic_status = "empty_dynamic_subset"
    elif dynamic_target > fixed_pre_pb.GetNumConformers():
        dynamic_status = "dynamic_target_capped"
    else:
        dynamic_status = "ok"
    if chembl_count_out.GetNumConformers() == 0:
        chembl_count_status = "empty_chembl_count_subset"
    elif chembl_count_target > fixed_pre_pb.GetNumConformers():
        chembl_count_status = "chembl_count_target_capped"
    else:
        chembl_count_status = "ok"

    write_molecule_sdf(
        fixed_out,
        input_mol,
        paths[fixed_method],
        fixed_method,
        "fixed",
        rotatable_bonds,
        fixed_target,
        minimization_applied,
    )
    write_indices_sidecar(
        indices_sidecar_path(paths[fixed_method]),
        input_mol.mol_id,
        family,
        "fixed",
        list(range(fixed_pre_pb.GetNumConformers())),
    )
    write_molecule_sdf(
        dynamic_out,
        input_mol,
        paths[dynamic_method],
        dynamic_method,
        "dynamic",
        rotatable_bonds,
        dynamic_target,
        minimization_applied,
    )
    write_indices_sidecar(
        indices_sidecar_path(paths[dynamic_method]),
        input_mol.mol_id,
        family,
        "dynamic",
        dynamic_indices,
    )
    write_molecule_sdf(
        chembl_count_out,
        input_mol,
        paths[chembl_count_method],
        chembl_count_method,
        "chembl_count",
        rotatable_bonds,
        chembl_count_target,
        minimization_applied,
    )
    write_indices_sidecar(
        indices_sidecar_path(paths[chembl_count_method]),
        input_mol.mol_id,
        family,
        "chembl_count",
        chembl_count_indices,
    )

    fixed_pb_pass = sum(fixed_passes)
    dynamic_pb_pass = sum(dynamic_passes)
    chembl_count_pb_pass = sum(chembl_count_passes)
    fixed_fail_counts_json, fixed_fail_rates_json = summarize_posebusters_checks(fixed_pb.check_rows)
    dynamic_fail_counts_json, dynamic_fail_rates_json = summarize_posebusters_checks(dynamic_check_rows)
    chembl_fail_counts_json, chembl_fail_rates_json = summarize_posebusters_checks(
        chembl_count_check_rows
    )
    walltime_seconds = time.perf_counter() - pipeline_start_time if pipeline_start_time else 0.0

    return [
        make_method_result(
            input_mol,
            fixed_method,
            "fixed",
            fixed_target,
            rotatable_bonds,
            stats,
            pre_pb_confs=fixed_pre_pb.GetNumConformers(),
            pb_input=len(fixed_passes),
            pb_pass=fixed_pb_pass,
            pb_fail=len(fixed_passes) - fixed_pb_pass,
            kept_confs=fixed_out.GetNumConformers(),
            status=pool_status_value if fixed_pre_pb.GetNumConformers() == fixed_target else "failed_to_fill_pool",
            sdf_path=paths[fixed_method],
            minimization_applied=minimization_applied,
            minimization_input_confs=minimization_input_confs,
            minimization_failed=minimization_failed,
            minimization_invalid=minimization_invalid,
            minimization_error=minimization_error,
            post_min_clash_rejected=post_min_clash_rejected,
            pb_check_fail_counts_json=fixed_fail_counts_json,
            pb_check_fail_rates_json=fixed_fail_rates_json,
            walltime_seconds=walltime_seconds,
        ),
        make_method_result(
            input_mol,
            dynamic_method,
            "dynamic",
            dynamic_target,
            rotatable_bonds,
            stats,
            pre_pb_confs=dynamic_pre_pb.GetNumConformers(),
            pb_input=len(dynamic_passes),
            pb_pass=dynamic_pb_pass,
            pb_fail=len(dynamic_passes) - dynamic_pb_pass,
            kept_confs=dynamic_out.GetNumConformers(),
            status=dynamic_status,
            sdf_path=paths[dynamic_method],
            minimization_applied=minimization_applied,
            minimization_input_confs=minimization_input_confs,
            minimization_failed=minimization_failed,
            minimization_invalid=minimization_invalid,
            minimization_error=minimization_error,
            post_min_clash_rejected=post_min_clash_rejected,
            pb_check_fail_counts_json=dynamic_fail_counts_json,
            pb_check_fail_rates_json=dynamic_fail_rates_json,
            walltime_seconds=walltime_seconds,
        ),
        make_method_result(
            input_mol,
            chembl_count_method,
            "chembl_count",
            chembl_count_target,
            rotatable_bonds,
            stats,
            pre_pb_confs=chembl_count_pre_pb.GetNumConformers(),
            pb_input=len(chembl_count_passes),
            pb_pass=chembl_count_pb_pass,
            pb_fail=len(chembl_count_passes) - chembl_count_pb_pass,
            kept_confs=chembl_count_out.GetNumConformers(),
            status=chembl_count_status,
            sdf_path=paths[chembl_count_method],
            minimization_applied=minimization_applied,
            minimization_input_confs=minimization_input_confs,
            minimization_failed=minimization_failed,
            minimization_invalid=minimization_invalid,
            minimization_error=minimization_error,
            post_min_clash_rejected=post_min_clash_rejected,
            pb_check_fail_counts_json=chembl_fail_counts_json,
            pb_check_fail_rates_json=chembl_fail_rates_json,
            walltime_seconds=walltime_seconds,
        ),
    ]


def process_rdkit_raw(
    input_mol: InputMolecule,
    base_mol: Chem.Mol,
    reference_mol: Chem.Mol,
    paths: dict[str, Path],
    rotatable_bonds: int,
    dynamic_target: int,
    chembl_count_target: int,
    args: argparse.Namespace,
) -> list[MethodResult]:
    pipeline_start_time = time.perf_counter()
    rdkit_args = clone_args(args, pre_clash_cutoff=False, post_clash_cutoff=False)
    pool, stats, status = accumulate_pre_clash_pool_rdkit(
        base_mol=base_mol,
        mol_id=input_mol.mol_id,
        target=args.fixed_set_size,
        seed=args.seed,
        num_threads=args.num_threads,
        args=rdkit_args,
    )
    if status != "ok":
        walltime_seconds = time.perf_counter() - pipeline_start_time
        return [
            empty_result(input_mol, "rdkit_random_raw_fixed", status, paths["rdkit_random_raw_fixed"], walltime_seconds),
            empty_result(input_mol, "rdkit_random_raw_dynamic", status, paths["rdkit_random_raw_dynamic"], walltime_seconds),
            empty_result(input_mol, "rdkit_random_raw_chembl_count", status, paths["rdkit_random_raw_chembl_count"], walltime_seconds),
        ]
    fixed_raw = build_raw_mol_from_pool(pool, base_mol)
    dynamic_indices = select_dynamic_indices(
        fixed_raw.GetNumConformers(),
        dynamic_target,
        args.seed,
        f"{input_mol.mol_id}:rdkit_random_raw:dynamic_indices",
    )
    chembl_count_indices = select_dynamic_indices(
        fixed_raw.GetNumConformers(),
        chembl_count_target,
        args.seed,
        f"{input_mol.mol_id}:rdkit_random_raw:chembl_count_indices",
    )
    return finalize_pipeline_pair(
        input_mol,
        "rdkit_random_raw_fixed",
        "rdkit_random_raw_dynamic",
        "rdkit_random_raw_chembl_count",
        "rdkit",
        fixed_raw,
        dynamic_indices,
        chembl_count_indices,
        base_mol,
        reference_mol,
        paths,
        rotatable_bonds,
        args.fixed_set_size,
        dynamic_target,
        chembl_count_target,
        stats,
        status,
        minimization_applied=False,
        posebusters_workers=int(args.minimize_workers or args.num_threads),
        posebusters_energy_threads=int(args.num_threads),
        pipeline_start_time=pipeline_start_time,
    )


def accumulate_minimized_pool_rdkit(
    base_mol: Chem.Mol,
    mol_id: str,
    target: int,
    seed: int,
    num_threads: int,
    args: argparse.Namespace,
) -> tuple[list[tuple[Chem.Mol, float, str]], FilterStats, str, int, int, int, int]:
    rdkit_args = clone_args(args, pre_clash_cutoff=False, post_clash_cutoff=False)
    min_args = clone_args(args, pre_clash_cutoff=False, post_clash_cutoff=False)

    final_pool: list[tuple[Chem.Mol, float, str]] = []
    total_generated = 0
    finite_rejected = 0
    pre_clash_rejected = 0
    minimization_input = 0
    min_failed = 0
    min_error = 0
    post_rejected = 0
    generation_batches = 0
    batch_index = 0
    embed_failed = False

    while len(final_pool) < target:
        raw_pool: list[tuple[Chem.Mol, float, str]] = []
        remaining = target - len(final_pool)
        pre_target = max(remaining, RDKIT_MIN_PRE_POOL_SIZE)
        pre_target = min(pre_target, int(args.generation_batch_size))
        while len(raw_pool) < pre_target:
            batch_size = min(int(args.generation_batch_size), pre_target - len(raw_pool))
            survivors, generated, batch_finite, batch_clash = generate_rdkit_batch(
                base_mol=base_mol,
                mol_id=mol_id,
                batch_size=batch_size,
                batch_index=batch_index,
                seed=seed,
                num_threads=num_threads,
                args=rdkit_args,
            )
            batch_index += 1
            generation_batches += 1
            total_generated += generated
            finite_rejected += batch_finite
            pre_clash_rejected += batch_clash
            if generated == 0:
                embed_failed = True
                break
            raw_pool.extend(survivors)

        if embed_failed and not raw_pool:
            break

        minimization_input += len(raw_pool)
        minimized, batch_min_failed, batch_min_error, batch_post_rejected = minimize_candidates_parallel(
            raw_pool,
            min_args,
        )
        min_failed += batch_min_failed
        min_error += batch_min_error
        post_rejected += batch_post_rejected
        for item in minimized:
            final_pool.append(item)
            if len(final_pool) == target:
                break
        del raw_pool, minimized
        gc.collect()

    stats = FilterStats(
        generated_candidates=total_generated,
        finite_rejected=finite_rejected,
        clash_rejected=pre_clash_rejected,
        pre_clash_passed=len(final_pool),
        generation_batches=generation_batches,
    )
    return final_pool, stats, pool_status(target, len(final_pool), embed_failed), minimization_input, min_failed, min_error, post_rejected


def process_rdkit_minimized(
    input_mol: InputMolecule,
    base_mol: Chem.Mol,
    reference_mol: Chem.Mol,
    paths: dict[str, Path],
    rotatable_bonds: int,
    dynamic_target: int,
    chembl_count_target: int,
    args: argparse.Namespace,
) -> list[MethodResult]:
    pipeline_start_time = time.perf_counter()
    minimized, stats, status, minimization_input, min_failed, min_error, post_rejected = accumulate_minimized_pool_rdkit(
        base_mol=base_mol,
        mol_id=input_mol.mol_id,
        target=args.fixed_set_size,
        seed=stable_seed(args.seed, f"{input_mol.mol_id}:rdkit_min"),
        num_threads=args.num_threads,
        args=args,
    )
    if status != "ok":
        walltime_seconds = time.perf_counter() - pipeline_start_time
        return [
            empty_result(input_mol, "rdkit_random_minimized_fixed", status, paths["rdkit_random_minimized_fixed"], walltime_seconds),
            empty_result(input_mol, "rdkit_random_minimized_dynamic", status, paths["rdkit_random_minimized_dynamic"], walltime_seconds),
            empty_result(input_mol, "rdkit_random_minimized_chembl_count", status, paths["rdkit_random_minimized_chembl_count"], walltime_seconds),
        ]

    fixed_min = build_raw_mol_from_pool(minimized, base_mol)
    dynamic_indices = select_dynamic_indices(
        fixed_min.GetNumConformers(),
        dynamic_target,
        args.seed,
        f"{input_mol.mol_id}:rdkit_random_min:dynamic_indices",
    )
    chembl_count_indices = select_dynamic_indices(
        fixed_min.GetNumConformers(),
        chembl_count_target,
        args.seed,
        f"{input_mol.mol_id}:rdkit_random_min:chembl_count_indices",
    )
    return finalize_pipeline_pair(
        input_mol,
        "rdkit_random_minimized_fixed",
        "rdkit_random_minimized_dynamic",
        "rdkit_random_minimized_chembl_count",
        "rdkit",
        fixed_min,
        dynamic_indices,
        chembl_count_indices,
        base_mol,
        reference_mol,
        paths,
        rotatable_bonds,
        args.fixed_set_size,
        dynamic_target,
        chembl_count_target,
        stats,
        "ok",
        minimization_applied=True,
        minimization_input_confs=minimization_input,
        minimization_failed=min_failed,
        minimization_invalid=0,
        minimization_error=min_error,
        post_min_clash_rejected=post_rejected,
        posebusters_workers=int(args.minimize_workers or args.num_threads),
        posebusters_energy_threads=int(args.num_threads),
        pipeline_start_time=pipeline_start_time,
    )


def process_torsion_raw(
    input_mol: InputMolecule,
    torsion_ref: Chem.Mol,
    paths: dict[str, Path],
    rotatable_bonds: int,
    dynamic_target: int,
    chembl_count_target: int,
    args: argparse.Namespace,
    validation_mol: Chem.Mol | None = None,
) -> list[MethodResult]:
    pipeline_start_time = time.perf_counter()
    reference_mol = torsion_ref if validation_mol is None else validation_mol
    torsion_args = clone_args(args, pre_clash_cutoff=0.7)
    pool, stats, status = accumulate_pre_clash_pool_torsion(
        ref_mol=torsion_ref,
        mol_id=input_mol.mol_id,
        target=args.fixed_set_size,
        seed=args.seed,
        args=torsion_args,
        seed_key="torsion_raw",
    )
    if status != "ok":
        walltime_seconds = time.perf_counter() - pipeline_start_time
        return [
            empty_result(input_mol, "torsion_raw_fixed", status, paths["torsion_raw_fixed"], walltime_seconds),
            empty_result(input_mol, "torsion_raw_dynamic", status, paths["torsion_raw_dynamic"], walltime_seconds),
            empty_result(input_mol, "torsion_raw_chembl_count", status, paths["torsion_raw_chembl_count"], walltime_seconds),
        ]
    fixed_raw = build_raw_mol_from_pool(pool, torsion_ref)
    dynamic_indices = select_dynamic_indices(
        fixed_raw.GetNumConformers(),
        dynamic_target,
        args.seed,
        f"{input_mol.mol_id}:torsion_raw:dynamic_indices",
    )
    chembl_count_indices = select_dynamic_indices(
        fixed_raw.GetNumConformers(),
        chembl_count_target,
        args.seed,
        f"{input_mol.mol_id}:torsion_raw:chembl_count_indices",
    )
    return finalize_pipeline_pair(
        input_mol,
        "torsion_raw_fixed",
        "torsion_raw_dynamic",
        "torsion_raw_chembl_count",
        "torsion",
        fixed_raw,
        dynamic_indices,
        chembl_count_indices,
        torsion_ref,
        reference_mol,
        paths,
        rotatable_bonds,
        args.fixed_set_size,
        dynamic_target,
        chembl_count_target,
        stats,
        status,
        minimization_applied=False,
        posebusters_workers=int(args.minimize_workers or args.num_threads),
        posebusters_energy_threads=int(args.num_threads),
        pipeline_start_time=pipeline_start_time,
    )


def process_torsion_minimized(
    input_mol: InputMolecule,
    torsion_ref: Chem.Mol,
    paths: dict[str, Path],
    rotatable_bonds: int,
    dynamic_target: int,
    chembl_count_target: int,
    args: argparse.Namespace,
    validation_mol: Chem.Mol | None = None,
) -> list[MethodResult]:
    pipeline_start_time = time.perf_counter()
    reference_mol = torsion_ref if validation_mol is None else validation_mol
    pool, stats, status, minimization_input, min_failed, min_error, post_clash_rejected = accumulate_torsion_minimized_pool(
        ref_mol=torsion_ref,
        mol_id=input_mol.mol_id,
        target=args.fixed_set_size,
        seed=stable_seed(args.seed, f"{input_mol.mol_id}:torsion_min"),
        args=args,
    )
    if status != "ok":
        walltime_seconds = time.perf_counter() - pipeline_start_time
        return [
            empty_result(input_mol, "torsion_minimized_fixed", status, paths["torsion_minimized_fixed"], walltime_seconds),
            empty_result(input_mol, "torsion_minimized_dynamic", status, paths["torsion_minimized_dynamic"], walltime_seconds),
            empty_result(input_mol, "torsion_minimized_chembl_count", status, paths["torsion_minimized_chembl_count"], walltime_seconds),
        ]
    fixed_min = build_raw_mol_from_pool(pool, torsion_ref)
    dynamic_indices = select_dynamic_indices(
        fixed_min.GetNumConformers(),
        dynamic_target,
        args.seed,
        f"{input_mol.mol_id}:torsion_min:dynamic_indices",
    )
    chembl_count_indices = select_dynamic_indices(
        fixed_min.GetNumConformers(),
        chembl_count_target,
        args.seed,
        f"{input_mol.mol_id}:torsion_min:chembl_count_indices",
    )
    return finalize_pipeline_pair(
        input_mol,
        "torsion_minimized_fixed",
        "torsion_minimized_dynamic",
        "torsion_minimized_chembl_count",
        "torsion",
        fixed_min,
        dynamic_indices,
        chembl_count_indices,
        torsion_ref,
        reference_mol,
        paths,
        rotatable_bonds,
        args.fixed_set_size,
        dynamic_target,
        chembl_count_target,
        stats,
        status,
        minimization_applied=True,
        minimization_input_confs=minimization_input,
        minimization_failed=min_failed,
        minimization_invalid=0,
        minimization_error=min_error,
        post_min_clash_rejected=post_clash_rejected,
        posebusters_workers=int(args.minimize_workers or args.num_threads),
        posebusters_energy_threads=int(args.num_threads),
        pipeline_start_time=pipeline_start_time,
    )


def load_intersection_molecules(
    chembl_map_csv: Path,
    ligand_dir: Path,
    limit: int | None,
    offset: int,
) -> list[InputMolecule]:
    molecules: list[InputMolecule] = []
    raw_rows: list[dict[str, str]] = []
    with chembl_map_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_rows.append(row)
            ligand_id = row["ligand_id"].strip()
            smiles = row.get("chembl3d_isomeric_smiles") or row.get("casf_heavy_isomeric_smiles") or ""
            source_file = row.get("source_file", f"{ligand_id}.mol2").strip()
            raw_conformer_count = str(row.get("conformer_count", "0")).strip()
            molecules.append(
                InputMolecule(
                    mol_id=safe_mol_id(ligand_id),
                    smiles=smiles.strip(),
                    source_input=str(ligand_dir / source_file),
                    chembl3d_group=str(row["chembl3d_group"]).strip(),
                    chembl3d_mol_id=str(row["chembl3d_mol_id"]).strip(),
                    chembl3d_conformer_count=(
                        int(raw_conformer_count) if raw_conformer_count.isdigit() else 0
                    ),
                    chembl3d_sdf_record_index=parse_optional_int(row.get("chembl3d_sdf_record_index")),
                )
            )
    validate_stereo_mapping_csv(raw_rows, path=chembl_map_csv)

    if offset:
        if offset >= len(molecules):
            return []
        molecules = molecules[offset:]
    if limit is not None:
        molecules = molecules[:limit]
    return molecules


def process_molecule(
    input_mol: InputMolecule,
    args: argparse.Namespace,
) -> list[MethodResult]:
    output_dir = Path(args.output_dir)
    paths = {method: output_dir / method / f"{input_mol.mol_id}.sdf" for method in ALL_METHODS}

    if not input_mol.chembl3d_group or not input_mol.chembl3d_mol_id or not input_mol.smiles:
        return [
            empty_result(input_mol, method, "missing_chembl3d_mapping", paths[method])
            for method in ALL_METHODS
        ]

    zarr_root = Path(getattr(args, "chembl3d_zarr_root", DEFAULT_ZARR_ROOT))
    torsion = resolve_torsion_generation(
        input_mol.chembl3d_group,
        input_mol.chembl3d_mol_id,
        Path(args.chembl3d_topology_root),
        Path(input_mol.source_input),
        expected_smiles=input_mol.smiles,
        sdf_record_index=input_mol.chembl3d_sdf_record_index,
        embed_seed=stable_seed(args.seed, f"{input_mol.mol_id}:torsion_seed"),
        zarr_root=zarr_root,
    )
    if torsion is None:
        return [
            empty_result(input_mol, method, "chembl3d_ref_load_failed", paths[method])
            for method in ALL_METHODS
        ]
    validation_mol = torsion.validation_mol
    torsion_seed = torsion.seed_mol

    base_mol = Chem.Mol(validation_mol)
    base_mol.RemoveAllConformers()
    rotatable_bonds = int(Descriptors.NumRotatableBonds(base_mol))
    dynamic_target = dynamic_candidate_count(rotatable_bonds)
    chembl_target = chembl_count_target(input_mol)
    report_dynamic_target_cap(
        input_mol,
        rotatable_bonds,
        dynamic_target,
        args.fixed_set_size,
    )
    report_chembl_count_target_cap(input_mol, chembl_target, args.fixed_set_size)

    results = [
        *process_rdkit_raw(
            input_mol,
            base_mol,
            validation_mol,
            paths,
            rotatable_bonds,
            dynamic_target,
            chembl_target,
            args,
        ),
        *process_rdkit_minimized(
            input_mol,
            base_mol,
            validation_mol,
            paths,
            rotatable_bonds,
            dynamic_target,
            chembl_target,
            args,
        ),
    ]

    if not get_rotatable_torsions(torsion_seed):
        for method in (
            "torsion_raw_fixed",
            "torsion_raw_dynamic",
            "torsion_raw_chembl_count",
            "torsion_minimized_fixed",
            "torsion_minimized_dynamic",
            "torsion_minimized_chembl_count",
        ):
            results.append(empty_result(input_mol, method, "no_rotatable_bonds", paths[method]))
        return results

    results.extend(
        process_torsion_raw(
            input_mol,
            torsion_seed,
            paths,
            rotatable_bonds,
            dynamic_target,
            chembl_target,
            args,
            validation_mol=validation_mol,
        )
    )
    results.extend(
        process_torsion_minimized(
            input_mol,
            torsion_seed,
            paths,
            rotatable_bonds,
            dynamic_target,
            chembl_target,
            args,
            validation_mol=validation_mol,
        )
    )
    return results


def write_manifest(path: Path, rows: list[MethodResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "mol_id",
        "input_smiles",
        "source_input",
        "generation_method",
        "set_tier",
        "num_target_confs",
        "rotatable_bonds",
        "generated_candidates",
        "kept_confs",
        "finite_rejected",
        "clash_rejected",
        "bond_rejected",
        "stereo_rejected",
        "rmsd_rejected",
        "pre_clash_passed",
        "generation_batches",
        "minimization_applied",
        "minimization_input_confs",
        "minimization_failed",
        "minimization_invalid",
        "minimization_error",
        "post_min_clash_rejected",
        "selected_confs",
        "waste_ratio",
        "minimization_failed_rate",
        "minimization_error_rate",
        "post_min_clash_rejected_rate",
        "pb_input_confs",
        "pb_pass_confs",
        "pb_fail_confs",
        "pb_pass_rate",
        "pb_check_fail_counts_json",
        "pb_check_fail_rates_json",
        "walltime_seconds",
        "status",
        "sdf_path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fields})


def active_filter_label(row: MethodResult) -> str:
    filters: list[str] = []
    if row.generation_method.startswith("torsion"):
        filters.append("pre-clash")
    if row.minimization_applied:
        filters.append("minimize")
        if row.generation_method.startswith("torsion"):
            filters.append("post-clash")
    filters.append("PoseBusters")
    return " -> ".join(filters)


def funnel_pool_created(row: MethodResult) -> int:
    if row.minimization_applied:
        return row.minimization_input_confs
    return row.pb_input_confs


def funnel_post_minimized(row: MethodResult) -> str:
    if not row.minimization_applied:
        return "-"
    return str(row.pb_input_confs)


def print_funnel_table(rows: list[MethodResult]) -> None:
    if not rows:
        return
    headers = [
        "method",
        "tier",
        "target",
        "generated",
        "pool_created",
        "post_min",
        "pb_pass",
        "kept",
        "active_filters",
    ]
    table_rows = [
        [
            row.generation_method,
            row.set_tier,
            str(row.num_target_confs),
            str(row.generated_candidates),
            str(funnel_pool_created(row)),
            funnel_post_minimized(row),
            f"{row.pb_pass_confs}/{row.pb_input_confs}",
            str(row.kept_confs),
            active_filter_label(row),
        ]
        for row in rows
    ]
    widths = [
        max(len(headers[idx]), *(len(table_row[idx]) for table_row in table_rows))
        for idx in range(len(headers))
    ]
    print("Generation funnel:")
    print("  " + "  ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)))
    print("  " + "  ".join("-" * width for width in widths))
    for table_row in table_rows:
        print("  " + "  ".join(value.ljust(widths[idx]) for idx, value in enumerate(table_row)))


def load_inputs(args: argparse.Namespace) -> list[InputMolecule]:
    limit = parse_limit_molecules(args.limit_molecules)
    offset = int(args.molecule_offset)
    return load_intersection_molecules(
        Path(args.chembl_map_csv),
        Path(args.ligand_dir),
        limit=limit,
        offset=offset,
    )


def expected_manifest_mol_ids(chembl_map_csv: Path) -> list[str]:
    mol_ids: list[str] = []
    with chembl_map_csv.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            mol_ids.append(safe_mol_id(row["ligand_id"].strip()))
    return mol_ids


def merge_manifest_parts(output_dir: Path, chembl_map_csv: Path) -> Path:
    parts_dir = output_dir / "manifest_parts"
    manifest_path = output_dir / "manifest.tsv"
    part_paths = sorted(parts_dir.glob("*.tsv")) if parts_dir.is_dir() else []
    if not part_paths:
        raise FileNotFoundError(f"No manifest parts found under {parts_dir}")

    expected_ids = expected_manifest_mol_ids(chembl_map_csv)
    if len(expected_ids) != len(set(expected_ids)):
        duplicates = sorted({mol_id for mol_id in expected_ids if expected_ids.count(mol_id) > 1})
        raise ValueError(f"Input CSV contains duplicate safe ligand IDs: {duplicates[:10]}")
    expected_set = set(expected_ids)
    found_stems = {path.stem for path in part_paths}
    missing = sorted(expected_set - found_stems)
    unexpected = sorted(found_stems - expected_set)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing {len(missing)} part(s): {', '.join(missing[:10])}")
        if unexpected:
            details.append(f"unexpected {len(unexpected)} part(s): {', '.join(unexpected[:10])}")
        raise ValueError("Manifest parts do not match input CSV: " + "; ".join(details))

    fields = None
    rows: list[dict[str, object]] = []
    rows_by_mol: dict[str, list[dict[str, object]]] = {}
    for part_path in part_paths:
        with part_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fields is None:
                fields = reader.fieldnames
            elif reader.fieldnames != fields:
                raise ValueError(f"Manifest part {part_path} has inconsistent columns")
            part_rows = list(reader)
        mol_ids = {str(row.get("mol_id", "")) for row in part_rows}
        if mol_ids != {part_path.stem}:
            raise ValueError(f"Manifest part {part_path} has mol_id values {sorted(mol_ids)}")
        methods = {str(row.get("generation_method", "")) for row in part_rows}
        if methods != set(ALL_METHODS):
            missing_methods = sorted(set(ALL_METHODS) - methods)
            extra_methods = sorted(methods - set(ALL_METHODS))
            details = []
            if missing_methods:
                details.append(f"missing methods {missing_methods}")
            if extra_methods:
                details.append(f"unexpected methods {extra_methods}")
            raise ValueError(f"Manifest part {part_path} is incomplete: {'; '.join(details)}")
        rows_by_mol[part_path.stem] = part_rows

    if not fields:
        raise ValueError("Manifest parts contain no header row")

    for mol_id in expected_ids:
        rows.extend(rows_by_mol[mol_id])

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return manifest_path


def run(args: argparse.Namespace) -> list[MethodResult]:
    args.output_dir = resolve_generation_dir(Path(args.output_dir))
    if args.merge_manifest:
        manifest_path = merge_manifest_parts(args.output_dir, Path(args.chembl_map_csv))
        print(f"Merged manifest: {manifest_path} ({manifest_path.stat().st_size} bytes)")
        return []

    if args.minimize_workers is None:
        args.minimize_workers = args.num_threads

    molecules = load_inputs(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[MethodResult] = []
    for input_mol in molecules:
        rows.extend(process_molecule(input_mol, args))

    if len(molecules) == 1:
        part_path = args.output_dir / "manifest_parts" / f"{molecules[0].mol_id}.tsv"
        write_manifest(part_path, rows)
    else:
        write_manifest(args.output_dir / "manifest.tsv", rows)
    return rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate RDKit/torsion conformer sets for ChEMBL3D intersection ligands."
    )
    parser.add_argument(
        "--chembl_map_csv",
        type=Path,
        default=DEFAULT_CHEMBL_MAP_CSV,
        help="CSV mapping ligand_id to chembl3d_group/mol_id.",
    )
    parser.add_argument(
        "--ligand_dir",
        type=Path,
        default=DEFAULT_LIGAND_DIR,
        help="Directory with intersection ligand mol2 files.",
    )
    parser.add_argument(
        "--chembl3d_topology_root",
        type=Path,
        default=DEFAULT_TOPOLOGY_ROOT,
        help="ChEMBL3D topologies/ directory with NNN.sdf files.",
    )
    parser.add_argument(
        "--chembl3d_zarr_root",
        type=Path,
        default=DEFAULT_ZARR_ROOT,
        help="ChEMBL3D zarr_database/ used when a topology SDF has no coordinates.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=DEFAULT_CORE_INTERSECTION_ROOT,
        help="conformer_sets root; SDFs and manifest are written under generation/.",
    )
    parser.add_argument("--limit_molecules", type=str, default="all")
    parser.add_argument(
        "--molecule_offset",
        type=nonnegative_int,
        default=0,
        help="Skip this many intersection CSV rows before applying --limit_molecules.",
    )
    parser.add_argument(
        "--merge_manifest",
        action="store_true",
        help="Merge manifest_parts/*.tsv into manifest.tsv and exit.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--num_threads", type=positive_int, default=1)
    parser.add_argument(
        "--minimize_workers",
        type=positive_int,
        default=None,
        help="Parallel MMFF94 workers inside one molecule process (defaults to --num_threads).",
    )
    parser.add_argument(
        "--candidates_per_rotbond",
        type=positive_int,
        default=17,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--fixed_set_size",
        type=positive_int,
        default=1000,
        help="Exact number of conformers in each fixed-tier pool before PoseBusters.",
    )
    parser.add_argument(
        "--generation_batch_size",
        type=positive_int,
        default=1000,
        help="Trials/embeddings per batch iteration for pool accumulation.",
    )
    parser.add_argument(
        "--torsion_min_pre_pool_size",
        type=positive_int,
        default=1500,
        help="Initial pre-minimize DG @0.7 passers before the first minimize tranche.",
    )
    parser.add_argument(
        "--max_torsion_delta_deg",
        type=nonnegative_float,
        default=120.0,
    )
    parser.add_argument(
        "--perturb_fraction",
        type=nonnegative_float,
        default=1.0,
    )
    parser.add_argument(
        "--pre_clash_cutoff",
        type=float_or_false,
        default=False,
        help="DG-bounds pre-clash cutoff; false disables it (overridden per pipeline).",
    )
    parser.add_argument(
        "--post_clash_cutoff",
        type=float_or_false,
        default=False,
        help="DG-bounds post-clash cutoff; false disables it (overridden per pipeline).",
    )
    parser.add_argument("--ff_variant", choices=["MMFF94", "MMFF94s"], default="MMFF94s")
    parser.add_argument("--max_minimize_iters", type=positive_int, default=500)
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    rows = run(args)
    if args.merge_manifest:
        return

    output_dir = resolve_generation_dir(Path(args.output_dir))
    processed = len({row.mol_id for row in rows})
    print(f"Processed molecules: {processed}")
    for method in ALL_METHODS:
        method_rows = [row for row in rows if row.generation_method == method]
        kept = sum(row.kept_confs for row in method_rows)
        pb_pass = sum(row.pb_pass_confs for row in method_rows)
        pb_input = sum(row.pb_input_confs for row in method_rows)
        if method in MINIMIZED_METHODS:
            print(
                f"{method}: {kept} kept after PB "
                f"({pb_pass}/{pb_input} PB pass pre-filter pool)"
            )
        else:
            trials = [row.generated_candidates for row in method_rows if row.generated_candidates > 0]
            pre_passed = sum(row.pre_clash_passed for row in method_rows)
            if method.endswith("_fixed") and trials:
                median_trials = statistics.median(trials)
                print(
                    f"{method}: {kept} kept after PB "
                    f"(median {median_trials:.0f} trials to fill {args.fixed_set_size}, "
                    f"{pre_passed} pre-PB pool, {pb_pass}/{pb_input} PB pass)"
                )
            else:
                print(
                    f"{method}: {kept} kept after PB "
                    f"({sum(trials) if trials else 0} total trials, "
                    f"{pb_pass}/{pb_input} PB pass)"
                )
    if processed == 1:
        print_funnel_table(rows)
    if len(rows) == 1 or (rows and len({row.mol_id for row in rows}) == 1):
        mol_id = rows[0].mol_id
        print(f"Manifest part: {output_dir / 'manifest_parts' / f'{mol_id}.tsv'}")
    else:
        print(f"Manifest: {output_dir / 'manifest.tsv'}")


if __name__ == "__main__":
    main()
