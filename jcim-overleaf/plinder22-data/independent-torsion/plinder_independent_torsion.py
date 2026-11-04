"""Reference-independent PLINDER torsion baselines.

Starting coordinates come from RDKit ETKDGv3 embeddings of a coordinate-free
molecule built from isomeric SMILES. Experimental conformers are not an input
to embedding or torsion perturbation. They may be compared afterwards for
graph identity only, after their coordinates have been removed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from casf_benchmark.chembl3d.identity import mol_matches_expected_smiles
from casf_benchmark.generation.conformer_sets import (
    build_arg_parser,
    etkdg_params,
    generate_torsion_batch,
    get_dg_bounds,
    get_forcefield,
    get_rotatable_torsions,
    has_finite_coordinates,
    minimize,
    stable_seed,
)

N_ETKDG_STARTS = 5
CONFORMERS_PER_START = 200
RAW_POOL_SIZE = N_ETKDG_STARTS * CONFORMERS_PER_START
MAX_EMBED_ATTEMPTS = 8
MAX_TORSION_BATCHES = 1000
TORSION_DELTA_DEG = 120.0
PERTURB_FRACTION = 1.0
CLASH_CUTOFF = 0.7
TOPOLOGY_SOURCE = "isomeric_smiles"
SEED_SOURCE = "rdkit_etkdg_v3"
RAW_METHOD = "torsion_independent_raw_druglike"
MINIMIZED_METHOD = "torsion_independent_minimized_druglike"


@dataclass
class StartRecord:
    etkdg_index: int
    etkdg_seed: int
    embed_attempts: int
    embed_seed_used: int | None
    attempted: int
    clash_rejected: int
    retained: int

    def as_dict(self) -> dict[str, int | None]:
        return {
            "etkdg_index": self.etkdg_index,
            "etkdg_seed": self.etkdg_seed,
            "embed_attempts": self.embed_attempts,
            "embed_seed_used": self.embed_seed_used,
            "attempted": self.attempted,
            "clash_rejected": self.clash_rejected,
            "retained": self.retained,
        }


@dataclass
class MoleculeGeneration:
    molecule_id: str
    input_smiles: str
    canonical_isomeric_smiles: str
    status: str
    global_seed: int
    etkdg_seeds: list[int]
    starts: list[StartRecord] = field(default_factory=list)
    raw_mols: list[Chem.Mol] = field(default_factory=list)
    minimized_mols: list[Chem.Mol] = field(default_factory=list)
    mmff_statuses: list[dict[str, object]] = field(default_factory=list)
    reference_topology_matches: bool | None = None
    experimental_coordinates_used: bool = False
    topology_source: str = TOPOLOGY_SOURCE
    seed_source: str = SEED_SOURCE

    @property
    def ok(self) -> bool:
        return self.status == "ok"


def canonical_isomeric_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Unparseable SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def topology_from_smiles(smiles: str) -> Chem.Mol:
    """Hydrogenated molecular graph with the requested stereochemistry and no coordinates."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Unparseable SMILES: {smiles!r}")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if mol.GetNumConformers():
        mol.RemoveAllConformers()
    hydrogenated = Chem.AddHs(mol)
    if hydrogenated.GetNumConformers():
        raise RuntimeError("SMILES topology acquired coordinates while adding hydrogens")
    return hydrogenated


def require_coordinate_free(mol: Chem.Mol) -> None:
    if mol.GetNumConformers():
        raise ValueError("ETKDG embedding received a molecule that already has coordinates")


def etkdg_seed(global_seed: int, identity: str, etkdg_index: int) -> int:
    return stable_seed(int(global_seed), f"{identity}:etkdg:{int(etkdg_index)}")


def etkdg_attempt_seed(global_seed: int, identity: str, etkdg_index: int, attempt: int) -> int:
    if attempt == 0:
        return etkdg_seed(global_seed, identity, etkdg_index)
    return stable_seed(
        int(global_seed),
        f"{identity}:etkdg:{int(etkdg_index)}:attempt:{int(attempt)}",
    )


def torsion_args(per_start: int) -> argparse.Namespace:
    args = build_arg_parser().parse_args([])
    args.max_torsion_delta_deg = TORSION_DELTA_DEG
    args.perturb_fraction = PERTURB_FRACTION
    args.pre_clash_cutoff = CLASH_CUTOFF
    args.generation_batch_size = max(1, int(per_start))
    args.minimize_workers = 1
    return args


def embed_etkdg_start(
    topology: Chem.Mol,
    attempt_seeds: list[int],
    *,
    expected_smiles: str,
) -> tuple[Chem.Mol | None, int, int | None]:
    """Embed one conformer. `topology` must not contain coordinates.

    Failed embeddings are retried with independent random-coordinate seeds.
    There is no coordinate fallback.
    """
    require_coordinate_free(topology)
    if not attempt_seeds:
        return None, 0, None
    for attempt, seed in enumerate(attempt_seeds, start=1):
        candidate = Chem.Mol(topology)
        candidate.RemoveAllConformers()
        params = etkdg_params(int(seed), num_threads=1)
        conf_id = AllChem.EmbedMolecule(candidate, params)
        if conf_id < 0 or candidate.GetNumConformers() == 0:
            continue
        if not has_finite_coordinates(candidate, 0):
            continue
        heavy = Chem.RemoveHs(candidate)
        if not mol_matches_expected_smiles(heavy, expected_smiles, from_3d=True):
            continue
        return candidate, attempt, int(seed)
    return None, len(attempt_seeds), None


def conformers_from_embedded_start(
    embedded: Chem.Mol,
    *,
    identity: str,
    global_seed: int,
    etkdg_index: int,
    per_start: int,
    max_batches: int = MAX_TORSION_BATCHES,
) -> tuple[list[Chem.Mol], int, int]:
    """Perturb every eligible torsion of one ETKDG conformer until `per_start` pass the clash filter."""
    if embedded.GetNumConformers() != 1:
        raise ValueError("torsion generation requires exactly one embedded conformer")
    torsions = get_rotatable_torsions(embedded)
    if not torsions:
        return [], 0, 0
    clash = get_dg_bounds(embedded)
    args = torsion_args(per_start)
    pool: list[Chem.Mol] = []
    attempted = 0
    rejected = 0
    batch_index = 0
    while len(pool) < per_start and batch_index < max_batches:
        remaining = per_start - len(pool)
        survivors, generated, batch_rejected = generate_torsion_batch(
            ref_mol=embedded,
            torsions=torsions,
            mol_id=identity,
            batch_size=remaining,
            batch_index=batch_index,
            seed=int(global_seed),
            args=args,
            seed_key=f"torsion_independent:etkdg:{int(etkdg_index)}",
            dg_clash_data=clash,
        )
        batch_index += 1
        attempted += generated
        rejected += batch_rejected
        pool.extend(Chem.RemoveHs(mol) for mol, _energy, _ff in survivors)
        if generated == 0:
            break
    return pool[:per_start], attempted, rejected


def _chiral_dihedrals(mol: Chem.Mol) -> list[tuple[int, int, int, int, float]]:
    conformer = mol.GetConformer()
    dihedrals = []
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue
        neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        if len(neighbors) < 3:
            continue
        first, second, third = neighbors[:3]
        angle = rdMolTransforms.GetDihedralDeg(conformer, first, second, atom.GetIdx(), third)
        dihedrals.append((first, second, atom.GetIdx(), third, float(angle)))
    return dihedrals


def _minimize_holding_chirality(mol: Chem.Mol, *, max_iters: int, ff_variant: str) -> tuple[int, float, str]:
    """MMFF94s from the current coordinates, with each chiral improper held near its start."""
    forcefield, ff_name = get_forcefield(mol, conf_id=0, ff_variant=ff_variant)
    for first, second, center, third, angle in _chiral_dihedrals(mol):
        forcefield.MMFFAddTorsionConstraint(first, second, center, third, False, angle - 30.0, angle + 30.0, 1000.0)
    status = forcefield.Minimize(maxIts=int(max_iters))
    return int(status), float(forcefield.CalcEnergy()), ff_name


def _mmff_ready(mol: Chem.Mol) -> bool:
    try:
        return bool(AllChem.MMFFHasAllMoleculeParams(Chem.AddHs(mol)))
    except Exception:
        return False


def minimize_raw_pool(
    raw_mols: list[Chem.Mol],
    *,
    max_iters: int,
    ff_variant: str,
    expected_smiles: str = "",
) -> tuple[list[Chem.Mol], list[dict[str, object]], str]:
    """Minimize each raw conformer in order and keep that order.

    Status 0 (converged) and 1 (iteration limit) are retained, matching the
    existing PLINDER MMFF94s policy. Any other outcome leaves the pool unpaired
    and is reported as a failure rather than a shifted list.
    """
    statuses: list[dict[str, object]] = []
    minimized: list[Chem.Mol] = []
    if raw_mols and not _mmff_ready(raw_mols[0]):
        for index, raw in enumerate(raw_mols):
            statuses.append(
                {
                    "raw_conformer_index": index,
                    "etkdg_index": int(raw.GetIntProp("etkdg_index")) if raw.HasProp("etkdg_index") else None,
                    "mmff_status": None,
                    "outcome": "missing_mmff_params",
                }
            )
        return [], statuses, "missing_mmff_params"

    failed = False
    for index, raw in enumerate(raw_mols):
        etkdg_index = int(raw.GetIntProp("etkdg_index")) if raw.HasProp("etkdg_index") else None
        record: dict[str, object] = {
            "raw_conformer_index": index,
            "etkdg_index": etkdg_index,
            "mmff_status": None,
            "outcome": "error",
        }
        try:
            hydrated = Chem.AddHs(raw, addCoords=True)
            starting = Chem.Mol(hydrated)
            status, energy, ff_name = minimize(
                hydrated,
                max_iters=int(max_iters),
                ff_variant=ff_variant,
            )
            heavy = Chem.RemoveHs(hydrated)
            stereo_ok = not expected_smiles or mol_matches_expected_smiles(heavy, expected_smiles, from_3d=True)
            if status in (0, 1) and not stereo_ok:
                record["unconstrained_mmff_status"] = int(status)
                record["unconstrained_stereo"] = "inverted"
                status, energy, ff_name = _minimize_holding_chirality(
                    starting,
                    max_iters=int(max_iters),
                    ff_variant=ff_variant,
                )
                heavy = Chem.RemoveHs(starting)
                stereo_ok = mol_matches_expected_smiles(heavy, expected_smiles, from_3d=True)
                record["chiral_constraint_retry"] = True
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
            statuses.append(record)
            failed = True
            continue
        record["mmff_status"] = int(status)
        record["energy"] = float(energy)
        record["force_field"] = ff_name
        if status not in (0, 1) or not stereo_ok:
            record["outcome"] = "stereo_inverted" if status in (0, 1) else "rejected"
            statuses.append(record)
            failed = True
            continue
        record["outcome"] = "retained_chiral_constraint" if record.get("chiral_constraint_retry") else "retained"
        heavy.SetIntProp("raw_conformer_index", index)
        if etkdg_index is not None:
            heavy.SetIntProp("etkdg_index", etkdg_index)
        heavy.SetIntProp("mmff_status", int(status))
        minimized.append(heavy)
        statuses.append(record)
    if failed:
        return [], statuses, "minimization_failed"
    return minimized, statuses, "ok"


def _tag_raw(mol: Chem.Mol, *, raw_index: int, etkdg_index: int) -> Chem.Mol:
    tagged = Chem.Mol(mol)
    tagged.SetIntProp("raw_conformer_index", int(raw_index))
    tagged.SetIntProp("etkdg_index", int(etkdg_index))
    return tagged


def identity_holds(mols: list[Chem.Mol], expected_smiles: str) -> bool:
    return all(mol_matches_expected_smiles(mol, expected_smiles, from_3d=True) for mol in mols)


def reference_topology_matches(reference_mol: Chem.Mol, expected_smiles: str) -> bool:
    """Compare the reference graph after discarding its coordinates."""
    bare = Chem.Mol(reference_mol)
    bare.RemoveAllConformers()
    return bool(mol_matches_expected_smiles(bare, expected_smiles, from_3d=False))


def generate_independent_pools(
    smiles: str,
    *,
    global_seed: int,
    molecule_id: str,
    n_starts: int = N_ETKDG_STARTS,
    per_start: int = CONFORMERS_PER_START,
    max_embed_attempts: int = MAX_EMBED_ATTEMPTS,
    minimize_pool: bool = True,
    max_minimize_iters: int = 500,
    ff_variant: str = "MMFF94s",
) -> MoleculeGeneration:
    """Build raw and minimized torsion pools from SMILES alone."""
    try:
        identity = canonical_isomeric_smiles(smiles)
        topology = topology_from_smiles(identity)
    except ValueError:
        return MoleculeGeneration(
            molecule_id=molecule_id,
            input_smiles=smiles,
            canonical_isomeric_smiles="",
            status="smiles_parse_failed",
            global_seed=int(global_seed),
            etkdg_seeds=[],
        )

    seeds = [etkdg_seed(global_seed, identity, index) for index in range(n_starts)]
    result = MoleculeGeneration(
        molecule_id=molecule_id,
        input_smiles=smiles,
        canonical_isomeric_smiles=identity,
        status="ok",
        global_seed=int(global_seed),
        etkdg_seeds=seeds,
    )
    if not get_rotatable_torsions(topology):
        result.status = "no_rotatable_bonds"
        return result

    raw: list[Chem.Mol] = []
    for index, primary_seed in enumerate(seeds):
        attempt_seeds = [
            etkdg_attempt_seed(global_seed, identity, index, attempt)
            for attempt in range(max_embed_attempts)
        ]
        embedded, attempts_used, seed_used = embed_etkdg_start(
            topology,
            attempt_seeds,
            expected_smiles=identity,
        )
        if embedded is None:
            result.starts.append(
                StartRecord(
                    etkdg_index=index,
                    etkdg_seed=primary_seed,
                    embed_attempts=attempts_used,
                    embed_seed_used=None,
                    attempted=0,
                    clash_rejected=0,
                    retained=0,
                )
            )
            result.status = "embedding_failed"
            result.raw_mols = []
            return result
        conformers, attempted, rejected = conformers_from_embedded_start(
            embedded,
            identity=identity,
            global_seed=int(global_seed),
            etkdg_index=index,
            per_start=per_start,
        )
        result.starts.append(
            StartRecord(
                etkdg_index=index,
                etkdg_seed=primary_seed,
                embed_attempts=attempts_used,
                embed_seed_used=seed_used,
                attempted=attempted,
                clash_rejected=rejected,
                retained=len(conformers),
            )
        )
        if len(conformers) != per_start:
            result.status = "failed_to_fill_pool"
            result.raw_mols = []
            return result
        base = len(raw)
        raw.extend(
            _tag_raw(mol, raw_index=base + offset, etkdg_index=index)
            for offset, mol in enumerate(conformers)
        )

    if not identity_holds(raw, identity):
        result.status = "stereo_mismatch"
        result.raw_mols = []
        return result
    result.raw_mols = raw
    if not minimize_pool:
        return result
    minimized, statuses, min_status = minimize_raw_pool(
        raw,
        max_iters=max_minimize_iters,
        ff_variant=ff_variant,
        expected_smiles=identity,
    )
    result.mmff_statuses = statuses
    if min_status != "ok" or len(minimized) != len(raw):
        result.status = min_status
        result.minimized_mols = []
        return result
    if not identity_holds(minimized, identity):
        result.status = "stereo_mismatch"
        result.minimized_mols = []
        return result
    result.minimized_mols = minimized
    return result


def generate_from_reference_row(
    smiles: str,
    reference_mol: Chem.Mol | None,
    *,
    global_seed: int,
    molecule_id: str,
    n_starts: int = N_ETKDG_STARTS,
    per_start: int = CONFORMERS_PER_START,
    max_embed_attempts: int = MAX_EMBED_ATTEMPTS,
    minimize_pool: bool = True,
    max_minimize_iters: int = 500,
    ff_variant: str = "MMFF94s",
) -> MoleculeGeneration:
    """Generate from SMILES, then compare the reference graph with its coordinates removed."""
    generated = generate_independent_pools(
        smiles,
        global_seed=global_seed,
        molecule_id=molecule_id,
        n_starts=n_starts,
        per_start=per_start,
        max_embed_attempts=max_embed_attempts,
        minimize_pool=minimize_pool,
        max_minimize_iters=max_minimize_iters,
        ff_variant=ff_variant,
    )
    if reference_mol is not None and generated.canonical_isomeric_smiles:
        generated.reference_topology_matches = reference_topology_matches(
            reference_mol,
            generated.canonical_isomeric_smiles,
        )
    generated.experimental_coordinates_used = False
    return generated


def mmff_status_counts(statuses: list[dict[str, object]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for record in statuses:
        value = record.get("mmff_status")
        key = "none" if value is None else str(value)
        counts[key] += 1
    return dict(sorted(counts.items()))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def provenance_document(
    generated: MoleculeGeneration,
    *,
    output_hashes: dict[str, str] | None = None,
) -> dict[str, object]:
    return {
        "molecule_id": generated.molecule_id,
        "input_smiles": generated.input_smiles,
        "canonical_isomeric_smiles": generated.canonical_isomeric_smiles,
        "topology_source": generated.topology_source,
        "seed_source": generated.seed_source,
        "global_seed": generated.global_seed,
        "etkdg_seeds": list(generated.etkdg_seeds),
        "starts": [start.as_dict() for start in generated.starts],
        "raw_conformer_count": len(generated.raw_mols),
        "minimized_conformer_count": len(generated.minimized_mols),
        "mmff_status_counts": mmff_status_counts(generated.mmff_statuses),
        "mmff_statuses": generated.mmff_statuses,
        "experimental_coordinates_used": False,
        "reference_topology_matches": generated.reference_topology_matches,
        "status": generated.status,
        "output_hashes": output_hashes or {},
    }


def _write_sdf(mols: list[Chem.Mol], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    try:
        for mol in mols:
            writer.write(mol)
    finally:
        writer.close()


def write_molecule_provenance(
    generated: MoleculeGeneration,
    provenance_dir: Path,
    *,
    output_hashes: dict[str, str] | None = None,
) -> Path:
    provenance_dir.mkdir(parents=True, exist_ok=True)
    slug = _slug(generated.molecule_id)
    hashes = dict(output_hashes or {})
    if generated.raw_mols:
        raw_sdf = provenance_dir / f"{slug}_raw.sdf"
        _write_sdf(generated.raw_mols, raw_sdf)
        hashes["raw_sdf"] = sha256_file(raw_sdf)
    if generated.minimized_mols:
        minimized_sdf = provenance_dir / f"{slug}_minimized.sdf"
        _write_sdf(generated.minimized_mols, minimized_sdf)
        hashes["minimized_sdf"] = sha256_file(minimized_sdf)
    path = provenance_dir / f"{slug}.json"
    path.write_text(
        json.dumps(provenance_document(generated, output_hashes=hashes), indent=2) + "\n",
        encoding="utf-8",
    )
    return path


def annotate_pickle_hashes(provenance_dir: Path, raw_pickle: Path, minimized_pickle: Path) -> None:
    hashes = {
        "raw_pickle": sha256_file(raw_pickle) if raw_pickle.is_file() else "",
        "minimized_pickle": sha256_file(minimized_pickle) if minimized_pickle.is_file() else "",
    }
    for path in sorted(provenance_dir.glob("*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        output_hashes = dict(document.get("output_hashes") or {})
        output_hashes.update(hashes)
        document["output_hashes"] = output_hashes
        document["experimental_coordinates_used"] = False
        path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def _slug(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "._-" else "_" for character in value.strip())
    return cleaned.strip("._") or "molecule"
