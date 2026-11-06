#!/usr/bin/env python3
"""Reproduce DODock's ligand-similarity filter on the CoordToken corpus.

The DODock rule is a 2048-bit, radius-2 Morgan fingerprint (ECFP4) and
Tanimoto similarity. A candidate training molecule is excluded when its
similarity to any benchmark molecule is strictly greater than 0.4.

The expensive stages are deliberately shardable for Slurm:

  prepare-queries  build-db  search-db  materialize-hits  attribute-sources
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import multiprocessing as mp
import re
import sys
from collections import Counter, defaultdict
from functools import partial
from pathlib import Path

import numpy as np


TRAINING_ROOT = Path("/mnt/weka/fgeikyan/fsq/3dbigdata_fsq_4096")
MERGED_TRAIN = Path("/mnt/weka/vtarasov/3DbigData/merged_train.csv")
FP_PARAMS = {
    "radius": 2,
    "fpSize": 2048,
    "countSimulation": False,
    "includeChirality": False,
    "useBondTypes": True,
    "onlyNonzeroInvariants": False,
    "includeRingMembership": True,
    "includeRedundantEnvironments": False,
}
THRESHOLD = 0.4
COORDINATE_RE = re.compile(r"<[^>]+>")
COORDINATE_BYTES_RE = re.compile(br"<[^>]+>")
SOURCE_NAMES = (
    "OMol25_small_mols",
    "OMol25_bio_mols",
    "CrossDocked2020",
    "BindingNet-High",
    "BindingNet-Low",
    "BindingNet-Mid",
    "Kinodata-3D",
    "BindingMoad",
    "DAVIS-3D",
    "KIBA-3D",
    "pubchem3d",
    "chembl3d",
    "nablaDFT",
    "HiQBind",
    "Plinder",
    "SPINDR",
    "SAIR",
    "geom",
    "zinc",
)
SOURCE_NAME_SET = set(SOURCE_NAMES)
DATASET_BITS = {
    "posebench": 1,
    "platinum": 2,
    "casf": 4,
    "plinder23": 8,
}
_SOURCE_HITS_BY_SMILES: dict[str, list[tuple[int, int]]] = {}
_SOURCE_UNATTRIBUTED_BY_NONISO: dict[str, list[tuple[int, int, str]]] = {}
_TRAINING_IDENTITY_HASHES: np.ndarray | None = None


def identity_digest(smiles: str) -> bytes:
    return hashlib.blake2b(smiles.encode("utf-8"), digest_size=16).digest()


def identity_is_in_training(digest: bytes) -> bool:
    if _TRAINING_IDENTITY_HASHES is None:
        return True
    key = np.asarray([digest], dtype="V16")
    position = int(np.searchsorted(_TRAINING_IDENTITY_HASHES, key)[0])
    return position < len(_TRAINING_IDENTITY_HASHES) and bool(
        _TRAINING_IDENTITY_HASHES[position] == key[0]
    )


def dataset_key(label: str) -> str:
    normalized = label.lower().replace("-", "").replace("_", "")
    if normalized.startswith("posebench"):
        return "posebench"
    if normalized.startswith("platinum"):
        return "platinum"
    if normalized.startswith("casf"):
        return "casf"
    if normalized.startswith("plinder23"):
        return "plinder23"
    raise ValueError(f"unrecognized benchmark dataset label: {label!r}")


def _canonical_source_chunk(rows: list[tuple[str, str]], collect_identities: bool = False):
    """Canonicalize source rows exactly as the original grouping pipeline did."""
    from rdkit import Chem

    molgen_root = Path("/mnt/weka/fgeikyan/3dmolgen_bigdata")
    sys.path.insert(0, str(molgen_root / "src"))
    from molgen3D.data_processing.smiles_encoder_decoder import (
        decode_cartesian_v2,
        strip_smiles,
    )

    hit_counts = Counter()
    identity_hashes: dict[str, bytearray] = defaultdict(bytearray)
    failures = 0
    for source, enriched_text in rows:
        try:
            mol = decode_cartesian_v2(enriched_text)
            try:
                smiles = Chem.MolToSmiles(
                    Chem.RemoveHs(mol), canonical=True, isomericSmiles=True
                )
            except Exception:  # noqa: BLE001
                stripped = strip_smiles(enriched_text)
                parsed = Chem.MolFromSmiles(stripped)
                if parsed is None:
                    raise ValueError("failed to canonicalize source graph")
                smiles = Chem.MolToSmiles(
                    parsed, canonical=True, isomericSmiles=True
                )
            direct_hits = _SOURCE_HITS_BY_SMILES.get(smiles, ())
            for shard, line_id in direct_hits:
                hit_counts[(source, shard, line_id)] += 1
            if _SOURCE_UNATTRIBUTED_BY_NONISO:
                noniso_mol = Chem.Mol(mol)
                Chem.RemoveStereochemistry(noniso_mol)
                noniso = Chem.MolToSmiles(
                    Chem.RemoveHs(noniso_mol), canonical=True, isomericSmiles=False
                )
                candidates = _SOURCE_UNATTRIBUTED_BY_NONISO.get(noniso, ())
                if candidates:
                    stereo_key = Chem.MolToInchiKey(Chem.RemoveHs(mol))
                    for shard, line_id, target_key in candidates:
                        if stereo_key == target_key:
                            hit_counts[(source, shard, line_id)] += 1
            if collect_identities:
                digest = identity_digest(smiles)
                if identity_is_in_training(digest):
                    identity_hashes[source].extend(digest)
        except Exception:  # noqa: BLE001
            failures += 1
    return hit_counts, failures, {key: bytes(value) for key, value in identity_hashes.items()}


def mapping_path(shard: int) -> Path:
    return TRAINING_ROOT / f"shard_{shard:04d}/train_geom_to_isomeric_smiles.jsonl"


def prepare_queries(
    output: Path,
    audit_module_dir: Path,
    extra_query_files: list[Path],
    casf_core_mapping: Path | None,
    plinder_shortlist: Path | None,
    include_standard_benchmarks: bool,
) -> None:
    sys.path.insert(0, str(audit_module_dir))
    from smiles_contamination_audit import platinum_records, posebench_records

    records = (
        posebench_records() + platinum_records()
        if include_standard_benchmarks
        else []
    )
    for path in extra_query_files:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                dataset_key(row["dataset"])
                records.append(
                    {
                        "dataset": row["dataset"],
                        "benchmark_id": row["benchmark_id"],
                        "canonical_isomeric_smiles": row["canonical_isomeric_smiles"],
                        "source_path": row.get("source_path", str(path)),
                    }
                )
    if casf_core_mapping is not None:
        with casf_core_mapping.open(newline="") as handle:
            for row in csv.DictReader(handle):
                smiles = row.get("casf_heavy_isomeric_smiles") or row.get(
                    "canonical_smiles"
                )
                if not smiles:
                    raise ValueError(
                        "CASF input must contain either "
                        "casf_heavy_isomeric_smiles or canonical_smiles"
                    )
                records.append(
                    {
                        "dataset": "casf",
                        "benchmark_id": row["ligand_id"],
                        "canonical_isomeric_smiles": smiles,
                        "structure_path": (
                            str(
                                casf_core_mapping.parent
                                / "ligands"
                                / row["source_file"]
                            )
                            if row.get("canonical_smiles") and row.get("source_file")
                            else ""
                        ),
                        "source_path": str(casf_core_mapping),
                    }
                )
    if plinder_shortlist is not None:
        with plinder_shortlist.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                records.append(
                    {
                        "dataset": "plinder23",
                        "benchmark_id": row.get("ccd") or row["name"],
                        "canonical_isomeric_smiles": row["smiles"],
                        "source_path": str(plinder_shortlist),
                    }
                )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "query_id",
                "dataset",
                "benchmark_id",
                "canonical_isomeric_smiles",
                "structure_path",
                "source_path",
            ),
            extrasaction="ignore",
        )
        writer.writeheader()
        for query_id, record in enumerate(records):
            writer.writerow({"query_id": query_id, **record})
    unique = len({row["canonical_isomeric_smiles"] for row in records})
    print(f"wrote {len(records):,} benchmark rows ({unique:,} unique SMILES) to {output}")


def iter_mapping_smiles(path: Path, limit: int | None = None):
    with path.open() as handle:
        for line_id, line in enumerate(handle, start=1):
            if limit is not None and line_id > limit:
                break
            row = json.loads(line)
            yield row["isomeric_smiles"], line_id


def build_training_identity_hashes(shard: int, output: Path) -> None:
    digests = bytearray()
    for smiles, line_id in iter_mapping_smiles(mapping_path(shard)):
        digests.extend(identity_digest(smiles))
        if line_id % 1_000_000 == 0:
            print(f"shard {shard}: hashed {line_id:,} training molecules", flush=True)
    values = np.frombuffer(digests, dtype="V16")
    unique_values = np.unique(values)
    output.parent.mkdir(parents=True, exist_ok=True)
    unique_values.tofile(output)
    print(f"shard {shard}: wrote {len(unique_values):,} identity hashes to {output}")


def merge_training_identity_hashes(inputs: list[Path], output: Path) -> None:
    arrays = [np.fromfile(path, dtype="V16") for path in inputs]
    values = np.unique(np.concatenate(arrays))
    output.parent.mkdir(parents=True, exist_ok=True)
    values.tofile(output)
    print(f"wrote {len(values):,} retained training identities to {output}")


def build_db(shard: int, output: Path, limit: int | None) -> None:
    from FPSim2.io.backends.pytables import create_db_file

    input_path = mapping_path(shard)
    output.parent.mkdir(parents=True, exist_ok=True)
    create_db_file(
        iter_mapping_smiles(input_path, limit),
        str(output),
        mol_format="smiles",
        fp_type="Morgan",
        fp_params=dict(FP_PARAMS),
        sort_by_popcnt=True,
        full_sanitization=True,
    )
    print(f"built {output} from {input_path}")


def load_queries(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def search_db(
    shard: int,
    database: Path,
    queries_path: Path,
    output: Path,
    workers: int,
) -> None:
    from FPSim2 import FPSim2Engine

    queries = load_queries(queries_path)
    engine = FPSim2Engine(str(database), in_memory_fps=True)
    # The in-memory PyTables backend exposes a dense uint64 matrix whose first
    # column is the original integer molecule ID.
    row_count = max(int(engine.fps[:, 0].max()), 0)
    max_similarity = np.zeros(row_count + 1, dtype=np.float32)
    dataset_flags = np.zeros(row_count + 1, dtype=np.uint8)
    closest_query = np.full(row_count + 1, -1, dtype=np.int32)

    def one_query(index: int):
        row = queries[index]
        result = engine.similarity(
            row["canonical_isomeric_smiles"], THRESHOLD, n_workers=1
        )
        # FPSim2 includes equality at the supplied threshold. DODock says
        # "exceeded 0.4", so enforce the paper's strict inequality here.
        return index, result[result["coeff"] > THRESHOLD]

    pending: set[concurrent.futures.Future] = set()
    next_index = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        while next_index < len(queries) or pending:
            while next_index < len(queries) and len(pending) < workers * 2:
                pending.add(executor.submit(one_query, next_index))
                next_index += 1
            done, pending = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                query_index, result = future.result()
                if not len(result):
                    continue
                ids = result["mol_id"].astype(np.int64, copy=False)
                coeffs = result["coeff"].astype(np.float32, copy=False)
                flag = DATASET_BITS[dataset_key(queries[query_index]["dataset"])]
                dataset_flags[ids] |= flag
                better = coeffs > max_similarity[ids]
                if np.any(better):
                    better_ids = ids[better]
                    max_similarity[better_ids] = coeffs[better]
                    closest_query[better_ids] = query_index
                if (query_index + 1) % 100 == 0:
                    print(
                        f"shard {shard}: completed query {query_index + 1}/{len(queries)}",
                        flush=True,
                    )

    hit_ids = np.flatnonzero(dataset_flags)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        shard=np.int16(shard),
        line_id=hit_ids.astype(np.uint32),
        max_similarity=max_similarity[hit_ids],
        dataset_flags=dataset_flags[hit_ids],
        closest_query=closest_query[hit_ids],
    )
    print(f"shard {shard}: wrote {len(hit_ids):,} removed molecule IDs to {output}")


def search_db_numkong(
    shard: int,
    database: Path,
    queries_path: Path,
    output: Path,
    workers: int,
    block_size: int,
) -> None:
    """Exact blocked MxN Jaccard search using NumKong's SIMD kernel."""
    import numkong as nk
    from FPSim2 import FPSim2Engine

    queries = load_queries(queries_path)
    engine = FPSim2Engine(str(database), in_memory_fps=True)
    if engine.fp_type != "Morgan" or engine.fp_params != FP_PARAMS:
        raise RuntimeError(
            f"unexpected fingerprint configuration: {engine.fp_type} {engine.fp_params}"
        )

    # FPSim2 stores [molecule ID, 32 x uint64 fingerprint words, popcount].
    # Copy the 32 words so every row is contiguous, then reinterpret as 256
    # packed bytes. Bit order does not affect Jaccard provided the same layout
    # is used for training and query fingerprints.
    ids = engine.fps[:, 0].astype(np.uint32, copy=True)
    train_packed = engine.fps[:, 1:-1].copy().view(np.uint8).reshape(len(ids), 256)
    from rdkit import Chem

    def query_fingerprint(row):
        structure_path = row.get("structure_path")
        if structure_path:
            mol = Chem.MolFromMol2File(structure_path, sanitize=True, removeHs=True)
            if mol is None:
                raise ValueError(f"failed to parse query structure: {structure_path}")
            fingerprint = engine.load_query(mol)
        else:
            fingerprint = engine.load_query(row["canonical_isomeric_smiles"])
        return fingerprint[1:-1].copy().view(np.uint8)

    query_packed = np.vstack([query_fingerprint(row) for row in queries])
    packed_queries = nk.hammings_pack(query_packed, dtype="uint1")
    query_indices_by_dataset = {
        key: np.asarray(
            [i for i, row in enumerate(queries) if dataset_key(row["dataset"]) == key],
            dtype=np.int32,
        )
        for key in DATASET_BITS
    }
    distance_cutoff = np.float32(1.0 - THRESHOLD)

    def one_block(start: int):
        end = min(start + block_size, len(ids))
        distances_tensor = nk.zeros((end - start, len(queries)), dtype="float32")
        nk.jaccards_packed(train_packed[start:end], packed_queries, out=distances_tensor)
        distances = np.asarray(distances_tensor)
        qualifies = distances < distance_cutoff
        flags = np.zeros(end - start, dtype=np.uint8)
        for key, bit in DATASET_BITS.items():
            indices = query_indices_by_dataset[key]
            if len(indices):
                flags |= np.any(qualifies[:, indices], axis=1).astype(np.uint8) * bit
        local_hits = np.flatnonzero(flags)
        if not len(local_hits):
            return start, None
        closest = np.argmin(distances[local_hits], axis=1).astype(np.int32)
        max_similarity = (1.0 - distances[local_hits, closest]).astype(np.float32)
        return start, (
            ids[start:end][local_hits],
            max_similarity,
            flags[local_hits],
            closest,
        )

    results = []
    starts = range(0, len(ids), block_size)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for block_index, (start, result) in enumerate(
            executor.map(one_block, starts), start=1
        ):
            if result is not None:
                results.append(result)
            if block_index % 20 == 0:
                print(
                    f"shard {shard}: completed {min(start + block_size, len(ids)):,}/"
                    f"{len(ids):,} fingerprints",
                    flush=True,
                )

    line_id = np.concatenate([result[0] for result in results])
    max_similarity = np.concatenate([result[1] for result in results])
    dataset_flags = np.concatenate([result[2] for result in results])
    closest_query = np.concatenate([result[3] for result in results])
    order = np.argsort(line_id)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        shard=np.int16(shard),
        line_id=line_id[order],
        max_similarity=max_similarity[order],
        dataset_flags=dataset_flags[order],
        closest_query=closest_query[order],
    )
    print(
        f"shard {shard}: wrote {len(line_id):,} exact NumKong hits to {output}",
        flush=True,
    )


def materialize_hits(
    shard: int, hits_path: Path, queries_path: Path, output: Path
) -> None:
    # Materialize every compressed NPZ member once. Indexing an NpzFile member
    # inside the loop would decompress the complete member for every hit.
    with np.load(hits_path) as hit_data:
        line_ids = hit_data["line_id"].astype(np.int64)
        max_similarity = hit_data["max_similarity"].astype(np.float32)
        dataset_flags = hit_data["dataset_flags"].astype(np.uint8)
        closest_query = hit_data["closest_query"].astype(np.int32)
    line_to_offset = {int(line_id): offset for offset, line_id in enumerate(line_ids)}
    queries = load_queries(queries_path)
    rows = []
    with mapping_path(shard).open(buffering=32 * 1024 * 1024) as handle:
        for line_id, line in enumerate(handle, start=1):
            offset = line_to_offset.get(line_id)
            if offset is None:
                continue
            mapping = json.loads(line)
            query = queries[int(closest_query[offset])]
            flag = int(dataset_flags[offset])
            rows.append(
                {
                    "shard": shard,
                    "line_id": line_id,
                    "geom_smiles": mapping["geom_smiles"],
                    "isomeric_smiles": mapping["isomeric_smiles"],
                    "max_similarity": float(max_similarity[offset]),
                    **{
                        f"{key}_hit": bool(flag & bit)
                        for key, bit in DATASET_BITS.items()
                    },
                    "closest_dataset": query["dataset"],
                    "closest_benchmark_id": query["benchmark_id"],
                    "closest_smiles": query["canonical_isomeric_smiles"],
                }
            )
    if len(rows) != len(line_ids):
        raise RuntimeError(f"materialized {len(rows):,} of {len(line_ids):,} hits")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows):,} materialized hits to {output}")


def source_from_name(name: str) -> str:
    # All ordinary source labels precede the first underscore. OMol25 is the
    # only family whose source label itself contains underscores.
    if name.startswith("OMol25_small_mols"):
        return "OMol25_small_mols"
    if name.startswith("OMol25_bio_mols"):
        return "OMol25_bio_mols"
    prefix = name.split("_", 1)[0]
    if prefix in SOURCE_NAME_SET:
        return prefix
    return "UNKNOWN"


def attribute_sources(hit_files: list[Path], output_dir: Path) -> None:
    hits_by_graph: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    hit_rows: dict[tuple[int, int], dict[str, str]] = {}
    for path in hit_files:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                key = (int(row["shard"]), int(row["line_id"]))
                flags = sum(
                    bit
                    for dataset, bit in DATASET_BITS.items()
                    if row.get(f"{dataset}_hit") == "True"
                )
                hits_by_graph[row["geom_smiles"]].append((*key, flags))
                hit_rows[key] = row

    sources_by_hit: dict[tuple[int, int], set[str]] = defaultdict(set)
    removed_conformers = Counter()
    total_conformers = Counter()
    unknown_examples = Counter()
    # The merged CSV has exactly two columns; the second is quoted because
    # coordinates contain commas. Splitting only at the first comma is much
    # faster than csv.DictReader across 190 million rows and remains exact.
    with MERGED_TRAIN.open() as handle:
        header = handle.readline().rstrip("\r\n")
        if header != "name,enriched_text":
            raise RuntimeError(f"unexpected merged CSV header: {header!r}")
        for row_number, line in enumerate(handle, start=1):
            name, enriched_text = line.rstrip("\r\n").split(",", 1)
            if enriched_text.startswith('"') and enriched_text.endswith('"'):
                enriched_text = enriched_text[1:-1]
            source = source_from_name(name)
            total_conformers[source] += 1
            if source == "UNKNOWN" and len(unknown_examples) < 100:
                unknown_examples[name] += 1
            graph = COORDINATE_RE.sub("", enriched_text)
            for shard, line_id, _flags in hits_by_graph.get(graph, ()):
                key = (shard, line_id)
                sources_by_hit[key].add(source)
                removed_conformers[(source, key)] += 1
            if row_number % 5_000_000 == 0:
                print(f"scanned {row_number:,} conformers", flush=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    attribution_path = output_dir / "removed_molecules_with_sources.csv"
    with attribution_path.open("w", newline="") as handle:
        fields = list(next(iter(hit_rows.values())).keys()) + (
            ["source_datasets", "source_count", "training_conformer_count"]
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key, row in sorted(hit_rows.items()):
            sources = sorted(sources_by_hit.get(key, {"UNATTRIBUTED"}))
            writer.writerow(
                {
                    **row,
                    "source_datasets": ";".join(sources),
                    "source_count": len(sources),
                    "training_conformer_count": sum(
                        removed_conformers[(source, key)] for source in sources
                    ),
                }
            )

    summary_rows = []
    all_sources = sorted(set(total_conformers) | {s for values in sources_by_hit.values() for s in values})
    for source in all_sources:
        keys = {key for key, sources in sources_by_hit.items() if source in sources}
        removed_by_dataset = {
            dataset: {
                key for key in keys if hit_rows[key].get(f"{dataset}_hit") == "True"
            }
            for dataset in DATASET_BITS
        }
        summary_rows.append(
            {
                "source_dataset": source,
                "total_conformers": total_conformers[source],
                "removed_unique_molecules_union": len(keys),
                **{
                    f"removed_unique_molecules_{dataset}": len(dataset_keys)
                    for dataset, dataset_keys in removed_by_dataset.items()
                },
                "removed_unique_molecules_posebench_and_platinum": len(
                    removed_by_dataset["posebench"] & removed_by_dataset["platinum"]
                ),
                "removed_conformers_union": sum(
                    removed_conformers[(source, key)] for key in keys
                ),
            }
        )
    summary_path = output_dir / "per_source_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    (output_dir / "source_attribution_metadata.json").write_text(
        json.dumps(
            {
                "merged_training_csv": str(MERGED_TRAIN),
                "hit_files": [str(path) for path in hit_files],
                "unattributed_hits": sum(not value for value in sources_by_hit.values()),
                "unknown_name_examples": dict(unknown_examples),
                "fingerprint": "RDKit Morgan/ECFP4",
                "fingerprint_parameters": FP_PARAMS,
                "metric": "Tanimoto",
                "removal_rule": "maximum similarity to either benchmark > 0.4",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote source attribution to {output_dir}")


def load_hits_by_graph(hit_files: list[Path]):
    hits_by_graph: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in hit_files:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                hits_by_graph[row["geom_smiles"]].append(
                    (int(row["shard"]), int(row["line_id"]))
                )
    return hits_by_graph


def load_unattributed_by_noniso(path: Path | None):
    if path is None:
        return {}
    from rdkit import Chem

    targets: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["source_datasets"] != "UNATTRIBUTED":
                continue
            mol = Chem.MolFromSmiles(row["geom_smiles"])
            if mol is None:
                raise ValueError(f"failed to parse unattributed hit: {row['geom_smiles']}")
            stereo_key = Chem.MolToInchiKey(mol)
            Chem.RemoveStereochemistry(mol)
            noniso = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
            targets[noniso].append(
                (int(row["shard"]), int(row["line_id"]), stereo_key)
            )
    return targets


def attribute_source_chunk(
    chunk_index: int,
    num_chunks: int,
    hit_files: list[Path],
    output: Path,
    workers: int,
    unattributed_file: Path | None,
    source_identity_dir: Path | None,
    training_identity_hashes: Path | None,
) -> None:
    global _SOURCE_HITS_BY_SMILES, _SOURCE_UNATTRIBUTED_BY_NONISO
    global _TRAINING_IDENTITY_HASHES
    _SOURCE_HITS_BY_SMILES = load_hits_by_graph(hit_files)
    _SOURCE_UNATTRIBUTED_BY_NONISO = load_unattributed_by_noniso(unattributed_file)
    _TRAINING_IDENTITY_HASHES = (
        np.memmap(training_identity_hashes, dtype="V16", mode="r")
        if training_identity_hashes is not None
        else None
    )
    file_size = MERGED_TRAIN.stat().st_size
    start = file_size * chunk_index // num_chunks
    end = file_size * (chunk_index + 1) // num_chunks
    hit_conformers: Counter = Counter()
    total_conformers = Counter()
    unknown_examples = Counter()

    def iter_work():
        nonlocal unknown_examples
        batch: list[tuple[str, str]] = []
        with MERGED_TRAIN.open("rb") as handle:
            if start == 0:
                handle.seek(0)
                header = handle.readline().rstrip(b"\r\n")
                if header != b"name,enriched_text":
                    raise RuntimeError(f"unexpected merged CSV header: {header!r}")
            else:
                handle.seek(start - 1)
                previous = handle.read(1)
                if previous != b"\n":
                    handle.readline()
            rows = 0
            while handle.tell() < end:
                line = handle.readline()
                if not line:
                    break
                rows += 1
                name_bytes, enriched = line.rstrip(b"\r\n").split(b",", 1)
                if enriched.startswith(b'"') and enriched.endswith(b'"'):
                    enriched = enriched[1:-1]
                name = name_bytes.decode("utf-8")
                source = source_from_name(name)
                total_conformers[source] += 1
                if source == "UNKNOWN" and len(unknown_examples) < 100:
                    unknown_examples[name] += 1
                batch.append((source, enriched.decode("utf-8")))
                if len(batch) == 500:
                    yield batch
                    batch = []
                if rows % 5_000_000 == 0:
                    print(
                        f"chunk {chunk_index}: fed {rows:,} rows at byte {handle.tell():,}",
                        flush=True,
                    )
            if batch:
                yield batch
        iter_work.rows = rows

    canonicalization_failures = 0
    identity_handles = {}
    if source_identity_dir is not None:
        source_identity_dir.mkdir(parents=True, exist_ok=True)
        identity_handles = {
            source: (source_identity_dir / f"{source}.raw.bin").open("wb")
            for source in (*SOURCE_NAMES, "UNKNOWN")
        }
    # Linux fork lets all workers share the read-only hit index copy-on-write.
    context = mp.get_context("fork")
    with context.Pool(processes=workers, maxtasksperchild=1000) as pool:
        worker = (
            partial(_canonical_source_chunk, collect_identities=True)
            if source_identity_dir is not None
            else _canonical_source_chunk
        )
        for partial_counts, failures, identity_hashes in pool.imap_unordered(
            worker, iter_work(), chunksize=4
        ):
            hit_conformers.update(partial_counts)
            canonicalization_failures += failures
            for source, hashes in identity_hashes.items():
                identity_handles[source].write(hashes)

    for handle in identity_handles.values():
        handle.close()

    unique_molecules_by_source = {}
    if source_identity_dir is not None:
        for source in (*SOURCE_NAMES, "UNKNOWN"):
            raw_path = source_identity_dir / f"{source}.raw.bin"
            values = np.fromfile(raw_path, dtype="V16")
            unique_values = np.unique(values)
            unique_values.tofile(source_identity_dir / f"{source}.unique.bin")
            unique_molecules_by_source[source] = int(len(unique_values))
            raw_path.unlink()

    rows = getattr(iter_work, "rows", sum(total_conformers.values()))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        fields = ("source_dataset", "shard", "line_id", "conformer_count")
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for (source, shard, line_id), count in sorted(hit_conformers.items()):
            writer.writerow(
                {
                    "source_dataset": source,
                    "shard": shard,
                    "line_id": line_id,
                    "conformer_count": count,
                }
            )
    output.with_suffix(".metadata.json").write_text(
        json.dumps(
            {
                "chunk_index": chunk_index,
                "num_chunks": num_chunks,
                "start_byte": start,
                "end_byte": end,
                "rows_scanned": rows,
                "canonicalization_failures": canonicalization_failures,
                "total_conformers_by_source": dict(total_conformers),
                "unique_molecules_by_source_in_chunk": unique_molecules_by_source,
                "unknown_name_examples": dict(unknown_examples),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"chunk {chunk_index}: wrote {len(hit_conformers):,} hit-source rows to {output}")


def merge_source_chunks(
    hit_files: list[Path],
    partial_files: list[Path],
    output_dir: Path,
    source_identity_dirs: list[Path],
    training_identity_hashes: Path | None,
) -> None:
    hit_rows: dict[tuple[int, int], dict[str, str]] = {}
    for path in hit_files:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                hit_rows[(int(row["shard"]), int(row["line_id"]))] = row

    conformers_by_source_hit = Counter()
    total_conformers = Counter()
    unknown_examples = Counter()
    for path in partial_files:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                key = (
                    row["source_dataset"],
                    int(row["shard"]),
                    int(row["line_id"]),
                )
                conformers_by_source_hit[key] += int(row["conformer_count"])
        metadata = json.loads(path.with_suffix(".metadata.json").read_text())
        total_conformers.update(metadata["total_conformers_by_source"])
        unknown_examples.update(metadata["unknown_name_examples"])

    total_unique_molecules = {}
    all_source_unique_arrays = []
    if source_identity_dirs:
        for source in (*SOURCE_NAMES, "UNKNOWN"):
            arrays = [
                np.fromfile(path / f"{source}.unique.bin", dtype="V16")
                for path in source_identity_dirs
                if (path / f"{source}.unique.bin").exists()
            ]
            source_values = np.unique(np.concatenate(arrays)) if arrays else np.empty(0, dtype="V16")
            total_unique_molecules[source] = int(len(source_values))
            if len(source_values):
                all_source_unique_arrays.append(source_values)

    sources_by_hit: dict[tuple[int, int], set[str]] = defaultdict(set)
    for source, shard, line_id in conformers_by_source_hit:
        sources_by_hit[(shard, line_id)].add(source)

    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "removed_molecules_with_sources.csv"
    fields = list(next(iter(hit_rows.values())).keys()) + [
        "source_datasets",
        "source_count",
        "training_conformer_count",
    ]
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key, row in sorted(hit_rows.items()):
            sources = sorted(sources_by_hit.get(key, {"UNATTRIBUTED"}))
            writer.writerow(
                {
                    **row,
                    "source_datasets": ";".join(sources),
                    "source_count": len(sources),
                    "training_conformer_count": sum(
                        conformers_by_source_hit[(source, *key)] for source in sources
                    ),
                }
            )

    summary_rows = []
    all_sources = sorted(
        set(total_conformers) | {source for sources in sources_by_hit.values() for source in sources}
    )
    for source in all_sources:
        keys = {key for key, sources in sources_by_hit.items() if source in sources}
        removed_by_dataset = {
            dataset: {
                key for key in keys if hit_rows[key].get(f"{dataset}_hit") == "True"
            }
            for dataset in DATASET_BITS
        }
        denominator = total_unique_molecules.get(source)
        summary_rows.append(
            {
                "source_dataset": source,
                "total_conformers": total_conformers[source],
                "total_unique_molecules": denominator,
                "removed_unique_molecules_union": len(keys),
                "removed_unique_molecules_union_percent": (
                    100.0 * len(keys) / denominator if denominator else None
                ),
                **{
                    f"removed_unique_molecules_{dataset}": len(dataset_keys)
                    for dataset, dataset_keys in removed_by_dataset.items()
                },
                **{
                    f"removed_unique_molecules_{dataset}_percent": (
                        100.0 * len(dataset_keys) / denominator if denominator else None
                    )
                    for dataset, dataset_keys in removed_by_dataset.items()
                },
                "removed_unique_molecules_posebench_and_platinum": len(
                    removed_by_dataset["posebench"] & removed_by_dataset["platinum"]
                ),
                "removed_conformers_union": sum(
                    conformers_by_source_hit[(source, *key)] for key in keys
                ),
            }
        )
    summary_path = output_dir / "per_source_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    overall_denominator = (
        training_identity_hashes.stat().st_size // 16
        if training_identity_hashes is not None
        else None
    )
    attributed_training_unique_molecules = (
        int(len(np.unique(np.concatenate(all_source_unique_arrays))))
        if all_source_unique_arrays
        else None
    )
    overall_by_dataset = {
        dataset: {
            key
            for key, row in hit_rows.items()
            if row.get(f"{dataset}_hit") == "True"
        }
        for dataset in DATASET_BITS
    }
    overall_row = {
        "total_unique_molecules": overall_denominator,
        "attributed_training_unique_molecules": attributed_training_unique_molecules,
        "unattributed_training_unique_molecules": (
            overall_denominator - attributed_training_unique_molecules
            if overall_denominator is not None
            and attributed_training_unique_molecules is not None
            else None
        ),
        "removed_unique_molecules_union": len(hit_rows),
        "removed_unique_molecules_union_percent": (
            100.0 * len(hit_rows) / overall_denominator
            if overall_denominator
            else None
        ),
        **{
            f"removed_unique_molecules_{dataset}": len(keys)
            for dataset, keys in overall_by_dataset.items()
        },
        **{
            f"removed_unique_molecules_{dataset}_percent": (
                100.0 * len(keys) / overall_denominator
                if overall_denominator
                else None
            )
            for dataset, keys in overall_by_dataset.items()
        },
    }
    with (output_dir / "overall_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=overall_row.keys())
        writer.writeheader()
        writer.writerow(overall_row)
    (output_dir / "source_attribution_metadata.json").write_text(
        json.dumps(
            {
                "merged_training_csv": str(MERGED_TRAIN),
                "hit_files": [str(path) for path in hit_files],
                "partial_files": [str(path) for path in partial_files],
                "unattributed_hits": sum(key not in sources_by_hit for key in hit_rows),
                "unknown_name_examples": dict(unknown_examples),
                "fingerprint": "RDKit Morgan/ECFP4",
                "fingerprint_parameters": FP_PARAMS,
                "metric": "Tanimoto",
                "removal_rule": "maximum similarity to either benchmark > 0.4",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"merged {len(partial_files)} source chunks into {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare-queries")
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--audit-module-dir", type=Path, required=True)
    prepare.add_argument("--extra-query-files", type=Path, nargs="*", default=[])
    prepare.add_argument("--casf-core-mapping", type=Path)
    prepare.add_argument("--plinder-shortlist", type=Path)
    prepare.add_argument("--exclude-standard-benchmarks", action="store_true")

    build = subparsers.add_parser("build-db")
    build.add_argument("--shard", type=int, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--limit", type=int)

    identity_build = subparsers.add_parser("build-training-identity-hashes")
    identity_build.add_argument("--shard", type=int, required=True)
    identity_build.add_argument("--output", type=Path, required=True)

    identity_merge = subparsers.add_parser("merge-training-identity-hashes")
    identity_merge.add_argument("--inputs", type=Path, nargs="+", required=True)
    identity_merge.add_argument("--output", type=Path, required=True)

    search = subparsers.add_parser("search-db")
    search.add_argument("--shard", type=int, required=True)
    search.add_argument("--database", type=Path, required=True)
    search.add_argument("--queries", type=Path, required=True)
    search.add_argument("--output", type=Path, required=True)
    search.add_argument("--workers", type=int, default=16)

    numkong_search = subparsers.add_parser("search-db-numkong")
    numkong_search.add_argument("--shard", type=int, required=True)
    numkong_search.add_argument("--database", type=Path, required=True)
    numkong_search.add_argument("--queries", type=Path, required=True)
    numkong_search.add_argument("--output", type=Path, required=True)
    numkong_search.add_argument("--workers", type=int, default=8)
    numkong_search.add_argument("--block-size", type=int, default=25_000)

    materialize = subparsers.add_parser("materialize-hits")
    materialize.add_argument("--shard", type=int, required=True)
    materialize.add_argument("--hits", type=Path, required=True)
    materialize.add_argument("--queries", type=Path, required=True)
    materialize.add_argument("--output", type=Path, required=True)

    attribute = subparsers.add_parser("attribute-sources")
    attribute.add_argument("--hit-files", type=Path, nargs="+", required=True)
    attribute.add_argument("--output-dir", type=Path, required=True)

    chunk = subparsers.add_parser("attribute-source-chunk")
    chunk.add_argument("--chunk-index", type=int, required=True)
    chunk.add_argument("--num-chunks", type=int, required=True)
    chunk.add_argument("--hit-files", type=Path, nargs="+", required=True)
    chunk.add_argument("--output", type=Path, required=True)
    chunk.add_argument("--workers", type=int, default=16)
    chunk.add_argument("--unattributed-file", type=Path)
    chunk.add_argument("--source-identity-dir", type=Path)
    chunk.add_argument("--training-identity-hashes", type=Path)

    merge = subparsers.add_parser("merge-source-chunks")
    merge.add_argument("--hit-files", type=Path, nargs="+", required=True)
    merge.add_argument("--partial-files", type=Path, nargs="+", required=True)
    merge.add_argument("--output-dir", type=Path, required=True)
    merge.add_argument("--source-identity-dirs", type=Path, nargs="*", default=[])
    merge.add_argument("--training-identity-hashes", type=Path)

    args = parser.parse_args()
    if args.command == "prepare-queries":
        prepare_queries(
            args.output,
            args.audit_module_dir,
            args.extra_query_files,
            args.casf_core_mapping,
            args.plinder_shortlist,
            not args.exclude_standard_benchmarks,
        )
    elif args.command == "build-db":
        build_db(args.shard, args.output, args.limit)
    elif args.command == "build-training-identity-hashes":
        build_training_identity_hashes(args.shard, args.output)
    elif args.command == "merge-training-identity-hashes":
        merge_training_identity_hashes(args.inputs, args.output)
    elif args.command == "search-db":
        search_db(args.shard, args.database, args.queries, args.output, args.workers)
    elif args.command == "search-db-numkong":
        search_db_numkong(
            args.shard,
            args.database,
            args.queries,
            args.output,
            args.workers,
            args.block_size,
        )
    elif args.command == "materialize-hits":
        materialize_hits(args.shard, args.hits, args.queries, args.output)
    elif args.command == "attribute-sources":
        attribute_sources(args.hit_files, args.output_dir)
    elif args.command == "attribute-source-chunk":
        attribute_source_chunk(
            args.chunk_index,
            args.num_chunks,
            args.hit_files,
            args.output,
            args.workers,
            args.unattributed_file,
            args.source_identity_dir,
            args.training_identity_hashes,
        )
    elif args.command == "merge-source-chunks":
        merge_source_chunks(
            args.hit_files,
            args.partial_files,
            args.output_dir,
            args.source_identity_dirs,
            args.training_identity_hashes,
        )


if __name__ == "__main__":
    main()
