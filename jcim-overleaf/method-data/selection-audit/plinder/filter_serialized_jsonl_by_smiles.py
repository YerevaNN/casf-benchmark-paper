#!/usr/bin/env python3
"""Stream-filter serialized conformer JSONL files by canonical SMILES.

The output keeps every retained JSONL line byte-for-byte unchanged and writes
the uint64 byte-offset ``.idx`` files expected by the MolGen3D dataloader.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from array import array
from pathlib import Path


CANONICAL_PREFIX = b'{"canonical_smiles":'
EMBEDDED_DELIMITER = b',"embedded_smiles":'


def json_string_token(value: str) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode()


def load_exclusions(path: Path) -> set[bytes]:
    exclusions: set[bytes] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "geom_smiles" not in (reader.fieldnames or ()):
            raise ValueError(f"missing geom_smiles column in {path}")
        for row in reader:
            exclusions.add(json_string_token(row["geom_smiles"]))
    return exclusions


def canonical_token(line: bytes) -> bytes:
    if line.startswith(CANONICAL_PREFIX):
        end = line.find(EMBEDDED_DELIMITER, len(CANONICAL_PREFIX))
        if end >= 0:
            return line[len(CANONICAL_PREFIX) : end]
    row = json.loads(line)
    return json_string_token(row["canonical_smiles"])


def filter_file(source: Path, destination: Path, exclusions: set[bytes]) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    partial_idx = Path(str(partial) + ".idx")
    offsets = array("Q")
    input_rows = removed_rows = retained_rows = output_bytes = 0

    with source.open("rb") as src, partial.open("wb") as dst:
        for line in src:
            input_rows += 1
            if canonical_token(line) in exclusions:
                removed_rows += 1
                continue
            offsets.append(output_bytes)
            dst.write(line)
            output_bytes += len(line)
            retained_rows += 1

    with partial_idx.open("wb") as handle:
        offsets.tofile(handle)
    os.replace(partial, destination)
    os.replace(partial_idx, Path(str(destination) + ".idx"))

    metadata = {
        "source": str(source),
        "destination": str(destination),
        "input_rows": input_rows,
        "removed_rows": removed_rows,
        "retained_rows": retained_rows,
        "output_bytes": output_bytes,
    }
    metadata_path = destination.with_name(destination.name + ".filter.json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--num-tasks", type=int, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("*.jsonl"))
    assigned = files[args.task_index :: args.num_tasks]
    exclusions = load_exclusions(args.exclusions)
    totals = {
        "task_index": args.task_index,
        "num_tasks": args.num_tasks,
        "exclusion_molecules": len(exclusions),
        "assigned_files": len(assigned),
        "processed_files": 0,
        "skipped_files": 0,
        "input_rows": 0,
        "removed_rows": 0,
        "retained_rows": 0,
        "output_bytes": 0,
    }

    for file_index, source in enumerate(assigned, start=1):
        destination = args.output_dir / source.name
        metadata_path = destination.with_name(destination.name + ".filter.json")
        if args.resume and destination.exists() and Path(str(destination) + ".idx").exists() and metadata_path.exists():
            metadata = json.loads(metadata_path.read_text())
            totals["skipped_files"] += 1
        else:
            metadata = filter_file(source, destination, exclusions)
            totals["processed_files"] += 1
        for key in ("input_rows", "removed_rows", "retained_rows", "output_bytes"):
            totals[key] += metadata[key]
        print(
            f"task {args.task_index}: {file_index}/{len(assigned)} {source.name} "
            f"removed={metadata['removed_rows']:,}",
            flush=True,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = args.output_dir / f"task_{args.task_index:02d}.summary.json"
    summary.write_text(json.dumps(totals, indent=2) + "\n")
    print(json.dumps(totals, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
