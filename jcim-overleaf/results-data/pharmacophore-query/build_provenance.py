"""Build the compact pharmacophore-query archive provenance record."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent

COPIED_SOURCES = {
    "pilot-aggregate-budget-curves.csv": "data/results/pharmacophore_plinder22/aggregate_budget_curves.csv",
    "pilot-per-molecule-budget-curves.csv": "data/results/pharmacophore_plinder22/per_molecule_budget_curves.csv",
    "pilot-endpoint-comparisons.csv": "data/results/pharmacophore_plinder22/endpoint_comparisons.csv",
    "pilot-coverage-breakdowns.csv": "data/results/pharmacophore_plinder22/coverage_breakdowns.csv",
    "pilot-query-manifest.csv": "data/results/pharmacophore_plinder22/query_manifest.csv",
    "pilot-query-freeze-manifest.json": "data/results/pharmacophore_plinder22/query_freeze_manifest.json",
    "pilot-match-manifest.json": "data/results/pharmacophore_plinder22/match_manifest.json",
    "pilot-summary-manifest.json": "data/results/pharmacophore_plinder22/summary_manifest.json",
    "pilot-pool-inventory.csv": "data/results/pharmacophore_plinder22/pool_inventory_completed.csv",
    "pilot-posebusters-summary.csv": "data/results/pharmacophore_plinder22/posebusters/summary.csv",
    "pilot-posebusters-manifest.json": "data/results/pharmacophore_plinder22/posebusters/manifest.json",
    "pilot-rmsd-molecule-summary.csv": "data/results/pharmacophore_plinder22/rmsd_agreement/molecule_weighted_summary.csv",
    "pilot-rmsd-summary.csv": "data/results/pharmacophore_plinder22/rmsd_agreement/summary.csv",
    "pilot-rmsd-manifest.json": "data/results/pharmacophore_plinder22/rmsd_agreement/manifest.json",
    "pilot-rmsd-cache-manifest.csv": "data/results/pharmacophore_plinder22/rmsd_agreement/cache_manifest.csv",
    "pilot-qwen-topup-manifest.json": "data/results/pharmacophore_plinder22/qwen_topup/topup_manifest.json",
    "pilot-config.yaml": "src/casf_benchmark/config/pharmacophore_plinder22.yaml",
    "pilot-findings.md": "data/results/pharmacophore_plinder22/FINDINGS.md",
    "atc48-aggregate-budget-curves.csv": "data/results/pharmacophore_atc48/aggregate_budget_curves.csv",
    "atc48-per-molecule-budget-curves.csv": "data/results/pharmacophore_atc48/per_molecule_budget_curves.csv",
    "atc48-endpoint-comparisons.csv": "data/results/pharmacophore_atc48/endpoint_comparisons.csv",
    "atc48-coverage-breakdowns.csv": "data/results/pharmacophore_atc48/coverage_breakdowns.csv",
    "atc48-query-manifest.csv": "data/results/pharmacophore_atc48/query_manifest.csv",
    "atc48-query-freeze-manifest.json": "data/results/pharmacophore_atc48/query_freeze_manifest.json",
    "atc48-match-manifest.json": "data/results/pharmacophore_atc48/match_manifest.json",
    "atc48-summary-manifest.json": "data/results/pharmacophore_atc48/summary_manifest.json",
    "atc48-pool-inventory.csv": "data/results/pharmacophore_atc48/pool_inventory_completed.csv",
    "atc48-inventory-manifest.json": "data/results/pharmacophore_atc48/inventory_manifest.json",
    "atc48-loqi-generation-summary.json": "data/results/pharmacophore_atc48/generation/loqi/generation_summary.json",
    "atc48-loqi-run-manifest.json": "data/results/pharmacophore_atc48/generation/loqi/run_manifest.json",
    "atc48-rdkit-generation-manifest.json": "data/results/pharmacophore_atc48/generation/rdkit/generation_manifest.json",
    "atc48-posebusters-summary.csv": "data/results/pharmacophore_atc48/posebusters/summary.csv",
    "atc48-posebusters-manifest.json": "data/results/pharmacophore_atc48/posebusters/manifest.json",
    "atc48-rmsd-molecule-summary.csv": "data/results/pharmacophore_atc48/rmsd_agreement/molecule_weighted_summary.csv",
    "atc48-rmsd-summary.csv": "data/results/pharmacophore_atc48/rmsd_agreement/summary.csv",
    "atc48-rmsd-manifest.json": "data/results/pharmacophore_atc48/rmsd_agreement/manifest.json",
    "atc48-rmsd-cache-manifest.csv": "data/results/pharmacophore_atc48/rmsd_agreement/cache_manifest.csv",
    "atc48-findings.md": "data/results/pharmacophore_atc48/FINDINGS.md",
    "atc48-cohort.csv": "data/pharmacophore_atc48/cohort.csv",
    "atc48-cohort-selection-audit.csv": "data/pharmacophore_atc48/cohort_selection_audit.csv",
    "atc48-reference-audit.csv": "data/pharmacophore_atc48/reference_audit.csv",
    "atc48-freeze-manifest.json": "data/pharmacophore_atc48/freeze_manifest.json",
    "atc48-exact-identity-overlap.csv": "data/results/pharmacophore_atc48/training_overlap/exact_identity_overlap.csv",
    "atc48-training-overlap-manifest.json": "data/results/pharmacophore_atc48/training_overlap/manifest.json",
    "atc48-config.yaml": "src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml",
}

CODE_SOURCES = [
    "src/casf_benchmark/analysis/pharmacophore.py",
    "src/casf_benchmark/cli/analyze_pharmacophores.py",
    "scripts/build_pharmit_pilot.sh",
    "scripts/pharmit_exact_search.patch",
    "scripts/build_pharmacophore_atc48.py",
    "scripts/audit_pharmacophore_training_overlap.py",
    "scripts/audit_pharmacophore_training_overlap.sbatch",
    "scripts/generate_pharmacophore_baselines.py",
    "scripts/generate_pharmacophore_baselines.sbatch",
    "scripts/topup_pharmacophore_qwen.py",
    "scripts/topup_pharmacophore_qwen.sbatch",
    "scripts/run_off_the_shelf_druglike.py",
    "scripts/audit_pharmacophore_loqi_generation.py",
    "scripts/run_pharmacophore_loqi_cpu_shard.sbatch",
    "scripts/run_pharmacophore_posebusters.sbatch",
    "docs/pharmacophore_pilot/REPRODUCING.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--benchmark-root",
        type=Path,
        default=Path("/mnt/weka/mbedrosian/code/casf-benchmark"),
    )
    args = parser.parse_args()
    benchmark = args.benchmark_root.resolve()

    copied = {}
    for archived_name, source_name in COPIED_SOURCES.items():
        archived = HERE / archived_name
        source = benchmark / source_name
        archived_record = file_record(archived)
        source_record = file_record(source)
        if archived_record != source_record:
            raise RuntimeError(f"Archived copy differs from source: {archived_name}")
        copied[archived_name] = {"source": str(source), **archived_record}

    code = {name: file_record(benchmark / name) for name in CODE_SOURCES}
    pools = {
        "qwen_atc48": Path(
            "/mnt/weka/mbedrosian/pharma_generation_analysis/"
            "pharmacophore_atc48_qwen17b/gen_results/"
            "qwen_1p7b_v2_coordtoken_step41054_atc48_plinder/generation_results.pickle"
        ),
        "loqi_atc48": benchmark
        / "data/results/pharmacophore_atc48/generation/loqi/generation_results.pickle",
        "rdkit_raw_atc48": benchmark
        / "data/results/pharmacophore_atc48/generation/rdkit/rdkit_etkdg_raw.pickle",
        "rdkit_mmff94s_atc48": benchmark
        / "data/results/pharmacophore_atc48/generation/rdkit/rdkit_etkdg_mmff94s.pickle",
        "pharmit_executable": Path(
            "/mnt/weka/mbedrosian/environments/casf-pharmit-pilot/bin/pharmit"
        ),
    }
    pool_records = {name: {"path": str(path), **file_record(path)} for name, path in pools.items()}
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=benchmark,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    working_tree_dirty = bool(
        subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=benchmark,
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
    )

    derived_files = [
        "atc48-training-overlap-endpoints.csv",
        "build_overlap_summary.py",
        "build_provenance.py",
        "render.py",
        "rendered-files.json",
        "verify.py",
    ]
    record = {
        "analysis_date": "2026-11-17",
        "benchmark_repository": {
            "path": str(benchmark),
            "base_commit": commit,
            "working_tree_dirty": working_tree_dirty,
            "code_source_hashes_are_authoritative": True,
        },
        "panels": {
            "PLINDER-22": {"molecules": 22, "references": 1067, "eligible_queries": 932},
            "PLINDER-ATC48": {
                "selected_molecules": 48,
                "references": 379,
                "eligible_molecules": 44,
                "eligible_queries": 311,
            },
        },
        "matching": {
            "feature_backend": "Pharmit",
            "pharmit_commit": "49d2c0a573f6e58bca62f674b33874bdd42c2274",
            "membership": "exhaustive typed injective correspondence plus proper Kabsch",
            "indexed_search_used_for_final_membership": False,
            "radii_angstrom": [0.5, 0.75, 1.0],
            "radius_monotonicity_violations": 0,
        },
        "sampling": {
            "budgets": [10, 20, 50, 100, 200, 500, 1000],
            "finite_pool_sampling": "exact uniform without replacement",
            "mixture": "equal Qwen/LoQI allocation; 500 + 500 at budget 1000",
            "bootstrap_replicates": 2000,
            "bootstrap_seed": 20261116,
            "posebusters_failures_consume_budget": True,
        },
        "training_identity_audit": {
            "definition": "exact isomeric InChIKey",
            "overall_exact_identities": "20/48",
            "eligible_exact_identities": "17/44",
            "eligible_without_exact_identity": "27/44",
            "limitation": "does not exclude scaffold or near-neighbor overlap",
        },
        "copied_results": copied,
        "code_sources": code,
        "external_artifacts": pool_records,
        "derived_archive_files": {name: file_record(HERE / name) for name in derived_files},
    }
    (HERE / "provenance.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
