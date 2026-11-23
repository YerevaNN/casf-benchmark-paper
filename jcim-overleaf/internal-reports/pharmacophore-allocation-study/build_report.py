#!/usr/bin/env python3
"""Build the internal report for the PLINDER-22 pharmacophore allocation study."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter


REPORT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = REPORT_DIR / "figures"
PAPER_ROOT = REPORT_DIR.parents[1]
ARCHIVE_DIR = PAPER_ROOT / "results-data" / "cdpkit-pharmacophore-mixture"
ANALYSIS_DIR = Path(
    "/mnt/weka/mbedrosian/code/casf-benchmark/data/results/"
    "cdpkit_plinder22_exploration_v1"
)
ANALYSIS_TABLES = ANALYSIS_DIR / "tables"

BUDGETS = [20, 100, 500, 1000]
FAMILIES = [
    "individual_ligand",
    "individual_interaction",
    "shared_ligand",
    "shared_protein",
]
FAMILY_LABELS = {
    "individual_ligand": "Individual ligand",
    "individual_interaction": "Individual interaction",
    "shared_ligand": "Shared ligand",
    "shared_protein": "Shared protein",
}
RECIPES = [
    "qwen",
    "loqi",
    "qwen50_loqi50",
    "qwen75_loqi25",
    "qwen25_loqi75",
    "rdkit_raw",
    "rdkit_minimized",
]
RECIPE_LABELS = {
    "qwen": "Qwen",
    "loqi": "LoQI",
    "qwen50_loqi50": "Qwen/LoQI 50:50",
    "qwen75_loqi25": "Qwen/LoQI 75:25",
    "qwen25_loqi75": "Qwen/LoQI 25:75",
    "rdkit_raw": "RDKit raw",
    "rdkit_minimized": "RDKit minimized",
}
COLORS = {
    "qwen": "#0072B2",
    "loqi": "#E69F00",
    "qwen50_loqi50": "#8E44AD",
    "qwen75_loqi25": "#009E73",
    "qwen25_loqi75": "#CC79A7",
    "rdkit_raw": "#6C757D",
    "rdkit_minimized": "#343A40",
}
LINESTYLES = {
    "qwen": "-",
    "loqi": "-",
    "qwen50_loqi50": "-",
    "qwen75_loqi25": "-",
    "qwen25_loqi75": "-",
    "rdkit_raw": "--",
    "rdkit_minimized": ":",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def use_report_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#D9DEE3",
            "grid.linewidth": 0.7,
            "grid.alpha": 0.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def primary_mask(frame: pd.DataFrame) -> pd.Series:
    individual = frame["family"].str.startswith("individual")
    shared = frame["family"].str.startswith("shared")
    return (
        (individual & (frame["construction_count"] == 1))
        | (
            shared
            & (frame["construction_count"] == 2)
            & (frame["support_rule"] == "all")
        )
    )


def plot_strict_coverage(strict: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4), sharex=True, sharey=True)
    for ax, family in zip(axes.flat, FAMILIES):
        family_data = strict[strict["family"] == family]
        for recipe in RECIPES:
            values = family_data[family_data["recipe"] == recipe].sort_values("budget")
            ax.plot(
                values["budget"],
                100 * values["coverage"],
                marker="o",
                markersize=4.2,
                linewidth=2,
                color=COLORS[recipe],
                linestyle=LINESTYLES[recipe],
                label=RECIPE_LABELS[recipe],
            )
            ax.fill_between(
                values["budget"].to_numpy(),
                (100 * values["ci_low"]).to_numpy(),
                (100 * values["ci_high"]).to_numpy(),
                color=COLORS[recipe],
                alpha=0.045,
                linewidth=0,
            )
        ax.set_title(FAMILY_LABELS[family], loc="left", fontweight="bold")
        ax.set_xscale("log")
        ax.set_xticks(BUDGETS)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_ylim(0, 102)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_ylabel("Coverage (%)")
        ax.set_xlabel("Cached-conformer budget")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.02),
        columnspacing=1.5,
        handlelength=2.8,
    )
    fig.suptitle(
        "Stringent pharmacophore coverage across all four query families",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.965,
        "Six of six features within 0.5 Å; shaded regions are 95% molecule-bootstrap intervals",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    save_figure(fig, "strict-coverage-by-budget")


def plot_full_grid_heatmap(primary: pd.DataFrame) -> None:
    means = (
        primary.groupby(["family", "recipe", "budget"], as_index=False)["coverage"]
        .mean()
        .assign(coverage=lambda x: 100 * x["coverage"])
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.5))
    norm = Normalize(vmin=70, vmax=100)
    image = None
    for panel_index, (ax, family) in enumerate(zip(axes.flat, FAMILIES)):
        matrix = (
            means[means["family"] == family]
            .pivot(index="recipe", columns="budget", values="coverage")
            .reindex(index=RECIPES, columns=BUDGETS)
        )
        image = ax.imshow(matrix.to_numpy(), cmap="YlGnBu", norm=norm, aspect="auto")
        ax.grid(False)
        ax.set_title(FAMILY_LABELS[family], loc="left", fontweight="bold")
        ax.set_xticks(np.arange(len(BUDGETS)), [str(x) for x in BUDGETS])
        ax.set_yticks(np.arange(len(RECIPES)))
        if panel_index % 2 == 0:
            ax.set_yticklabels([RECIPE_LABELS[recipe] for recipe in RECIPES])
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("Cached-conformer budget")
        for row_index in range(matrix.shape[0]):
            for col_index in range(matrix.shape[1]):
                value = matrix.iloc[row_index, col_index]
                ax.text(
                    col_index,
                    row_index,
                    f"{value:.1f}",
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    color="white" if value >= 88 else "#17212B",
                )
    assert image is not None
    color_axis = fig.add_axes([0.925, 0.14, 0.016, 0.72])
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label("Mean coverage across 21 query definitions (%)")
    fig.suptitle(
        "Complete primary-grid sensitivity summary",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.965,
        "Each cell averages seven feature requirements and three positional tolerances",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.subplots_adjust(
        top=0.91,
        bottom=0.08,
        left=0.17,
        right=0.90,
        hspace=0.30,
        wspace=0.18,
    )
    save_figure(fig, "full-grid-mean-coverage")


def plot_allocation_difference(primary: pd.DataFrame) -> None:
    index = [
        "family",
        "condition",
        "radius_angstrom",
        "construction_count",
        "support_rule",
        "budget",
    ]
    pivot = primary.pivot_table(index=index, columns="recipe", values="coverage").reset_index()
    pivot["difference_pp"] = 100 * (
        pivot["qwen75_loqi25"] - pivot["qwen50_loqi50"]
    )
    rng = np.random.default_rng(20261123)
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.8), sharex=True, sharey=True)
    for panel_index, (ax, family) in enumerate(zip(axes.flat, FAMILIES)):
        family_data = pivot[pivot["family"] == family]
        values = [
            family_data[family_data["budget"] == budget]["difference_pp"].to_numpy()
            for budget in BUDGETS
        ]
        boxes = ax.boxplot(
            values,
            positions=np.arange(1, 5),
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#17212B", "linewidth": 1.4},
            whiskerprops={"color": "#6C757D"},
            capprops={"color": "#6C757D"},
        )
        for box in boxes["boxes"]:
            box.set_facecolor("#A8D5C2")
            box.set_edgecolor(COLORS["qwen75_loqi25"])
        for position, points in enumerate(values, start=1):
            jitter = rng.uniform(-0.16, 0.16, size=len(points))
            ax.scatter(
                position + jitter,
                points,
                s=13,
                color="#4F5B66",
                alpha=0.45,
                linewidth=0,
                zorder=2,
            )
            ax.scatter(
                position,
                np.mean(points),
                marker="D",
                s=32,
                color="#D55E00",
                edgecolor="white",
                linewidth=0.5,
                zorder=4,
            )
        ax.axhline(0, color="#17212B", linewidth=1, linestyle="--")
        ax.set_title(FAMILY_LABELS[family], loc="left", fontweight="bold")
        ax.set_xticks(np.arange(1, 5), [str(x) for x in BUDGETS])
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_xlabel("Cached-conformer budget")
        if panel_index % 2 == 0:
            ax.set_ylabel("75:25 minus 50:50 coverage (percentage points)")
        ax.set_ylim(-7, 3)
    fig.suptitle(
        "The 50:50 and 75:25 ordering changes across query definitions",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.965,
        "Each box contains all 21 feature-requirement and tolerance combinations; diamonds are means",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, "allocation-difference-full-grid")


def plot_diversity(diversity: pd.DataFrame) -> None:
    selected = diversity[diversity["cluster_radius"] == 0.5]
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), sharex=True)
    measures = [
        ("clusters", "Mean greedy cluster count"),
        ("effective_clusters", "Effective cluster count"),
    ]
    for ax, (column, label) in zip(axes, measures):
        for recipe in RECIPES:
            values = selected[selected["recipe"] == recipe].sort_values("budget")
            ax.plot(
                values["budget"],
                values[column],
                marker="o",
                markersize=4.2,
                linewidth=2,
                color=COLORS[recipe],
                linestyle=LINESTYLES[recipe],
                label=RECIPE_LABELS[recipe],
            )
        ax.set_xscale("log")
        ax.set_xticks(BUDGETS)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.set_xlabel("Cached-conformer budget")
        ax.set_ylabel(label)
        ax.set_title(label, loc="left", fontweight="bold")
        ax.set_ylim(bottom=0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.08),
        columnspacing=1.5,
        handlelength=2.8,
    )
    fig.suptitle(
        "Conformational diversity across the same allocation budgets",
        fontsize=15,
        fontweight="bold",
        y=1.03,
    )
    fig.text(
        0.5,
        0.955,
        "Symmetry-aware heavy-atom RMSD clustering at 0.5 Å; 20 fixed draws across 22 molecules",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    save_figure(fig, "diversity-by-budget")


def runtime_frame() -> pd.DataFrame:
    rows = []
    for recipe, qwen_fraction in [
        ("qwen", 1.00),
        ("qwen75_loqi25", 0.75),
        ("qwen50_loqi50", 0.50),
        ("qwen25_loqi75", 0.25),
        ("loqi", 0.00),
    ]:
        relative_qwen_cost = qwen_fraction + 2 * (1 - qwen_fraction)
        rows.append(
            {
                "recipe": recipe,
                "qwen_fraction": qwen_fraction,
                "cost_qwen_units": relative_qwen_cost,
                "time_vs_50": 100 * relative_qwen_cost / 1.5,
                "throughput_vs_50": 100 * 1.5 / relative_qwen_cost,
            }
        )
    return pd.DataFrame(rows)


def plot_runtime(runtime: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    labels = [RECIPE_LABELS[recipe] for recipe in runtime["recipe"]]
    colors = [COLORS[recipe] for recipe in runtime["recipe"]]
    positions = np.arange(len(runtime))
    for ax, column, title, ylabel in [
        (
            axes[0],
            "time_vs_50",
            "Estimated generation time",
            "Time relative to 50:50 (%)",
        ),
        (
            axes[1],
            "throughput_vs_50",
            "Estimated equal-time throughput",
            "Candidates relative to 50:50 (%)",
        ),
    ]:
        bars = ax.bar(positions, runtime[column], color=colors, width=0.68)
        ax.axhline(100, color="#17212B", linewidth=1, linestyle="--")
        ax.set_xticks(positions, labels, rotation=25, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_ylim(0, max(runtime[column]) * 1.18)
        for bar, value in zip(bars, runtime[column]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 2.5,
                f"{value:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9,
            )
    fig.suptitle(
        "Runtime sensitivity if Qwen throughput is twice LoQI throughput",
        fontsize=15,
        fontweight="bold",
        y=1.03,
    )
    fig.text(
        0.5,
        0.95,
        "Calculation only; the pharmacophore experiment did not run a controlled timing benchmark",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    save_figure(fig, "runtime-sensitivity")


def pairwise_grid_summary(paired: pd.DataFrame) -> pd.DataFrame:
    selected = paired[primary_mask(paired)].copy()
    rows = []
    for comparator in ["loqi", "rdkit_raw", "rdkit_minimized"]:
        qwen_left = selected[
            (selected["left"] == "qwen") & (selected["right"] == comparator)
        ].copy()
        qwen_right = selected[
            (selected["left"] == comparator) & (selected["right"] == "qwen")
        ].copy()
        if not qwen_right.empty:
            original_low = qwen_right["ci_low"].copy()
            qwen_right["difference"] = -qwen_right["difference"]
            qwen_right["ci_low"] = -qwen_right["ci_high"]
            qwen_right["ci_high"] = -original_low
        subset = pd.concat([qwen_left, qwen_right], ignore_index=True)
        assert len(subset) == 336, (comparator, len(subset))
        difference = subset["difference"].to_numpy()
        equality_tolerance = 1e-8
        rows.append(
            {
                "comparator": comparator,
                "cells": len(subset),
                "qwen_higher": int(np.sum(difference > equality_tolerance)),
                "equal": int(np.sum(np.abs(difference) <= equality_tolerance)),
                "comparator_higher": int(np.sum(difference < -equality_tolerance)),
                "interval_qwen": int(np.sum(subset["ci_low"] > 0)),
                "interval_zero": int(
                    np.sum((subset["ci_low"] <= 0) & (subset["ci_high"] >= 0))
                ),
                "interval_comparator": int(np.sum(subset["ci_high"] < 0)),
            }
        )
    return pd.DataFrame(rows)


def build_report(
    coverage: pd.DataFrame,
    full_coverage: pd.DataFrame,
    paired: pd.DataFrame,
    diversity: pd.DataFrame,
    scope: pd.DataFrame,
    validity: pd.DataFrame,
    runtime: pd.DataFrame,
) -> str:
    primary = coverage.copy()
    strict = primary[
        (primary["condition"] == "6/6") & (primary["radius_angstrom"] == 0.5)
    ].copy()
    strict_scope = strict[(strict["recipe"] == "qwen") & (strict["budget"] == 1000)][
        ["family", "molecules", "parents", "queries"]
    ]
    scope_table = scope.merge(strict_scope, on="family", how="left")
    scope_table["family_order"] = scope_table["family"].map(
        {family: index for index, family in enumerate(FAMILIES)}
    )
    scope_table = scope_table.sort_values("family_order")

    query_rows = []
    for row in scope_table.itertuples(index=False):
        query_rows.append(
            [
                FAMILY_LABELS[row.family],
                f"{row.queries_x:,}",
                f"{row.parents_x:,}",
                str(row.molecules_with_queries),
                f"{row.queries_y:,}",
                str(row.molecules),
            ]
        )

    validity = validity.copy()
    validity["valid_yield"] = 100 * validity["valid_for_coverage"] / validity["supplied"]
    validity_rows = []
    validity_order = ["qwen", "loqi", "rdkit_raw", "rdkit_minimized"]
    for method in validity_order:
        row = validity[validity["method"] == method].iloc[0]
        validity_rows.append(
            [
                RECIPE_LABELS[method],
                f"{int(row.supplied):,}",
                f"{int(row.pb_pass):,}",
                f"{int(row.identity_pass):,}",
                f"{int(row.valid_for_coverage):,}",
                f"{row.valid_yield:.2f}%",
            ]
        )

    diversity_endpoint = diversity[
        (diversity["cluster_radius"] == 0.5) & (diversity["budget"] == 1000)
    ].set_index("recipe")
    diversity_rows = []
    for recipe in RECIPES:
        row = diversity_endpoint.loc[recipe]
        diversity_rows.append(
            [
                RECIPE_LABELS[recipe],
                f"{row.clusters:.2f}",
                f"{row.effective_clusters:.2f}",
                f"{100 * row.largest_fraction:.1f}%",
                f"{100 * row.valid_yield:.2f}%",
            ]
        )

    runtime_rows = []
    for row in runtime.itertuples(index=False):
        runtime_rows.append(
            [
                RECIPE_LABELS[row.recipe],
                f"{100 * row.qwen_fraction:.0f}%",
                f"{100 * (1 - row.qwen_fraction):.0f}%",
                f"{row.cost_qwen_units:.2f}",
                f"{row.time_vs_50:.1f}%",
                f"{row.throughput_vs_50:.1f}%",
            ]
        )

    full_rank = primary.copy()
    scenario_keys = [
        "family",
        "condition",
        "radius_angstrom",
        "construction_count",
        "support_rule",
        "budget",
    ]
    full_rank["rank"] = full_rank.groupby(scenario_keys)["coverage"].rank(
        ascending=False, method="average"
    )
    rank_summary = (
        full_rank.groupby("recipe", as_index=False)
        .agg(mean_coverage=("coverage", "mean"), mean_rank=("rank", "mean"))
        .set_index("recipe")
    )
    rank_rows = []
    for recipe in RECIPES:
        row = rank_summary.loc[recipe]
        rank_rows.append(
            [
                RECIPE_LABELS[recipe],
                "336",
                f"{100 * row.mean_coverage:.2f}%",
                f"{row.mean_rank:.2f}",
            ]
        )

    pair_summary = pairwise_grid_summary(paired)
    pair_rows = []
    for row in pair_summary.itertuples(index=False):
        pair_rows.append(
            [
                RECIPE_LABELS[row.comparator],
                str(row.cells),
                str(row.qwen_higher),
                str(row.equal),
                str(row.comparator_higher),
                str(row.interval_qwen),
                str(row.interval_zero),
                str(row.interval_comparator),
            ]
        )

    secondary = full_coverage[
        full_coverage["family"].str.startswith("shared")
        & (full_coverage["condition"] == "6/6")
        & (full_coverage["radius_angstrom"] == 0.5)
        & (full_coverage["budget"] == 1000)
    ]
    secondary_rows = []
    for (family, count, support, molecules), group in secondary.groupby(
        ["family", "construction_count", "support_rule", "molecules"], sort=False
    ):
        values = group.set_index("recipe")["coverage"] * 100
        support_label = "all poses" if support == "all" else "at least two thirds"
        secondary_rows.append(
            [
                FAMILY_LABELS[family],
                f"{count} poses, {support_label}",
                str(molecules),
                *[f"{values.loc[recipe]:.2f}%" for recipe in RECIPES],
            ]
        )

    breadth_density = full_coverage[
        primary_mask(full_coverage)
        & (full_coverage["condition"] == "4/4")
        & (full_coverage["radius_angstrom"] == 0.5)
        & (full_coverage["budget"] == 1000)
        & full_coverage["recipe"].isin(["qwen", "loqi"])
    ]
    breadth_rows = []
    for family in FAMILIES:
        family_data = breadth_density[breadth_density["family"] == family].set_index(
            "recipe"
        )
        breadth_rows.append(
            [
                FAMILY_LABELS[family],
                f"{100 * family_data.loc['qwen', 'coverage']:.2f}%",
                f"{100 * family_data.loc['loqi', 'coverage']:.2f}%",
                f"{100 * family_data.loc['qwen', 'gated_match_fraction']:.2f}%",
                f"{100 * family_data.loc['loqi', 'gated_match_fraction']:.2f}%",
            ]
        )

    strict_tables = []
    for family in FAMILIES:
        family_data = strict[strict["family"] == family]
        rows = []
        for recipe in RECIPES:
            values = (
                family_data[family_data["recipe"] == recipe]
                .set_index("budget")["coverage"]
                .reindex(BUDGETS)
            )
            rows.append(
                [RECIPE_LABELS[recipe], *[f"{100 * value:.2f}%" for value in values]]
            )
        strict_tables.append(
            f"### {FAMILY_LABELS[family]}\n\n"
            + markdown_table(
                ["Recipe", "Budget 20", "Budget 100", "Budget 500", "Budget 1,000"],
                rows,
            )
        )

    source_coverage = ANALYSIS_TABLES / "coverage_summary.csv"
    source_paired = ANALYSIS_TABLES / "paired_recipe_differences.csv"
    strict_tables_text = "\n\n".join(strict_tables)

    report = f"""# PLINDER-22 pharmacophore allocation experiment

**Internal laboratory report · 23 November 2026**

This experiment compared fixed conformer-library recipes built from Qwen, LoQI,
three Qwen/LoQI mixtures, and RDKit. The main question was whether a generated library
contains at least one conformer that satisfies pharmacophore arrangements observed
in experimental structures. The generators received the ligand only; protein
information was used only afterward to construct two of the evaluation-query families.

The broad result is straightforward. LoQI is strong at small library sizes, Qwen
improves more as the library grows, and the learned generators separate clearly
from RDKit in the stringent large-library comparison. At 1,000 conformers, 50:50
and 75:25 have almost identical coverage. Qwen and 75:25 retain more geometric
diversity than 50:50. If Qwen is assumed to generate accepted conformers twice as
fast as LoQI, 75:25 also has an estimated 16.7% generation-time advantage over
50:50. This timing result is a sensitivity calculation, not a controlled benchmark.

# Experiment design

## How the four query groups were constructed

- **Individual ligand:** Features were read from one experimental ligand pose. No
  protein information was used. The query asks whether a generated conformer can
  reproduce that pose-specific ligand feature arrangement.

- **Individual interaction:** One experimental protein-ligand complex was used.
  Only ligand features involved in the native interaction pattern were retained.
  The query asks whether a ligand-only generated conformer can reproduce those
  protein-selected ligand features.

- **Shared ligand:** Several experimental poses of the same ligand were aligned by
  their ligand features. Features supported across the construction poses were
  retained. A separate PDB entry was held out from construction. The query asks
  whether the library covers a ligand feature arrangement shared across structures.

- **Shared protein:** Binding pockets were aligned using protein C-alpha atoms, and
  interaction features were transferred into the common protein-site frame. They
  were not realigned afterward by ligand features. A separate PDB entry was held
  out. The query asks whether the library covers an interaction arrangement shared
  in the aligned binding site.

![Examples of the four query constructions](figures/query-family-examples.png)

The colored points are pharmacophore features. Gray points show the surrounding
reference geometry. These are examples for explanation, not performance results.

## Query and panel size

{markdown_table(
    [
        "Query family",
        "All queries",
        "Parent hypotheses",
        "Molecules with queries",
        "Strict 6/6 queries",
        "Strict molecules",
    ],
    query_rows,
)}

- The panel contained **22 ligand identities and 1,067 experimental references**.
- The full query construction produced **57,459 unique queries from 809 parent hypotheses**.
- The primary allocation grid contained **2,352 coverage cells**: four families,
  seven feature requirements, three tolerances, seven recipes, and four budgets.
- Secondary shared-query constructions added another **2,352 cells**, producing
  **4,704 coverage cells** in the complete analysis.
- Feature requirements were 3/3, 4/4, 5/5, 6/6, 4/5, 5/6, and the complete available
  hypothesis. Position tolerances were 0.5, 0.75, and 1.0 Å.

## How candidate libraries were tested

- Each source supplied a cache of 1,000 conformers per molecule. No new conformers
  were generated for this analysis.
- Candidate budgets were 20, 100, 500, and 1,000 cached conformers. Mixed libraries
  used exact source quotas. For example, 75:25 at budget 1,000 used 750 Qwen and
  250 LoQI cache entries.
- A query was covered when at least one sampled conformer passed the native CDPKit
  match, PoseBusters, stored-graph identity, and coordinate-derived identity checks.
- Rejected entries still consumed the candidate budget and were not replaced.
- Coverage was averaged within each parent hypothesis, then within each molecule,
  and finally across molecules with equal molecule weight.
- Intervals used 5,000 molecule-level bootstrap resamples. Queries from one molecule
  were therefore not treated as independent observations.
- Query features used positions only. Feature directions and exclusion volumes were
  disabled, and candidate conformers were allowed to contain additional features.

**Coverage in simple terms:** the expected percentage of pharmacophore queries for
which a library of the requested size contains at least one valid matching conformer.

## Cache validity before coverage

{markdown_table(
    ["Source", "Supplied", "PoseBusters pass", "Identity pass", "Final valid", "Valid yield"],
    validity_rows,
)}

The learned sources lost about 9–10% of their supplied entries to the validity and
identity gates. RDKit retained essentially every entry. Since rejected entries were
not replaced, the learned methods achieved their reported coverage with smaller
effective valid pools.

# Coverage results

## Stringent comparison used for the main interpretation

The following comparison requires all six selected features to match within 0.5 Å.
It is deliberately difficult and avoids the strong saturation seen for easier or
looser queries. The tables report means; the figure also shows the molecule-bootstrap
intervals.

![Coverage by family and budget](figures/strict-coverage-by-budget.png)

{strict_tables_text}

The budget pattern is consistent across the four families:

- LoQI is stronger than Qwen at budget 20.
- Qwen gains more from increasing the library and is above LoQI at budget 1,000.
- At budget 1,000, Qwen exceeds the stronger RDKit control by 23.65–32.67 percentage
  points across the four families. LoQI also exceeds the stronger RDKit control by
  18.10–21.02 points.
- The 50:50 and 75:25 mixtures are nearly tied at budget 1,000. The largest difference
  between them is 1.52 points in the shared-protein task; all four paired intervals
  include zero.
- The protein-informed families support the learned-versus-RDKit separation, but the
  ligand-only families remain the direct evidence for a reusable protein-agnostic library.

## Complete primary-grid summary

The stringent tables above are one predeclared slice, not an average over the full
experiment. The next figure averages all seven feature requirements and all three
tolerances within each family, recipe, and budget. This descriptive average is included
to show the complete pattern; it is not used as a separate optimization score.

![Full primary-grid mean coverage](figures/full-grid-mean-coverage.png)

{markdown_table(
    ["Recipe", "Primary-grid cells", "Mean coverage", "Mean rank"],
    rank_rows,
)}

Easy and loose queries often approach 100% coverage, so the full-grid means compress
the differences. Even in this saturated summary, the Qwen-heavy mixtures occupy the
top average ranks, while RDKit raw is lowest overall.

The paired table below uses every primary-grid scenario rather than selected cells.
Intervals are exploratory and are not adjusted for the number of grid comparisons.

{markdown_table(
    [
        "Comparator",
        "Cells",
        "Qwen higher",
        "Equal",
        "Comparator higher",
        "Interval favors Qwen",
        "Interval includes zero",
        "Interval favors comparator",
    ],
    pair_rows,
)}

Qwen is numerically above LoQI in most cells, but most Qwen-LoQI intervals include
zero. Separation from RDKit raw is more consistent. This is why the experiment is
better described as a budget- and task-dependent Qwen-LoQI profile than as a universal
winner-loser result.

## Qwen breadth and LoQI density

Coverage gives the same success score to a query with one matching conformer and a
query with hundreds of matching conformers. The table below uses the 4/4, 0.5 Å,
budget-1,000 setting to compare those two views. The valid-hit fraction is the fraction
of all requested cache slots that are valid matches, averaged with the same molecule
weighting as coverage.

{markdown_table(
    [
        "Family",
        "Qwen coverage",
        "LoQI coverage",
        "Qwen valid-hit fraction",
        "LoQI valid-hit fraction",
    ],
    breadth_rows,
)}

Qwen reaches more distinct queries, while LoQI supplies a larger fraction of matching
conformers for the queries it reaches. This breadth-versus-density difference explains
why LoQI can be efficient at small budgets even though Qwen has greater large-budget
coverage and geometric diversity.

## Complete 50:50 versus 75:25 sensitivity

![Allocation differences across the primary grid](figures/allocation-difference-full-grid.png)

Each box contains all 21 combinations of feature requirement and positional tolerance.
Positive values favor 75:25; negative values favor 50:50. The distributions cross zero
in every family. Across these definitions, the mean difference is smaller than one
percentage point in all 16 family-budget groups. The allocation ordering is therefore
sensitive to the query definition and is not a reliable coverage distinction.

## Secondary shared-query constructions

The primary shared queries used two construction poses with support in both poses.
The complete experiment also tested three construction poses with either support in
all three or support in at least two thirds. The table shows the stringent 6/6,
0.5 Å endpoint at budget 1,000.

{markdown_table(
    [
        "Family",
        "Construction",
        "N",
        *[RECIPE_LABELS[recipe] for recipe in RECIPES],
    ],
    secondary_rows,
)}

The Qwen-heavy recipes remain strong in the secondary constructions. The shared-protein
rows have only 6–11 eligible molecules, so small differences should not be treated as
stable allocation rankings.

# Geometric diversity

Coverage asks whether at least one conformer matches. Diversity asks how many distinct
geometries the library contains. Valid conformers were clustered by symmetry-aware
heavy-atom RMSD. The analysis used 20 fixed draws per allocation and budget across all
22 molecules.

![Diversity by budget](figures/diversity-by-budget.png)

At budget 1,000 and a 0.5 Å clustering radius:

{markdown_table(
    ["Recipe", "Mean clusters", "Effective clusters", "Largest cluster fraction", "Valid yield"],
    diversity_rows,
)}

- Pure Qwen has the highest diversity at every budget among these allocation recipes.
- Increasing the LoQI fraction progressively reduces the number of distinct clusters.
- At budget 1,000, 75:25 has 256.74 mean clusters, compared with 220.46 for 50:50.
  This is a 16.5% increase over 50:50.
- RDKit raw is more diverse than LoQI and RDKit minimized, but its pharmacophore
  coverage is substantially lower at the stringent large-library endpoint.
- Effective cluster counts give the same overall ordering and show that repeated
  occupancy within common clusters is strongest for LoQI-heavy libraries.

# Generation-time sensitivity

Generation time was not measured by this cached-conformer experiment. The following
calculation uses the laboratory assumption that **Qwen generates accepted conformers
at twice the LoQI rate**. One Qwen conformer is assigned one time unit and one LoQI
conformer two time units. The calculation represents accumulated generation work under
comparable conditions; exact elapsed wall time can differ with batching, memory, and
parallel execution.

![Runtime sensitivity](figures/runtime-sensitivity.png)

{markdown_table(
    [
        "Recipe",
        "Qwen share",
        "LoQI share",
        "Relative cost, Qwen = 1",
        "Time versus 50:50",
        "Equal-time throughput versus 50:50",
    ],
    runtime_rows,
)}

Under this assumption, a 1,000-conformer 75:25 library costs 1,250 relative time units,
whereas 50:50 costs 1,500. The 75:25 recipe therefore uses **16.7% less accumulated
generation time** or produces **20% more candidates in the same generation time**.
RDKit is excluded from this estimate because its CPU implementation and cost model are
not comparable to the Qwen-LoQI assumption.

# Combined interpretation for internal discussion

1. **The strongest result is learned versus RDKit at scale.** Qwen, LoQI, and their
   mixtures cover substantially more stringent queries than RDKit at budget 1,000,
   despite the learned pools having lower valid yield.

2. **Qwen and LoQI contribute different sampling profiles.** LoQI is efficient at
   small budgets and concentrates more samples in covered regions. Qwen expands to
   broader pharmacophore and geometric coverage as the budget increases.

3. **Coverage alone does not distinguish 50:50 from 75:25 at the production budget.**
   Their budget-1,000 ligand-only means differ by 0.26 and 0.01 percentage points,
   and both paired intervals include zero.

4. **Diversity favors the Qwen-heavy allocation.** The 75:25 recipe retains about
   36 more mean clusters per molecule than 50:50 at budget 1,000.

5. **The twofold Qwen-speed assumption also favors 75:25.** It gives an estimated
   16.7% time saving relative to 50:50. This direction is clear, but the magnitude
   should remain labeled as an estimate until both generators are timed under one
   production protocol.

6. **Practical internal recommendation:** if the twofold speed difference is confirmed
   on the intended production hardware, 75:25 is the stronger operational choice for
   a 1,000-conformer release. It has indistinguishable coverage, higher diversity, and
   lower estimated generation cost. A 50:50 library remains a simple equal-source
   baseline, but the present evidence does not make it the best-performing production
   allocation.

# Scope and limitations

- The study contains 22 selected ligand identities and is exploratory.
- Coverage means proximity to the tested pharmacophore arrangements. It does not
  demonstrate binding affinity, docking enrichment, or downstream model improvement.
- Protein-informed queries evaluate binding-site relevance after ligand-only generation;
  they do not make the generators protein-conditioned.
- The same molecule contributes many related queries. Molecules, not queries, are the
  uncertainty unit.
- Easy queries and larger tolerances saturate. The stringent 6/6, 0.5 Å condition is
  more discriminating but has fewer eligible molecules.
- The runtime section is an explicit sensitivity calculation. It is not a measured
  hardware-matched benchmark.

# Data and reproduction

- Compact paper archive: `../../results-data/cdpkit-pharmacophore-mixture/`
- Complete analysis tables: `{ANALYSIS_TABLES}`
- Complete coverage grid SHA-256: `{sha256(source_coverage)}`
- Paired-difference grid SHA-256: `{sha256(source_paired)}`
- Rebuild this report and all figures with `python build_report.py` from this directory.
"""
    return report


def main() -> None:
    use_report_style()
    coverage = pd.read_csv(ARCHIVE_DIR / "coverage-primary.csv")
    diversity = pd.read_csv(ARCHIVE_DIR / "diversity-summary.csv")
    scope = pd.read_csv(ARCHIVE_DIR / "query-scope.csv")
    validity = pd.read_csv(ARCHIVE_DIR / "validity-method-counts.csv")
    full_coverage = pd.read_csv(ANALYSIS_TABLES / "coverage_summary.csv")
    paired = pd.read_csv(ANALYSIS_TABLES / "paired_recipe_differences.csv")

    assert len(coverage) == 2352
    assert len(full_coverage) == 4704
    assert int(scope["queries"].sum()) == 57459
    assert int(scope["parents"].sum()) == 809
    assert set(coverage["recipe"]) == set(RECIPES)

    primary = coverage.copy()
    strict = primary[
        (primary["condition"] == "6/6") & (primary["radius_angstrom"] == 0.5)
    ]
    assert len(strict) == 112

    plot_strict_coverage(strict)
    plot_full_grid_heatmap(primary)
    plot_allocation_difference(primary)
    plot_diversity(diversity)
    runtime = runtime_frame()
    plot_runtime(runtime)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        ANALYSIS_DIR / "figures" / "query_family_examples.png",
        FIGURE_DIR / "query-family-examples.png",
    )

    report = build_report(
        coverage,
        full_coverage,
        paired,
        diversity,
        scope,
        validity,
        runtime,
    )
    (REPORT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
