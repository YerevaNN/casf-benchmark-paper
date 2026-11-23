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
# Put ligand-only tasks first because they are the primary evidence for selecting
# a reusable, protein-agnostic library.
FAMILIES = [
    "individual_ligand",
    "shared_ligand",
    "individual_interaction",
    "shared_protein",
]
FAMILY_LABELS = {
    "individual_ligand": "Individual ligand",
    "shared_ligand": "Shared ligand",
    "individual_interaction": "Individual interaction",
    "shared_protein": "Shared protein",
}
FAMILY_INFORMATION = {
    "individual_ligand": "One ligand pose",
    "shared_ligand": "Two or three ligand poses",
    "individual_interaction": "One protein–ligand complex",
    "shared_protein": "Two or three aligned protein sites",
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
LEARNED_RECIPES = RECIPES[:5]
RDKIT_RECIPES = RECIPES[-2:]
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
    "rdkit_raw": "#7B838A",
    "rdkit_minimized": "#343A40",
}
LINESTYLES = {
    "qwen": "-",
    "loqi": "-",
    "qwen50_loqi50": "-",
    "qwen75_loqi25": "-",
    "qwen25_loqi75": "-.",
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
    svg_path = FIGURE_DIR / f"{stem}.svg"
    fig.savefig(svg_path, bbox_inches="tight")
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
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
                linewidth=1.6 if recipe == "qwen25_loqi75" else 2.2,
                alpha=0.75 if recipe == "qwen25_loqi75" else 1.0,
                color=COLORS[recipe],
                linestyle=LINESTYLES[recipe],
                label=RECIPE_LABELS[recipe],
            )
        role = (
            "ligand-only selection task"
            if family in {"individual_ligand", "shared_ligand"}
            else "protein-informed check"
        )
        ax.set_title(f"{FAMILY_LABELS[family]}\n{role}", loc="left", fontweight="bold")
        ax.set_xscale("log")
        ax.set_xticks(BUDGETS)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_ylim(0, 102)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_ylabel("Queries covered (%)")
        ax.set_xlabel("Cached candidates allowed")
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
        "Hard task: all six features must match within 0.5 Å",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.965,
        "Each panel is a separate task; values are not averaged across task families",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
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
        ax.set_xlabel("Cached candidates allowed")
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
    colorbar.set_label("Descriptive mean coverage (%)")
    fig.suptitle("Full-grid sensitivity check", fontsize=15, fontweight="bold", y=1.01)
    fig.text(
        0.5,
        0.965,
        "Each cell averages seven feature rules and three distance tolerances",
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


def plot_diversity(diversity: pd.DataFrame) -> None:
    selected = diversity[diversity["cluster_radius"] == 0.5]
    fig, ax = plt.subplots(figsize=(9.6, 5.5))
    for recipe in RECIPES:
        values = selected[selected["recipe"] == recipe].sort_values("budget")
        ax.plot(
            values["budget"],
            values["clusters"],
            marker="o",
            markersize=4.5,
            linewidth=1.6 if recipe == "qwen25_loqi75" else 2.2,
            alpha=0.75 if recipe == "qwen25_loqi75" else 1.0,
            color=COLORS[recipe],
            linestyle=LINESTYLES[recipe],
            label=RECIPE_LABELS[recipe],
        )
    ax.set_xscale("log")
    ax.set_xticks(BUDGETS)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xlabel("Cached candidates allowed")
    ax.set_ylabel("Mean number of distinct geometry clusters")
    ax.set_ylim(bottom=0)
    ax.legend(ncol=2, loc="upper left")
    fig.suptitle(
        "Geometric diversity at a 0.5 Å clustering radius",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.945,
        "Symmetry-aware heavy-atom RMSD; 20 fixed draws across 22 molecules",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
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
        relative_cost = qwen_fraction + 2 * (1 - qwen_fraction)
        rows.append(
            {
                "recipe": recipe,
                "time_vs_50": 100 * relative_cost / 1.5,
                "throughput_vs_50": 100 * 1.5 / relative_cost,
            }
        )
    return pd.DataFrame(rows)


def plot_runtime(runtime: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    labels = [RECIPE_LABELS[recipe] for recipe in runtime["recipe"]]
    colors = [COLORS[recipe] for recipe in runtime["recipe"]]
    positions = np.arange(len(runtime))
    bars = ax.bar(positions, runtime["time_vs_50"], color=colors, width=0.68)
    ax.axhline(100, color="#17212B", linewidth=1, linestyle="--")
    ax.set_xticks(positions, labels, rotation=20, ha="right")
    ax.set_ylabel("Estimated time relative to 50:50 (%)")
    ax.set_ylim(0, max(runtime["time_vs_50"]) * 1.15)
    for bar, value in zip(bars, runtime["time_vs_50"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 2.5,
            f"{value:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    fig.suptitle(
        "Simple timing estimate if Qwen is twice as fast as LoQI",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        0.94,
        "Calculation only; no controlled generator timing was performed in this experiment",
        ha="center",
        fontsize=10,
        color="#4F5B66",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save_figure(fig, "runtime-sensitivity")


def hard_coverage_rows(strict: pd.DataFrame, budget: int) -> list[list[str]]:
    selected = strict[strict["budget"] == budget].pivot(
        index="recipe", columns="family", values="coverage"
    )
    rows = []
    for recipe in RECIPES:
        rows.append(
            [
                RECIPE_LABELS[recipe],
                *[f"{100 * selected.loc[recipe, family]:.2f}%" for family in FAMILIES],
            ]
        )
    return rows


def signed_points(value: float) -> str:
    return f"{value:+.2f} points"


def build_report(
    coverage: pd.DataFrame,
    full_coverage: pd.DataFrame,
    diversity: pd.DataFrame,
    scope: pd.DataFrame,
    validity: pd.DataFrame,
    runtime: pd.DataFrame,
) -> str:
    strict = coverage[
        (coverage["condition"] == "6/6") & (coverage["radius_angstrom"] == 0.5)
    ].copy()

    strict_scope = strict[(strict["recipe"] == "qwen") & (strict["budget"] == 1000)][
        ["family", "molecules", "queries"]
    ].rename(columns={"molecules": "strict_molecules", "queries": "strict_queries"})
    scope_table = scope.merge(strict_scope, on="family", how="left").set_index("family")
    query_rows = []
    for family in FAMILIES:
        row = scope_table.loc[family]
        query_rows.append(
            [
                FAMILY_LABELS[family],
                FAMILY_INFORMATION[family],
                f"{int(row.queries):,}",
                f"{int(row.strict_queries):,}",
                str(int(row.strict_molecules)),
            ]
        )

    validity = validity.copy()
    validity["valid_yield"] = 100 * validity["valid_for_coverage"] / validity["supplied"]
    validity_rows = []
    for method in ["qwen", "loqi", "rdkit_raw", "rdkit_minimized"]:
        row = validity[validity["method"] == method].iloc[0]
        validity_rows.append(
            [
                RECIPE_LABELS[method],
                f"{int(row.valid_for_coverage):,} of {int(row.supplied):,}",
                f"{row.valid_yield:.2f}%",
            ]
        )

    strict_1000 = strict[strict["budget"] == 1000].pivot(
        index="recipe", columns="family", values="coverage"
    )
    direct_rows = []
    for family in FAMILIES:
        learned = strict_1000.loc[LEARNED_RECIPES, family]
        rdkit = strict_1000.loc[RDKIT_RECIPES, family]
        best_learned = learned.idxmax()
        best_rdkit = rdkit.idxmax()
        direct_rows.append(
            [
                FAMILY_LABELS[family],
                signed_points(
                    100
                    * (
                        strict_1000.loc["qwen", family]
                        - strict_1000.loc["loqi", family]
                    )
                ),
                signed_points(
                    100
                    * (
                        strict_1000.loc["qwen75_loqi25", family]
                        - strict_1000.loc["qwen50_loqi50", family]
                    )
                ),
                f"{RECIPE_LABELS[best_learned]} ({100 * learned.max():.2f}%)",
                f"{RECIPE_LABELS[best_rdkit]} ({100 * rdkit.max():.2f}%)",
                signed_points(100 * (learned.max() - rdkit.max())),
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

    secondary = full_coverage[
        full_coverage["family"].str.startswith("shared")
        & (full_coverage["construction_count"] == 3)
        & (full_coverage["condition"] == "6/6")
        & (full_coverage["radius_angstrom"] == 0.5)
        & (full_coverage["budget"] == 1000)
    ]
    secondary_rows = []
    for (family, support, molecules), group in secondary.groupby(
        ["family", "support_rule", "molecules"], sort=False
    ):
        values = group.set_index("recipe")["coverage"] * 100
        best_rdkit = values.loc[RDKIT_RECIPES].idxmax()
        secondary_rows.append(
            [
                FAMILY_LABELS[family],
                "All three poses" if support == "all" else "At least two of three",
                str(molecules),
                f"{values.loc['qwen']:.2f}%",
                f"{values.loc['loqi']:.2f}%",
                f"{values.loc['qwen50_loqi50']:.2f}%",
                f"{values.loc['qwen75_loqi25']:.2f}%",
                f"{RECIPE_LABELS[best_rdkit]} ({values.loc[best_rdkit]:.2f}%)",
            ]
        )

    diversity_endpoint = diversity[
        (diversity["cluster_radius"] == 0.5) & (diversity["budget"] == 1000)
    ].set_index("recipe")
    baseline_clusters = diversity_endpoint.loc["qwen50_loqi50", "clusters"]
    diversity_rows = []
    for recipe in RECIPES:
        clusters = diversity_endpoint.loc[recipe, "clusters"]
        diversity_rows.append(
            [
                RECIPE_LABELS[recipe],
                f"{clusters:.2f}",
                f"{clusters - baseline_clusters:+.2f}",
            ]
        )

    runtime_rows = []
    for row in runtime.itertuples(index=False):
        runtime_rows.append(
            [
                RECIPE_LABELS[row.recipe],
                f"{row.time_vs_50:.1f}%",
                f"{row.throughput_vs_50:.1f}%",
            ]
        )

    hard_headers = [
        "Recipe",
        "Individual ligand",
        "Shared ligand",
        "Individual interaction",
        "Shared protein",
    ]
    source_coverage = ANALYSIS_TABLES / "coverage_summary.csv"

    report = f"""# PLINDER-22 pharmacophore allocation experiment

**Internal laboratory report · 23 November 2026**

This study asks a practical question: when we build a reusable conformer library,
which pure source or fixed Qwen/LoQI allocation should fill its candidate slots?
RDKit provides the conventional controls. The experiment
does **not** confirm that 50:50, 75:25, or any other mixture is the single best policy.
It shows which conclusions are clear, which depend on library size, and where the
two main mixtures remain almost tied.

# Start with the query hypothesis

A pharmacophore query is a **small 3D checklist**. It says that particular chemical
features—for example, an aromatic center, one donor, and two acceptors—should occur
at particular positions. A generated conformer covers the query if it places every
requested feature close enough to those positions.

Imagine that the same ligand has experimental poses A, B, and C:

- **Individual ligand:** take the feature arrangement from pose A. Example: can a
  generated conformer reproduce A's ring center, donor, and acceptors?
- **Shared ligand:** align poses A and B by the ligand and keep the features found
  in both. Example: can a conformer reproduce the ligand geometry that remains
  common across the two structures?
- **Individual interaction:** start from pose A in its protein complex and keep only
  ligand features involved in contacts. Example: if one acceptor in A does not contact
  the protein, it is left out of this query.
- **Shared protein:** align the protein sites for A and B and keep contact features
  supported in the common site frame. Example: can a conformer place its donor and
  acceptors where both aligned binding sites use them?

The first two families are ligand-only evaluations and are the primary evidence for
selecting a protein-agnostic library. The last two use protein information only to
evaluate the already generated conformers; they are useful binding-site relevance
checks, but they should not by themselves determine a protein-agnostic generation
policy.

![Examples of the four query constructions](figures/query-family-examples.png)

The colored points are the feature positions in the query. The molecular and protein
geometry around them is shown only to explain where the hypotheses came from.

# What selection and coverage mean

## Selection

Selection here means **assigning library slots to a source before testing the
queries**. It does not mean examining all conformers and keeping the best matches.

- Qwen at budget 1,000 uses 1,000 cached Qwen candidates.
- LoQI at budget 1,000 uses 1,000 cached LoQI candidates.
- 50:50 uses 500 Qwen and 500 LoQI candidates.
- 75:25 uses 750 Qwen and 250 LoQI candidates.
- 25:75 was also run as an unequal-allocation sensitivity check.
- Each RDKit control uses 1,000 candidates from its own cache.

Rejected candidates still consume their assigned slots. The budget is therefore a
**cached-candidate budget**, not the number of valid conformers and not the number of
generation attempts.

## Coverage

One query is covered when at least one candidate:

1. is chemically valid and has the correct molecular identity; and
2. matches every feature required by that query within the chosen distance.

If a library covers 80% of a set of queries, an intuitive reading is that it satisfies
about 8 of every 10 3D checklists. The formal calculation first averages related
queries within their parent hypothesis, then within each molecule, and finally gives
each molecule equal weight.

The main hard task in this report is **6/6 at 0.5 Å**: all six selected features must
match, and every feature must be within half an ångström of its target position.

# Experiment setup

## Query scale

{markdown_table(
    [
        "Query family",
        "Information used",
        "Total queries",
        "Hard 6/6 queries",
        "Molecules in hard task",
    ],
    query_rows,
)}

- The panel contained 22 ligand identities and 1,067 experimental references.
- Query construction produced 57,459 unique queries from 809 parent hypotheses.
- The primary grid tested seven feature rules, three distance tolerances, seven
  recipes, and four candidate budgets: 20, 100, 500, and 1,000.
- Additional three-pose constructions were tested for the two shared-query families.

## Candidate validity

{markdown_table(
    ["Source", "Final valid candidates", "Valid share"],
    validity_rows,
)}

The learned sources lost about 9–10% of their cache entries to the validity and
identity checks, while RDKit retained essentially all entries. Failed entries were
not replaced, so the learned methods reached their reported coverage with fewer
valid candidates inside the same nominal budget.

# Results

## Hard-task coverage at each library size

The figure keeps the four tasks separate. The top row contains the two ligand-only
tasks used for selection; the bottom row contains the protein-informed checks.

![Hard-task coverage by family and budget](figures/strict-coverage-by-budget.png)

Budget 20 shows small-library behavior, and budget 1,000 shows the intended
large-library comparison. The figure shows the intermediate budgets.

### Budget 20

{markdown_table(hard_headers, hard_coverage_rows(strict, 20))}

### Budget 1,000

{markdown_table(hard_headers, hard_coverage_rows(strict, 1000))}

At budget 20, LoQI and the LoQI-heavy mixtures are usually strongest. As the library
grows, Qwen gains more coverage. At budget 1,000, Qwen is above LoQI in every hard
task, and every learned recipe is clearly above the stronger RDKit control.

## Direct hard-task differences at budget 1,000

Each row below is one task under the same 6/6, 0.5 Å rule. There is **no averaging
across the four task families** in this table. Positive differences favor Qwen over
LoQI or 75:25 over 50:50.

{markdown_table(
    [
        "Hard task",
        "Qwen − LoQI",
        "75:25 − 50:50",
        "Numerically highest learned recipe",
        "Stronger RDKit control",
        "Learned − RDKit",
    ],
    direct_rows,
)}

This table supports three direct observations:

- The learned-versus-RDKit advantage is large in every hard task: about 25–33 points
  for the numerically highest learned recipe.
- Qwen's large-budget advantage over LoQI is task dependent: it ranges from 2.85
  points for shared ligand to 12.21 points for shared protein.
- The 50:50 and 75:25 coverage values are almost the same in the two ligand-only
  tasks. Their differences are +0.26 points for individual ligand and −0.01 points
  for shared ligand. These numbers do not establish either mixture as the structural
  winner.

## Full-grid sensitivity check

The main comparison above uses the deliberately difficult 6/6, 0.5 Å task. For
completeness, the next figure summarizes the entire primary grid. Each cell averages
the seven feature rules and three distance tolerances for one task, recipe, and
budget. This descriptive average is a check against relying on one hard setting; it
is not the basis for choosing an allocation.

![Full primary-grid sensitivity summary](figures/full-grid-mean-coverage.png)

The easier rules often approach complete coverage, which compresses differences.
The broad pattern remains the same: learned libraries separate from RDKit as the
budget grows, and no single Qwen/LoQI allocation dominates every task and budget.

## Qwen breadth and LoQI density

Coverage records whether a query is reached at least once. It does not say whether
one or hundreds of candidates match that query. The table below uses the 4/4,
0.5 Å, budget-1,000 setting to show both ideas. “Matching candidate share” is the
percentage of all requested cache slots that are both valid and matching.

{markdown_table(
    [
        "Task",
        "Qwen coverage",
        "LoQI coverage",
        "Qwen matching candidate share",
        "LoQI matching candidate share",
    ],
    breadth_rows,
)}

Qwen reaches more different queries at this large budget, while LoQI places a larger
share of its candidates in regions that it already reaches. In simple terms, Qwen
has greater breadth and LoQI has greater local density. That difference is a
reasonable motivation for testing mixtures, but it does not by itself prove that
one mixture ratio is best.

## Additional shared-query constructions

The primary shared queries were built from two poses with support required in both.
The complete experiment also built queries from three poses. The table below keeps
the hard 6/6, 0.5 Å, budget-1,000 setting and reports the main decision candidates.
The 25:75 sensitivity results remain in the complete data archive.

{markdown_table(
    [
        "Task",
        "Three-pose rule",
        "Molecules",
        "Qwen",
        "LoQI",
        "50:50",
        "75:25",
        "Stronger RDKit control",
    ],
    secondary_rows,
)}

The additional constructions preserve the large learned-versus-RDKit separation.
The shared-protein rows contain only 6–9 eligible molecules, so their small differences
between learned recipes should be treated as exploratory.

# Geometric diversity

Coverage asks whether at least one candidate passes a checklist. Diversity asks how
many distinct 3D geometries are present. Valid conformers were grouped by
symmetry-aware heavy-atom RMSD at a 0.5 Å clustering radius. More clusters mean that
the library explores more distinct shapes.

![Geometric diversity by candidate budget](figures/diversity-by-budget.png)

At budget 1,000:

{markdown_table(
    ["Recipe", "Mean distinct clusters", "Difference from 50:50"],
    diversity_rows,
)}

Among the tested allocation recipes, pure Qwen has the most clusters. The 75:25
mixture retains 256.74 clusters on average, compared with 220.46 for 50:50 and
100.68 for LoQI. Thus, 75:25 gives up some of Qwen's breadth but preserves more of
it than 50:50.

# Simple generation-time estimate

Generation time was not measured in this cached-conformer experiment. The following
calculation uses one explicit assumption: **one cached Qwen candidate takes one
time unit, while one cached LoQI candidate takes two time units**. This represents
the laboratory expectation that Qwen is about twice as fast. Actual wall time can
change with batching, GPU memory, and parallel execution.

![Estimated generation time by allocation](figures/runtime-sensitivity.png)

{markdown_table(
    [
        "Recipe",
        "Estimated time versus 50:50",
        "Candidates in the same time versus 50:50",
    ],
    runtime_rows,
)}

For 1,000 candidates, 50:50 costs 1,500 relative time units and 75:25 costs 1,250.
Under the twofold-speed assumption, 75:25 therefore uses 16.7% less generation time,
or produces 20% more candidates in the same time. This gives 75:25 a practical edge
over 50:50 **if the speed assumption holds on the production hardware**. It is not
evidence that 75:25 has better pharmacophore coverage.

# Selection interpretation

The ligand-only hard tasks are the appropriate starting point for a protein-agnostic
library decision:

- At budget 1,000, 75:25 and 50:50 differ by only +0.26 points on individual ligand
  and −0.01 points on shared ligand.
- Their hard-task coverage is therefore effectively tied at the precision useful for
  this exploratory panel.
- Diversity favors 75:25 over 50:50.
- The simple timing estimate also favors 75:25 if Qwen is truly twice as fast.
- 50:50 remains an easy-to-explain equal-source baseline.

The result is a choice between two reasonable policies, not confirmation of one
optimum. **50:50 is the simpler balanced baseline. 75:25 is the more attractive
operational candidate when diversity and the current timing assumption are included.**
A controlled timing run and a larger molecular panel would be needed before calling
either allocation the confirmed best policy.

# Conclusions

The clearest result is that learned conformer libraries outperform both RDKit controls
on the difficult pharmacophore tasks at large scale. Qwen and LoQI also show different
sampling behavior: Qwen expands to broader coverage and greater geometric diversity,
whereas LoQI supplies denser matching populations and stronger small-budget coverage.
This complementarity makes a mixture reasonable. The present experiment does not
separate 50:50 from 75:25 by hard-task coverage, but diversity and the assumed Qwen
speed give 75:25 a practical—not yet confirmed—advantage.

# Scope and data

- The study is exploratory and contains 22 selected ligand identities.
- Coverage measures agreement with the tested 3D feature arrangements. It does not
  demonstrate binding affinity, docking enrichment, or downstream model improvement.
- Protein-informed queries are after-the-fact relevance checks; the generators still
  receive ligand information only.
- The runtime result is a calculation under one speed assumption, not a measured
  hardware-matched benchmark.
- Compact paper archive: ../../results-data/cdpkit-pharmacophore-mixture/
- Complete analysis tables: {ANALYSIS_TABLES}
- Complete coverage grid SHA-256: {sha256(source_coverage)}
- Rebuild the report and figures by running python build_report.py in this directory.
"""
    return report


def main() -> None:
    use_report_style()
    coverage = pd.read_csv(ARCHIVE_DIR / "coverage-primary.csv")
    diversity = pd.read_csv(ARCHIVE_DIR / "diversity-summary.csv")
    scope = pd.read_csv(ARCHIVE_DIR / "query-scope.csv")
    validity = pd.read_csv(ARCHIVE_DIR / "validity-method-counts.csv")
    full_coverage = pd.read_csv(ANALYSIS_TABLES / "coverage_summary.csv")

    assert len(coverage) == 2352
    assert len(full_coverage) == 4704
    assert int(scope["queries"].sum()) == 57459
    assert int(scope["parents"].sum()) == 809
    assert set(coverage["recipe"]) == set(RECIPES)

    strict = coverage[
        (coverage["condition"] == "6/6") & (coverage["radius_angstrom"] == 0.5)
    ]
    assert len(strict) == 112

    plot_strict_coverage(strict)
    plot_full_grid_heatmap(coverage)
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
        diversity,
        scope,
        validity,
        runtime,
    )
    (REPORT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
