#!/usr/bin/env python3
"""Render the CDPKit pharmacophore allocation figure and manuscript tables."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIGURE_BASENAME = ROOT / "figures" / "figure-pharmacophore-mixture"
SUPPORTING = ROOT / "figures" / "supporting"

REFERENCE_CONDITION = "6/6"
REFERENCE_RADIUS = 0.5
RECIPES = ["qwen", "qwen75_loqi25", "qwen50_loqi50", "qwen25_loqi75", "loqi"]
MAIN_RECIPES = ["qwen", "loqi", "qwen50_loqi50", "rdkit_raw", "rdkit_minimized"]
FAMILIES = [
    "individual_ligand",
    "individual_interaction",
    "shared_ligand",
    "shared_protein",
]
LABELS = {
    "qwen": "Qwen",
    "qwen75_loqi25": "Qwen 75 / LoQI 25",
    "qwen50_loqi50": "Qwen 50 / LoQI 50",
    "qwen25_loqi75": "Qwen 25 / LoQI 75",
    "loqi": "LoQI",
    "rdkit_raw": "RDKit raw",
    "rdkit_minimized": "RDKit minimized",
}
SHORT_LABELS = {
    "qwen": "Qwen",
    "qwen75_loqi25": "75/25",
    "qwen50_loqi50": "50/50",
    "qwen25_loqi75": "25/75",
    "loqi": "LoQI",
    "rdkit_raw": "RDKit raw",
    "rdkit_minimized": "RDKit min.",
}
FAMILY_LABELS = {
    "individual_ligand": "Individual ligand",
    "individual_interaction": "Individual interaction",
    "shared_ligand": "Shared ligand",
    "shared_protein": "Shared protein",
}
COLORS = {
    "qwen": "#0072B2",
    "qwen75_loqi25": "#56B4E9",
    "qwen50_loqi50": "#7B3294",
    "qwen25_loqi75": "#D55E00",
    "loqi": "#E69F00",
    "rdkit_raw": "#6F7A83",
    "rdkit_minimized": "#A7AFB5",
}
MARKERS = {
    "qwen": "o",
    "loqi": "s",
    "qwen50_loqi50": "D",
    "rdkit_raw": "^",
    "rdkit_minimized": "v",
}
LINESTYLES = {
    "qwen": "-",
    "loqi": "-",
    "qwen50_loqi50": "-",
    "rdkit_raw": "--",
    "rdkit_minimized": ":",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_coverage(frame: pd.DataFrame) -> pd.DataFrame:
    return frame[
        (frame["condition"] == REFERENCE_CONDITION)
        & np.isclose(frame["radius_angstrom"], REFERENCE_RADIUS)
    ].copy()


def style_axis(axis: plt.Axes) -> None:
    axis.grid(color="#E7EBEE", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def plot_coverage(
    axis: plt.Axes,
    frame: pd.DataFrame,
    family: str,
    letter: str,
    *,
    secondary: bool = False,
) -> None:
    subset = frame[frame["family"] == family]
    for recipe in MAIN_RECIPES:
        curve = subset[subset["recipe"] == recipe].sort_values("budget")
        if len(curve) != 4:
            raise ValueError(f"Expected four budgets for {family}/{recipe}")
        highlighted = recipe == "qwen50_loqi50"
        axis.plot(
            curve["budget"],
            100 * curve["coverage"],
            color=COLORS[recipe],
            marker=MARKERS[recipe],
            linestyle=LINESTYLES[recipe],
            markersize=5.2 if highlighted else 4.2,
            linewidth=2.8 if highlighted else 1.7,
            alpha=1.0 if highlighted else 0.92,
            zorder=4 if highlighted else 2,
        )
        axis.fill_between(
            curve["budget"],
            100 * curve["ci_low"],
            100 * curve["ci_high"],
            color=COLORS[recipe],
            alpha=0.045 if secondary else 0.07,
            linewidth=0,
        )
    if secondary:
        axis.set_facecolor("#F3F5F6")
    axis.set_xscale("log")
    axis.set_xlim(17, 1180)
    axis.set_xticks([20, 100, 500, 1000])
    axis.set_xticklabels(["20", "100", "500", "1,000"])
    axis.set_ylim(18, 101)
    axis.set_yticks([20, 40, 60, 80, 100])
    axis.set_xlabel("Attempt budget")
    axis.set_ylabel("Expected query coverage (%)")
    axis.set_title(f"{letter}   {FAMILY_LABELS[family]}", loc="left", pad=9, fontweight="bold")
    if secondary:
        axis.text(
            0.97,
            0.95,
            "SECONDARY VALIDATION",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=6.7,
            fontweight="bold",
            color="#66717A",
        )
    style_axis(axis)


def normalized_difference(row: pd.Series) -> tuple[float, float, float]:
    if row["left"] == "qwen50_loqi50":
        return row["difference"], row["ci_low"], row["ci_high"]
    if row["right"] == "qwen50_loqi50":
        return -row["difference"], -row["ci_high"], -row["ci_low"]
    raise ValueError("Contrast does not contain the 50/50 recipe")


def plot_differences(axis: plt.Axes, frame: pd.DataFrame) -> None:
    comparator_order = ["qwen", "qwen75_loqi25", "loqi"]
    offsets = {"qwen": -0.22, "qwen75_loqi25": 0.0, "loqi": 0.22}
    markers = {"qwen": "o", "qwen75_loqi25": "s", "loqi": "^"}
    positions = np.arange(len(FAMILIES))
    for comparator in comparator_order:
        values = []
        lows = []
        highs = []
        for family in FAMILIES:
            rows = frame[
                (frame["family"] == family)
                & (
                    ((frame["left"] == "qwen50_loqi50") & (frame["right"] == comparator))
                    | ((frame["right"] == "qwen50_loqi50") & (frame["left"] == comparator))
                )
            ]
            if len(rows) != 1:
                raise ValueError(f"Expected one 50/50 contrast for {family}/{comparator}")
            value, low, high = normalized_difference(rows.iloc[0])
            values.append(100 * value)
            lows.append(100 * low)
            highs.append(100 * high)
        values_array = np.asarray(values)
        axis.errorbar(
            values_array,
            positions + offsets[comparator],
            xerr=np.vstack([values_array - lows, np.asarray(highs) - values_array]),
            fmt=markers[comparator],
            markersize=5,
            color=COLORS[comparator],
            ecolor=COLORS[comparator],
            elinewidth=1.2,
            capsize=2,
            label=f"50/50 minus {SHORT_LABELS[comparator]}",
        )
    axis.axvline(0, color="#66717C", linewidth=0.9, linestyle="--")
    axis.set_yticks(positions, [FAMILY_LABELS[family] for family in FAMILIES])
    axis.invert_yaxis()
    axis.set_xlim(-1.7, 11.5)
    axis.set_xlabel("Coverage difference at budget 20 (percentage points)\n95% molecule-bootstrap interval")
    axis.set_title("E   Paired small-budget differences", loc="left", pad=9, fontweight="bold")
    style_axis(axis)


def plot_diversity(axis: plt.Axes, frame: pd.DataFrame) -> None:
    subset = frame[np.isclose(frame["cluster_radius"], REFERENCE_RADIUS)]
    for recipe in RECIPES:
        curve = subset[subset["recipe"] == recipe].sort_values("budget")
        if len(curve) != 4:
            raise ValueError(f"Expected four diversity budgets for {recipe}")
        highlighted = recipe == "qwen50_loqi50"
        axis.plot(
            curve["budget"],
            curve["effective_clusters"],
            color=COLORS[recipe],
            marker="o",
            markersize=4.5 if highlighted else 3.5,
            linewidth=2.7 if highlighted else 1.35,
            alpha=1.0 if highlighted else 0.85,
            zorder=4 if highlighted else 2,
        )
    axis.set_xscale("log")
    axis.set_xlim(17, 1180)
    axis.set_xticks([20, 100, 500, 1000])
    axis.set_xticklabels(["20", "100", "500", "1,000"])
    axis.set_xlabel("Attempt budget")
    axis.set_ylabel("Mean effective clusters at 0.5 Å")
    axis.set_title("F   Geometric breadth", loc="left", pad=9, fontweight="bold")
    style_axis(axis)


def latex_method(recipe: str) -> str:
    return {
        "qwen": "Qwen",
        "qwen75_loqi25": "Qwen/LoQI 75:25",
        "qwen50_loqi50": "\\textbf{Qwen/LoQI 50:50}",
        "qwen25_loqi75": "Qwen/LoQI 25:75",
        "loqi": "LoQI",
        "rdkit_raw": "RDKit raw",
        "rdkit_minimized": "RDKit minimized",
    }[recipe]


def paired_contrast(
    paired: pd.DataFrame,
    family: str,
    budget: int,
    left: str,
    right: str,
) -> tuple[float, float, float]:
    rows = paired[
        (paired["family"] == family)
        & (paired["budget"] == budget)
        & (
            ((paired["left"] == left) & (paired["right"] == right))
            | ((paired["left"] == right) & (paired["right"] == left))
        )
    ]
    if len(rows) != 1:
        raise ValueError(f"Expected one contrast for {family}/{budget}/{left}/{right}")
    row = rows.iloc[0]
    if row["left"] == left:
        values = row["difference"], row["ci_low"], row["ci_high"]
    else:
        values = -row["difference"], -row["ci_high"], -row["ci_low"]
    return tuple(100 * float(value) for value in values)


def format_contrast(values: tuple[float, float, float]) -> str:
    difference, low, high = values
    text = f"{difference:+.2f} [{low:.2f}, {high:.2f}]"
    if low > 0 or high < 0:
        return f"\\textbf{{{text}}}"
    return text


def build_main_table(coverage: pd.DataFrame, paired: pd.DataFrame) -> None:
    reference = reference_coverage(coverage)
    ligand_families = ["individual_ligand", "shared_ligand"]
    means = reference.pivot_table(
        index=["family", "budget"], columns="recipe", values="coverage"
    )

    rows = []
    specifications = [
        (20, "LoQI minus Qwen", "loqi", "qwen"),
        (100, "50:50 minus stronger pure$^{a}$", "qwen50_loqi50", "best_pure"),
        (500, "50:50 minus stronger pure$^{a}$", "qwen50_loqi50", "best_pure"),
        (1000, "Qwen minus LoQI", "qwen", "loqi"),
        (1000, "50:50 minus best RDKit$^{b}$", "qwen50_loqi50", "best_rdkit"),
    ]
    for budget, label, left, right_rule in specifications:
        cells = []
        for family in ligand_families:
            if right_rule == "best_pure":
                row = means.loc[(family, budget)]
                right = max(["qwen", "loqi"], key=lambda recipe: row[recipe])
            elif right_rule == "best_rdkit":
                row = means.loc[(family, budget)]
                right = max(
                    ["rdkit_raw", "rdkit_minimized"], key=lambda recipe: row[recipe]
                )
            else:
                right = right_rule
            cells.append(format_contrast(paired_contrast(paired, family, budget, left, right)))
        rows.append(f"{budget:,} & {label} & {cells[0]} & {cells[1]} \\\\")

    content = r"""% Rebuild with results-data/cdpkit-pharmacophore-mixture/render.py.
\begin{table*}[htbp]
\centering\footnotesize
\setlength{\tabcolsep}{5.5pt}
\renewcommand{\arraystretch}{1.22}
\caption{Key ligand-only contrasts in the stringent pharmacophore comparison.}
\label{tab:pharmacophore-mixture}
\par\smallskip
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}rlll@{}}
\toprule
\textbf{Budget} & \textbf{Contrast} & \textbf{Individual ligand (pp [95\% interval])} & \textbf{Shared ligand (pp [95\% interval])} \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular*}
\par\smallskip
\begin{minipage}{\linewidth}
Queries require six of six features within 0.5~\AA{}. Intervals resample molecules (individual ligand, $N=21$; shared ligand, $N=19$). Bold values have intervals excluding zero. $^{a}$The pure generator with the higher family-level mean at that budget, not a molecule-specific oracle. $^{b}$The RDKit control with the higher family-level mean at budget 1,000.
\end{minipage}
\end{table*}
"""
    (ROOT / "tables" / "pharmacophore-mixture.tex").write_text(content)


def build_supporting_tables(coverage: pd.DataFrame, scope: pd.DataFrame) -> None:
    reference = reference_coverage(coverage)
    primary_n = (
        reference[reference["budget"] == 20]
        .groupby("family", as_index=True)["molecules"]
        .first()
    )
    scope_rows = []
    for family in FAMILIES:
        row = scope.set_index("family").loc[family]
        scope_rows.append(
            f"{FAMILY_LABELS[family]} & {int(row['queries']):,} & {int(row['parents']):,} & "
            f"{int(row['molecules_with_queries'])} & {int(primary_n.loc[family])} \\\\"
        )

    k1000 = reference[reference["budget"] == 1000].pivot(index="recipe", columns="family", values="coverage")
    all_recipes = [
        "qwen", "qwen75_loqi25", "qwen50_loqi50", "qwen25_loqi75", "loqi",
        "rdkit_raw", "rdkit_minimized",
    ]
    endpoint_rows = []
    for index, recipe in enumerate(all_recipes):
        shade = "\\rowcolor{tableShade}\n" if index % 2 == 0 else ""
        values = " & ".join(f"{100 * k1000.loc[recipe, family]:.2f}" for family in FAMILIES)
        endpoint_rows.append(f"{shade}{latex_method(recipe)} & {values} \\\\")

    content = r"""% Rebuild with results-data/cdpkit-pharmacophore-mixture/render.py.
\begin{table}[htbp]
\centering\footnotesize
\caption{Scope of the four CDPKit pharmacophore-query families.}
\label{tab:cdpkit-query-scope}
\begin{tabular}{lrrrr}
\toprule
\textbf{Family} & \textbf{Queries} & \textbf{Parents} & \textbf{Any-query molecules} & \textbf{Reference-cell molecules} \\
\midrule
""" + "\n".join(scope_rows) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{0.95\linewidth}\footnotesize
Query counts span all feature requirements and radii. The manuscript display uses the primary construction rule, six of six features, and a 0.5~\AA{} tolerance. Shared-family eligibility is lower because a shared hypothesis and held-out comparison could not be constructed for every molecule.
\end{minipage}
\end{table}

\begin{table}[htbp]
\centering\footnotesize
\caption{Expected stringent pharmacophore-query coverage at budget 1,000.}
\label{tab:cdpkit-endpoint-coverage}
\begin{tabular}{lrrrr}
\toprule
\textbf{Allocation} & \textbf{Ind. ligand} & \textbf{Ind. interaction} & \textbf{Shared ligand} & \textbf{Shared protein} \\
\midrule
""" + "\n".join(endpoint_rows) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{0.95\linewidth}\footnotesize
Values are percentages for the six-of-six, 0.5~\AA{} query. Qwen is listed first in mixed allocations. Coverage averages subsets within a parent, parents within a molecule, and then molecules equally. Unequal allocations are retained here as exploratory results and were not used to fit a production ratio.
\end{minipage}
\end{table}
"""
    (ROOT / "tables" / "supporting-pharmacophore-mixture.tex").write_text(content)


def main() -> None:
    coverage = pd.read_csv(HERE / "coverage-primary.csv")
    paired = pd.read_csv(HERE / "paired-k20.csv")
    paired_ligand = pd.read_csv(HERE / "paired-ligand-only.csv")
    diversity = pd.read_csv(HERE / "diversity-summary.csv")
    scope = pd.read_csv(HERE / "query-scope.csv")
    reference = reference_coverage(coverage)

    if len(reference) != 4 * 7 * 4:
        raise ValueError("Reference coverage archive is incomplete")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 8.7,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.edgecolor": "#8B929A",
            "axes.labelcolor": "#29323C",
            "text.color": "#29323C",
            "xtick.color": "#535D68",
            "ytick.color": "#535D68",
            "axes.linewidth": 0.65,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "cdpkit-pharmacophore-mixture",
            "savefig.facecolor": "white",
        }
    )
    figure = plt.figure(figsize=(12.4, 6.1))
    grid = GridSpec(
        2,
        3,
        figure=figure,
        width_ratios=[1.0, 1.0, 0.66],
        hspace=0.48,
        wspace=0.30,
        left=0.065,
        right=0.985,
        bottom=0.12,
        top=0.84,
    )
    axes = [
        figure.add_subplot(grid[:, 0]),
        figure.add_subplot(grid[:, 1]),
        figure.add_subplot(grid[0, 2]),
        figure.add_subplot(grid[1, 2]),
    ]
    plot_coverage(axes[0], reference, "individual_ligand", "A")
    plot_coverage(axes[1], reference, "shared_ligand", "B")
    plot_coverage(
        axes[2], reference, "individual_interaction", "C", secondary=True
    )
    plot_coverage(axes[3], reference, "shared_protein", "D", secondary=True)
    axes[1].tick_params(labelleft=False)
    axes[1].set_ylabel("")
    axes[2].tick_params(labelbottom=False)
    axes[2].set_xlabel("")
    axes[2].set_ylabel("Coverage (%)")
    axes[3].set_ylabel("Coverage (%)")

    axes[0].annotate(
        "LoQI leads Qwen at $K=20$",
        xy=(20, 49.65),
        xytext=(28, 69),
        color="#9A6500",
        fontsize=8.2,
        arrowprops=dict(arrowstyle="-", color="#C08300", linewidth=0.9),
    )
    axes[1].annotate(
        "LoQI leads Qwen at $K=20$",
        xy=(20, 64.52),
        xytext=(28, 83),
        color="#9A6500",
        fontsize=8.2,
        arrowprops=dict(arrowstyle="-", color="#C08300", linewidth=0.9),
    )
    axes[0].annotate(
        "50:50 remains at or above\nthe stronger pure mean",
        xy=(500, 82.21),
        xytext=(180, 94),
        color="#67267D",
        fontsize=8.2,
        arrowprops=dict(arrowstyle="-", color="#7B3294", linewidth=0.9),
    )

    handles = [
        Line2D(
            [],
            [],
            color=COLORS[recipe],
            marker=MARKERS[recipe],
            linestyle=LINESTYLES[recipe],
            markersize=4.5,
            linewidth=2.7 if recipe == "qwen50_loqi50" else 1.7,
            label=LABELS[recipe],
        )
        for recipe in MAIN_RECIPES
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=8.2,
        bbox_to_anchor=(0.5, 0.985),
    )

    FIGURE_BASENAME.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    metadata = {
        "png": {"Software": "matplotlib CDPKit mixture renderer"},
        "pdf": {"Creator": "matplotlib CDPKit mixture renderer", "CreationDate": None, "ModDate": None},
        "svg": {"Date": None},
    }
    for extension in ("png", "pdf", "svg"):
        destination = FIGURE_BASENAME.with_suffix(f".{extension}")
        figure.savefig(destination, dpi=300, metadata=metadata[extension])
        outputs.append(destination)
    plt.close(figure)
    svg = FIGURE_BASENAME.with_suffix(".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")

    SUPPORTING.mkdir(parents=True, exist_ok=True)
    for source_name, destination_name in [
        ("query-family-examples.pdf", "cdpkit-query-family-examples.pdf"),
        ("generator-separation-grid.pdf", "cdpkit-generator-separation-grid.pdf"),
    ]:
        destination = SUPPORTING / destination_name
        shutil.copy2(HERE / source_name, destination)
        outputs.append(destination)

    build_main_table(coverage, paired_ligand)
    build_supporting_tables(coverage, scope)
    outputs.extend(
        [
            ROOT / "tables" / "pharmacophore-mixture.tex",
            ROOT / "tables" / "supporting-pharmacophore-mixture.tex",
        ]
    )

    inputs = [
        HERE / "coverage-primary.csv",
        HERE / "paired-k20.csv",
        HERE / "paired-ligand-only.csv",
        HERE / "diversity-summary.csv",
        HERE / "query-scope.csv",
        HERE / "query-family-examples.pdf",
        HERE / "generator-separation-grid.pdf",
        Path(__file__),
    ]
    (HERE / "rendered-files.json").write_text(
        json.dumps(
            {
                "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
                "outputs": {str(path.relative_to(ROOT)): sha256(path) for path in outputs},
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
