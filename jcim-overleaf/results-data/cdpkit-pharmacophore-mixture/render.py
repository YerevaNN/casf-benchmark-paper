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
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIGURE_BASENAME = ROOT / "figures" / "figure-pharmacophore-mixture"
SUPPORTING = ROOT / "figures" / "supporting"

REFERENCE_CONDITION = "4/4"
REFERENCE_RADIUS = 0.5
RECIPES = ["qwen", "qwen75_loqi25", "qwen50_loqi50", "qwen25_loqi75", "loqi"]
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


def plot_coverage(axis: plt.Axes, frame: pd.DataFrame, family: str, letter: str) -> None:
    subset = frame[frame["family"] == family]
    for recipe in RECIPES:
        curve = subset[subset["recipe"] == recipe].sort_values("budget")
        if len(curve) != 4:
            raise ValueError(f"Expected four budgets for {family}/{recipe}")
        highlighted = recipe == "qwen50_loqi50"
        axis.plot(
            curve["budget"],
            100 * curve["coverage"],
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
    axis.set_ylim(69, 100.8)
    axis.set_yticks([70, 80, 90, 100])
    axis.set_xlabel("Attempt budget")
    axis.set_ylabel("Expected query coverage (%)")
    axis.set_title(f"{letter}   {FAMILY_LABELS[family]}", loc="left", pad=9, fontweight="bold")
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


def build_main_table(coverage: pd.DataFrame, diversity: pd.DataFrame) -> None:
    reference = reference_coverage(coverage)
    k20 = reference[reference["budget"] == 20].pivot(index="recipe", columns="family", values="coverage")
    k1000 = reference[reference["budget"] == 1000].pivot(index="recipe", columns="family", values="coverage")
    deficits = (k1000.max(axis=0) - k1000).max(axis=1)
    endpoint = diversity[
        (diversity["budget"] == 1000) & np.isclose(diversity["cluster_radius"], REFERENCE_RADIUS)
    ].set_index("recipe")
    rows = []
    for recipe in RECIPES:
        cells = []
        for family in FAMILIES:
            value = 100 * k20.loc[recipe, family]
            formatted = f"{value:.2f}"
            if recipe == "qwen50_loqi50":
                formatted = f"\\textbf{{{formatted}}}"
            cells.append(formatted)
        shade = "\\rowcolor{tableShade}\n" if recipe == "qwen50_loqi50" else ""
        rows.append(
            shade
            + f"{latex_method(recipe)} & "
            + " & ".join(cells)
            + f" & {100 * deficits.loc[recipe]:.2f} & {endpoint.loc[recipe, 'effective_clusters']:.1f}"
            + f" & {100 * endpoint.loc[recipe, 'valid_yield']:.1f} \\\\"
        )
    content = r"""% Rebuild with results-data/cdpkit-pharmacophore-mixture/render.py.
\begin{table*}[htbp]
\centering\footnotesize
\setlength{\tabcolsep}{3.2pt}
\renewcommand{\arraystretch}{1.25}
\caption{Fixed Qwen/LoQI allocations in the four-family pharmacophore experiment.}
\label{tab:pharmacophore-mixture}
\par\smallskip
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrrrr@{}}
\toprule
 & \multicolumn{4}{c}{Expected coverage at budget 20 (\%)} & \multicolumn{3}{c}{Budget 1,000} \\
\cmidrule(lr){2-5}\cmidrule(l){6-8}
\textbf{Allocation} & Ind. ligand & Ind. interaction & Shared ligand & Shared protein & Max. deficit$^{a}$ & Eff. clusters$^{b}$ & Valid yield (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular*}
\par\smallskip
\begin{minipage}{\linewidth}
The reference query requires four of four features within 0.5~\AA{}. Qwen is listed first in each allocation. The 50:50 allocation is the numerical leader in all four budget-20 columns; its advantages over 75:25 are 0.02--0.19 percentage point, with paired intervals spanning zero. $^{a}$Largest shortfall from the best allocation in any of the four query families at budget 1,000, in percentage points. $^{b}$Mean exponential-Shannon effective cluster count at 0.5~\AA{}; all diversity summaries use the same 22 molecules.
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
Queries counts span all feature requirements and radii. The reference cell uses the primary construction rule, four of four features, and a 0.5~\AA{} tolerance. Shared-family eligibility is lower because a shared hypothesis and held-out comparison could not be constructed for every molecule.
\end{minipage}
\end{table}

\begin{table}[htbp]
\centering\footnotesize
\caption{Expected pharmacophore-query coverage at budget 1,000.}
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
Values are percentages for the four-of-four, 0.5~\AA{} reference query. Qwen is listed first in mixed allocations. Coverage averages subsets within a parent, parents within a molecule, and then molecules equally.
\end{minipage}
\end{table}
"""
    (ROOT / "tables" / "supporting-pharmacophore-mixture.tex").write_text(content)


def main() -> None:
    coverage = pd.read_csv(HERE / "coverage-primary.csv")
    paired = pd.read_csv(HERE / "paired-k20.csv")
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
    figure, axes = plt.subplots(2, 3, figsize=(12.4, 7.8))
    figure.subplots_adjust(left=0.075, right=0.985, top=0.89, bottom=0.145, hspace=0.47, wspace=0.36)
    for axis, family, letter in zip(axes.flat[:4], FAMILIES, "ABCD", strict=True):
        plot_coverage(axis, reference, family, letter)
    plot_differences(axes[1, 1], paired)
    plot_diversity(axes[1, 2], diversity)
    figure.suptitle(
        "Choosing a fixed Qwen/LoQI allocation for pharmacophore coverage",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=13,
    )
    handles = [
        Line2D(
            [], [], color=COLORS[recipe], marker="o", markersize=4,
            linewidth=2.7 if recipe == "qwen50_loqi50" else 1.5,
            label=LABELS[recipe],
        )
        for recipe in RECIPES
    ]
    figure.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8.2)

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

    build_main_table(coverage, diversity)
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
