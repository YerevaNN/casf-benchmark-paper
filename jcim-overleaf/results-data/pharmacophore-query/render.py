"""Render the pharmacophore-query Supporting Figure from archived summaries."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIGURE_BASENAME = ROOT / "figures" / "supporting" / "pharmacophore-query-coverage"

LABELS = {
    "qwen_1p7b_v2_coordtoken_step41054": "Qwen 1.7B FSQ",
    "qwen_1p7b_v2_coordtoken_step41054_atc48": "Qwen 1.7B FSQ",
    "loqi_druglike": "LoQI",
    "loqi_atc48": "LoQI",
    "flowr_druglike": "FlowR",
    "mcf_drugs_l_druglike": "MCF",
    "nextmol_dmt_l_druglike": "NExT-Mol",
    "torsional_diffusion_druglike": "Torsional Diffusion",
    "rdkit_random_raw_druglike": "RDKit raw",
    "rdkit_random_minimized_22_druglike": "RDKit minimized",
    "rdkit_etkdg_raw_atc48": "RDKit raw",
    "rdkit_etkdg_mmff94s_atc48": "RDKit minimized",
    "qwen_loqi_equal_mixture": "Qwen/LoQI 50:50",
}
COLORS = {
    "Qwen 1.7B FSQ": "#0072B2",
    "LoQI": "#E69F00",
    "FlowR": "#009E73",
    "MCF": "#CC79A7",
    "NExT-Mol": "#56B4E9",
    "Torsional Diffusion": "#D55E00",
    "RDKit raw": "#F0D629",
    "RDKit minimized": "#666666",
    "Qwen/LoQI 50:50": "#8E44AD",
}
CURVE_LABELS = [
    "Qwen 1.7B FSQ",
    "LoQI",
    "RDKit raw",
    "RDKit minimized",
    "Qwen/LoQI 50:50",
]


def primary(frame: pd.DataFrame) -> pd.DataFrame:
    selected = frame[
        (frame["analysis"] == "primary") & np.isclose(frame["radius"], 0.5)
    ].copy()
    selected["label"] = selected["method"].map(LABELS)
    if selected["label"].isna().any():
        raise ValueError("An archived method has no publication label")
    return selected


def style_axis(axis: plt.Axes) -> None:
    axis.grid(color="#E7EBEE", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def plot_curves(axis: plt.Axes, frame: pd.DataFrame, title: str) -> None:
    frame = primary(frame)
    for label in CURVE_LABELS:
        curve = frame[frame["label"] == label].sort_values("budget")
        if len(curve) != 7:
            raise ValueError(f"Expected seven archived budgets for {label} in {title}")
        axis.plot(
            curve["budget"],
            curve["coverage"],
            marker="o",
            markersize=3.5,
            linewidth=1.8,
            color=COLORS[label],
            label=label,
        )
    axis.set_xscale("log")
    axis.set_xlim(9, 1100)
    axis.set_xticks([10, 20, 50, 100, 200, 500, 1000])
    axis.set_xticklabels(["10", "20", "50", "100", "200", "500", "1,000"])
    axis.set_ylim(0, 1.02)
    axis.set_xlabel("Requested conformer budget")
    axis.set_ylabel("Expected query coverage")
    axis.set_title(title, loc="left", pad=10, fontweight="bold")
    style_axis(axis)


def plot_pilot_endpoint(axis: plt.Axes, frame: pd.DataFrame) -> None:
    endpoint = primary(frame)
    endpoint = endpoint[endpoint["budget"] == 1000].copy()
    order = [
        "Qwen/LoQI 50:50",
        "LoQI",
        "Qwen 1.7B FSQ",
        "Torsional Diffusion",
        "MCF",
        "FlowR",
        "NExT-Mol",
        "RDKit minimized",
        "RDKit raw",
    ]
    endpoint = endpoint.set_index("label").loc[order].reset_index()
    positions = np.arange(len(endpoint))
    axis.errorbar(
        endpoint["coverage"],
        positions,
        xerr=np.vstack(
            [
                endpoint["coverage"] - endpoint["ci_low"],
                endpoint["ci_high"] - endpoint["coverage"],
            ]
        ),
        fmt="none",
        ecolor="#899296",
        elinewidth=1.2,
        capsize=2,
    )
    for position, row in zip(positions, endpoint.itertuples(index=False), strict=True):
        axis.scatter(row.coverage, position, s=34, color=COLORS[row.label], zorder=3)
    axis.set_yticks(positions, order)
    axis.invert_yaxis()
    axis.set_xlim(0.3, 1.0)
    axis.set_xlabel("Coverage at budget 1,000 (95% molecule bootstrap CI)")
    axis.set_title("B   PLINDER-22 endpoint", loc="left", pad=10, fontweight="bold")
    style_axis(axis)


def plot_overlap(axis: plt.Axes, frame: pd.DataFrame) -> None:
    methods = ["Qwen 1.7B FSQ", "LoQI", "Qwen/LoQI 50:50"]
    groups = ["No exact identity", "Exact identity"]
    offsets = {"No exact identity": -0.13, "Exact identity": 0.13}
    markers = {"No exact identity": "o", "Exact identity": "s"}
    positions = np.arange(len(methods))
    for group in groups:
        subset = frame[frame["training_group"] == group].set_index("method_label").loc[methods]
        y = positions + offsets[group]
        axis.errorbar(
            subset["coverage"],
            y,
            xerr=np.vstack(
                [
                    subset["coverage"] - subset["ci_low"],
                    subset["ci_high"] - subset["coverage"],
                ]
            ),
            fmt=markers[group],
            markersize=5,
            color="#3D4650" if group == "No exact identity" else "#9B4D96",
            ecolor="#899296",
            elinewidth=1.1,
            capsize=2,
            label=f"{group} (n={int(subset['molecule_count'].iloc[0])})",
        )
    axis.set_yticks(positions, methods)
    axis.invert_yaxis()
    axis.set_xlim(0, 1.0)
    axis.set_xlabel("Coverage at budget 1,000")
    axis.set_title("D   Extension by exact Qwen training identity", loc="left", pad=10, fontweight="bold")
    axis.legend(frameon=False, fontsize=8, loc="upper left")
    style_axis(axis)


def main() -> None:
    pilot = pd.read_csv(HERE / "pilot-aggregate-budget-curves.csv")
    extension = pd.read_csv(HERE / "atc48-aggregate-budget-curves.csv")
    overlap = pd.read_csv(HERE / "atc48-training-overlap-endpoints.csv")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
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
            "svg.hashsalt": "pharmacophore-query-coverage",
            "savefig.facecolor": "white",
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(11.6, 8.0))
    figure.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.14, hspace=0.45, wspace=0.38)
    plot_curves(axes[0, 0], pilot, "A   PLINDER-22 budget curves")
    plot_pilot_endpoint(axes[0, 1], pilot)
    plot_curves(axes[1, 0], extension, "C   PLINDER-ATC48 extension")
    plot_overlap(axes[1, 1], overlap)
    figure.suptitle(
        "Protein-context pharmacophore coverage from ligand-only conformer libraries",
        x=0.10,
        y=0.975,
        ha="left",
        fontsize=13,
    )
    handles = [
        Line2D([], [], color=COLORS[label], marker="o", markersize=3.5, linewidth=1.8, label=label)
        for label in CURVE_LABELS
    ]
    figure.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8)
    FIGURE_BASENAME.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    output_metadata = {
        "png": {"Software": "matplotlib pharmacophore-query renderer"},
        "pdf": {
            "Creator": "matplotlib pharmacophore-query renderer",
            "CreationDate": None,
            "ModDate": None,
        },
        "svg": {"Date": None},
    }
    for extension_name in ("png", "pdf", "svg"):
        destination = FIGURE_BASENAME.with_suffix(f".{extension_name}")
        figure.savefig(
            destination,
            dpi=300,
            metadata=output_metadata[extension_name],
        )
        outputs.append(destination)
    plt.close(figure)
    svg = FIGURE_BASENAME.with_suffix(".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")

    inputs = [
        HERE / "pilot-aggregate-budget-curves.csv",
        HERE / "atc48-aggregate-budget-curves.csv",
        HERE / "atc48-training-overlap-endpoints.csv",
        Path(__file__),
    ]
    (HERE / "rendered-files.json").write_text(
        json.dumps(
            {
                "inputs": {
                    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in inputs
                },
                "outputs": {
                    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in outputs
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
