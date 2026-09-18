#!/usr/bin/env python3
"""Create collaborator-facing figures from the completed N=100 analysis.

This script is deliberately read-only with respect to simulation outputs.  It
loads the existing aggregate tables and saved spatial snapshots, then writes a
new figure set without replacing the original analysis.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import numpy as np

from analysis.load_results import SeedRun, discover_seeds, load_snapshot
from analysis.morphometry import metrics as morphology_metrics


NO_SOFTENING = "#4C78A8"
WITH_SOFTENING = "#D55E00"
DIFFERENCE = "#6A3D9A"
REGION_NAMES = ("Contracted", "Preserved", "Enlarged")
REGION_COLORS = ("#0072B2", "#009E73", "#D55E00")
LOWER = 0.75
UPPER = 1.25


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 400,
        }
    )


def save(fig: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def quantiles_by_phase(
    rows: list[dict[str, str]], key: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phases = np.array(sorted({int(float(row["phase"])) for row in rows}))
    values = np.array(
        [[float(row[key]) for row in rows if int(float(row["phase"])) == phase] for phase in phases]
    )
    return phases, np.median(values, axis=1), np.percentile(values, 25, axis=1), np.percentile(values, 75, axis=1)


def main_01(rows: list[dict[str, str]], output: Path) -> None:
    panels = (
        ("delta_operational_preserved_cell_fraction", "Preserved alveoli"),
        ("delta_normalized_area_cv", "Variation in alveolar area"),
        ("endpoint_delta_mean_stress", "Average spring stress"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 2.75))
    for label, (ax, (key, title)) in zip("ABC", zip(axes, panels)):
        values = np.sort(np.array([float(row[key]) for row in rows if row.get("run_status") == "success"]))
        ax.scatter(np.arange(1, len(values) + 1), values, color=DIFFERENCE, s=15, alpha=0.82, edgecolors="none")
        ax.axhline(0, color="black", lw=1.0, zorder=0)
        ax.set_title(title)
        ax.set_xlabel("Simulations ordered by effect")
        ax.set_ylabel("Difference: With softening − No softening")
        ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    fig.suptitle("Final effects of softening", fontsize=12, y=1.02)
    fig.tight_layout()
    save(fig, output / "main_01_final_effects")


def spatial_timecourse_rows(runs: list[SeedRun]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for run in runs:
        if run.status != "success":
            continue
        for condition in ("control", "softening"):
            for path in sorted((run.path / condition).glob("spatial_phase_*.npz")):
                phase = int(path.stem.rsplit("_", 1)[1])
                values = morphology_metrics(load_snapshot(path), LOWER, UPPER)
                rows.append(
                    {
                        "seed": run.seed,
                        "phase": phase,
                        "condition": condition,
                        "preserved": values["operational_preserved_cell_fraction"],
                        "relative_area": values["normalized_area_mean"],
                        "area_variation": values["normalized_area_cv"],
                        "neighbor_difference": values["neighbor_normalized_area_contrast"],
                    }
                )
    return rows


def main_02(spatial_rows: list[dict[str, float | int]], output: Path) -> None:
    panels = (
        ("preserved", "Preserved alveoli", "Fraction of alveoli"),
        ("relative_area", "Relative alveolar area", "Relative area"),
        ("area_variation", "Variation in alveolar area", "Coefficient of variation"),
        ("neighbor_difference", "Difference between neighboring\nalveolar areas", "Mean absolute difference"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 5.6), sharex=True)
    for panel_label, (ax, (key, title, ylabel)) in zip("ABCD", zip(axes.flat, panels)):
        for condition, color, legend_label in (
            ("control", NO_SOFTENING, "No softening"),
            ("softening", WITH_SOFTENING, "With softening"),
        ):
            subset = [row for row in spatial_rows if row["condition"] == condition]
            phases = np.array(sorted({int(row["phase"]) for row in subset}))
            samples = np.array(
                [[float(row[key]) for row in subset if int(row["phase"]) == phase] for phase in phases]
            )
            median = np.median(samples, axis=1)
            lower = np.percentile(samples, 25, axis=1)
            upper = np.percentile(samples, 75, axis=1)
            ax.plot(phases, median, color=color, lw=1.7, label=legend_label)
            ax.fill_between(phases, lower, upper, color=color, alpha=0.18, linewidth=0)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Simulation phase")
        ax.text(-0.12, 1.06, panel_label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Evolution across simulations (median and middle 50%)", fontsize=12, y=1.01)
    fig.tight_layout()
    save(fig, output / "main_02_evolution")


def detailed_timecourse(time_rows: list[dict[str, str]], output: Path) -> None:
    panels = (
        ("mean_normalized_cell_area", "Relative alveolar area"),
        ("mean_stress", "Average spring stress"),
        ("mean_strain", "Average spring strain"),
        ("mean_A", "Average remodeling level"),
        ("mean_activation", "Average biological activation"),
        ("mean_D", "Average agent density"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(8.5, 5.1), sharex=True)
    for ax, (metric, title) in zip(axes.flat, panels):
        for prefix, color, label in (
            ("control", NO_SOFTENING, "No softening"),
            ("softening", WITH_SOFTENING, "With softening"),
        ):
            phases, median, lower, upper = quantiles_by_phase(time_rows, f"{prefix}_{metric}")
            ax.plot(phases, median, color=color, lw=1.5, label=label)
            ax.fill_between(phases, lower, upper, color=color, alpha=0.18, linewidth=0)
        ax.set_title(title)
        ax.set_xlabel("Simulation phase")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Detailed evolution across simulations (median and middle 50%)", fontsize=11)
    fig.tight_layout()
    save(fig, output / "supplementary_detailed_evolution")


def main_03(rows: list[dict[str, str]], output: Path) -> None:
    variables = (
        ("endpoint_delta_mean_E", "Change in\nYoung's modulus"),
        ("endpoint_delta_mean_stiffness", "Change in\nspring stiffness"),
        ("endpoint_delta_mean_stress", "Change in average\nspring stress"),
        ("endpoint_delta_mean_activation", "Change in biological\nactivation"),
        ("endpoint_delta_mean_D", "Change in\nagent density"),
        ("endpoint_delta_mean_normalized_cell_area", "Change in\nalveolar area"),
    )
    successful = [row for row in rows if row.get("run_status") == "success"]
    fig, axes = plt.subplots(1, 5, figsize=(11.5, 2.55))
    for panel_label, (ax, ((x_key, x_label), (y_key, y_label))) in zip("ABCDE", zip(axes, zip(variables, variables[1:]))):
        x = np.array([float(row[x_key]) for row in successful])
        y = np.array([float(row[y_key]) for row in successful])
        ax.scatter(x, y, s=15, color=DIFFERENCE, alpha=0.72, edgecolors="none")
        ax.axhline(0, color="0.55", lw=0.7, zorder=0)
        ax.axvline(0, color="0.55", lw=0.7, zorder=0)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.text(-0.18, 1.08, panel_label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    fig.suptitle("Associations along the proposed softening pathway", fontsize=12, y=1.04)
    fig.tight_layout()
    save(fig, output / "main_03_proposed_pathway")


def snapshot_file(run: SeedRun, condition: str, phase: int) -> Path:
    path = run.path / condition / f"spatial_phase_{phase:04d}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def main_04(runs: list[SeedRun], representatives: list[dict[str, str]], output: Path) -> None:
    category_titles = {
        "least_protective": "Smallest effect",
        "median_like": "Median-like effect",
        "strongest_protective": "Strongest effect",
    }
    chosen = [row for row in representatives if row["category"] in category_titles]
    chosen.sort(key=lambda row: ("least_protective", "median_like", "strongest_protective").index(row["category"]))
    run_by_seed = {run.seed: run for run in runs}
    phases = (0, 200, 400, 600)
    phase_titles = ("Initial", "Early/mid", "Late", "Final")
    cmap = ListedColormap(REGION_COLORS)
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5], cmap.N)
    fig = plt.figure(figsize=(10.0, 8.7))
    grid = fig.add_gridspec(
        2 * len(chosen), len(phases) + 1, width_ratios=(0.78, 1, 1, 1, 1), hspace=0.15, wspace=0.025
    )
    axes = np.empty((2 * len(chosen), len(phases)), dtype=object)
    for row_index in range(2 * len(chosen)):
        for column in range(len(phases)):
            axes[row_index, column] = fig.add_subplot(grid[row_index, column + 1])
    for block, representative in enumerate(chosen):
        seed = int(representative["seed"])
        run = run_by_seed[seed]
        snapshots = {
            (condition, phase): load_snapshot(snapshot_file(run, condition, phase))
            for condition in ("control", "softening")
            for phase in phases
        }
        all_points = np.concatenate([snapshot["points"] for snapshot in snapshots.values()])
        xmin, ymin = all_points.min(axis=0)
        xmax, ymax = all_points.max(axis=0)
        padx = 0.02 * (xmax - xmin)
        pady = 0.02 * (ymax - ymin)
        for condition_index, (condition, condition_title) in enumerate(
            (("control", "No softening"), ("softening", "With softening"))
        ):
            row_index = 2 * block + condition_index
            for column, (phase, phase_title) in enumerate(zip(phases, phase_titles)):
                ax = axes[row_index, column]
                snapshot = snapshots[(condition, phase)]
                normalized_area = np.asarray(snapshot["normalized_cell_area"], float)
                classes = np.where(normalized_area < LOWER, 0, np.where(normalized_area <= UPPER, 1, 2))
                polygons = snapshot["points"][snapshot["hex_cells"]]
                collection = PolyCollection(
                    polygons, array=classes, cmap=cmap, norm=norm, edgecolor="white", linewidth=0.055
                )
                ax.add_collection(collection)
                ax.set_xlim(xmin - padx, xmax + padx)
                ax.set_ylim(ymin - pady, ymax + pady)
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                if block == 0 and condition_index == 0:
                    ax.set_title(f"{phase_title}\nPhase {phase}", fontsize=9)
                if column == 0:
                    ax.set_ylabel(condition_title, fontsize=9, labelpad=8)
        label_axis = fig.add_subplot(grid[2 * block:2 * block + 2, 0])
        label_axis.set_axis_off()
        label_axis.text(
            0.04,
            0.5,
            f"{category_titles[representative['category']]}\nSimulation {seed}",
            fontsize=9,
            fontweight="bold",
            va="center",
            ha="left",
        )
    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=8, color=color, label=name)
        for name, color in zip(REGION_NAMES, REGION_COLORS)
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, 0.012), ncol=3, frameon=False)
    fig.suptitle("Representative simulated morphology", fontsize=12, y=0.995)
    fig.subplots_adjust(left=0.015, right=0.995, top=0.93, bottom=0.06)
    save(fig, output / "main_04_representative_morphology")


def main_05(convergence_rows: list[dict[str, str]], output: Path) -> None:
    n = np.array([int(row["N"]) for row in convergence_rows])
    panels = (
        ("preserved_effect", "Difference in preserved alveoli"),
        ("area_cv_effect", "Difference in variation in alveolar area"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 2.9))
    for panel_label, (ax, (key, title)) in zip("AB", zip(axes, panels)):
        median = np.array([float(row[f"{key}_median"]) for row in convergence_rows])
        lower = np.array([float(row[f"{key}_q25"]) for row in convergence_rows])
        upper = np.array([float(row[f"{key}_q75"]) for row in convergence_rows])
        ax.plot(n, median, "o-", color=DIFFERENCE, lw=1.5, ms=4)
        ax.fill_between(n, lower, upper, color=DIFFERENCE, alpha=0.18, linewidth=0)
        ax.axhline(0, color="black", lw=0.9, zorder=0)
        ax.set_xticks(n)
        ax.set_xlabel("Number of simulations")
        ax.set_ylabel("Difference: With softening − No softening")
        ax.set_title(title)
        ax.text(-0.14, 1.08, panel_label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    fig.suptitle("Stability of the softening effect with increasing simulations", fontsize=11.5, y=1.02)
    fig.tight_layout()
    save(fig, output / "main_05_convergence")


def write_text(output: Path) -> None:
    (output / "README.md").write_text(
        """# Simplified collaborator-facing N=100 figures

These figures were generated from the completed N=100 paired simulations. No simulations were rerun, and the model equations, thresholds, statistics, and numerical results were not changed. In every comparison, **Difference = With softening − No softening**.

Alveolar regions are classified from area relative to the same region's initial area:

- Contracted: relative area < 0.75
- Preserved: 0.75 <= relative area <= 1.25
- Enlarged: relative area > 1.25

These are operational, model-based geometric definitions and have not yet been validated histologically.
""",
        encoding="utf-8",
    )
    (output / "figure_captions.md").write_text(
        """# Figure captions

**Main 01 — Final effects of softening.** Each point is one paired stochastic simulation. Values are calculated as With softening minus No softening. The zero line marks no difference. Simulations are ordered independently from the smallest to the largest paired effect in each panel.

**Main 02 — Evolution across simulations.** Lines show the median across simulations; shaded areas show the middle 50% of results. Preserved alveoli are the fraction satisfying the fixed preserved criterion. Relative alveolar area is area divided by initial area. Variation is the existing normalized-area coefficient of variation, and neighboring-area difference is the existing mean absolute contrast between adjacent alveoli.

**Main 03 — Associations along the proposed softening pathway.** Each point is one paired simulation at the final phase, and each change is With softening minus No softening. These are paired associations consistent with a proposed mechanical-biological pathway; they do not demonstrate a causal sequence.

**Main 04 — Representative simulated morphology.** Three preselected paired simulations show the smallest, median-like, and strongest preserved-alveoli effects at four saved phases. Colors show Contracted, Preserved, and Enlarged regions only; no mechanical or biological variable is overlaid.

**Main 05 — Stability of the softening effect with increasing simulations.** Lines show the median paired effect and shaded areas show the middle 50% at N = 24, 48, 72, and 100. Difference = With softening − No softening.

The Contracted, Preserved, and Enlarged classes use operational, model-based geometric definitions (relative area < 0.75, 0.75–1.25, and > 1.25, respectively) and have not yet been histologically validated.
""",
        encoding="utf-8",
    )


def validate_final_morphology(runs: list[SeedRun], seed_rows: list[dict[str, str]]) -> float:
    """Confirm recomputed phase-600 morphology matches the frozen master table."""
    row_by_seed = {int(row["seed"]): row for row in seed_rows if row.get("run_status") == "success"}
    keys = (
        "operational_preserved_cell_fraction",
        "normalized_area_mean",
        "normalized_area_cv",
        "neighbor_normalized_area_contrast",
    )
    differences = []
    for run in runs:
        if run.status != "success":
            continue
        for condition in ("control", "softening"):
            recomputed = morphology_metrics(load_snapshot(snapshot_file(run, condition, 600)), LOWER, UPPER)
            for key in keys:
                stored = float(row_by_seed[run.seed][f"{condition}_{key}"])
                differences.append(abs(recomputed[key] - stored))
    maximum = max(differences, default=float("nan"))
    if not np.isfinite(maximum) or maximum > 1e-12:
        raise RuntimeError(f"Final morphology validation failed; maximum absolute difference = {maximum}")
    return maximum


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, default=Path("experiment_results/pilot_m005_n100"))
    parser.add_argument("--analysis-root", type=Path, default=Path("analysis_results/pilot_m005_n100/ensemble"))
    parser.add_argument("--output-root", type=Path, default=Path("analysis_results/pilot_m005_n100/ensemble_simple"))
    parser.add_argument(
        "--convergence-table",
        type=Path,
        default=Path("collaborator_package/FF_softening_pilot_N100/03_numerical_summary/convergence_summary.csv"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output_root
    seed_rows = read_csv(args.analysis_root / "tables" / "master_seed_summary.csv")
    time_rows = read_csv(args.analysis_root / "tables" / "master_timecourse_summary.csv")
    representative_rows = read_csv(args.analysis_root / "tables" / "representative_seed_catalog.csv")
    convergence_rows = read_csv(args.convergence_table)
    runs = discover_seeds(args.experiment_root)
    if len([run for run in runs if run.status == "success"]) != 100:
        raise RuntimeError("Expected exactly 100 complete paired simulations")

    style()
    main_01(seed_rows, output)
    spatial_rows = spatial_timecourse_rows(runs)
    main_02(spatial_rows, output)
    detailed_timecourse(time_rows, output)
    main_03(seed_rows, output)
    main_04(runs, representative_rows, output)
    main_05(convergence_rows, output)
    write_text(output)
    maximum_difference = validate_final_morphology(runs, seed_rows)
    print(f"Wrote revised figures to {output}")
    print(f"Final morphology agrees with the frozen master table; maximum absolute difference = {maximum_difference:.3g}")


if __name__ == "__main__":
    main()
