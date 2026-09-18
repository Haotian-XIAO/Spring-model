"""Raw quantitative output and directly supported QoI extraction."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .core import EquilibriumInfo, ModelState, polygon_areas


def _safe_mean(values: np.ndarray) -> float:
    return float(np.mean(values)) if len(values) else float("nan")


def component_metrics(edges: np.ndarray, selected: np.ndarray) -> dict[str, float]:
    """Connected components of selected springs, connected by a shared node."""

    selected_indices = np.flatnonzero(selected)
    if not len(selected_indices):
        return {
            "component_count": 0.0,
            "largest_component_springs": 0.0,
            "largest_component_fraction": 0.0,
        }
    parent = np.arange(len(selected_indices), dtype=np.int32)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    first_at_node: dict[int, int] = {}
    for local_index, edge_index in enumerate(selected_indices):
        for node in edges[edge_index]:
            node = int(node)
            if node in first_at_node:
                union(local_index, first_at_node[node])
            else:
                first_at_node[node] = local_index
    sizes: dict[int, int] = {}
    for index in range(len(selected_indices)):
        root = find(index)
        sizes[root] = sizes.get(root, 0) + 1
    largest = max(sizes.values())
    return {
        "component_count": float(len(sizes)),
        "largest_component_springs": float(largest),
        "largest_component_fraction": float(largest / len(selected_indices)),
    }


def state_summary(
    state: ModelState,
    phase: int,
    equilibrium: EquilibriumInfo | None = None,
) -> dict[str, float | int]:
    """Compact phase-level table; no preserved/remodeled cell threshold is invented."""

    areas = polygon_areas(state.points, state.hex_cells)
    reference = state.initial_cell_areas
    normalized = np.divide(areas, reference, out=np.full_like(areas, np.nan), where=reference != 0)
    row: dict[str, float | int] = {
        "phase": phase,
        "initial_fibrotic_fraction": float(np.mean(state.initial_fibrotic)),
        "mean_stiffness": float(np.mean(state.spring_constants)),
        "min_stiffness": float(np.min(state.spring_constants)),
        "max_stiffness": float(np.max(state.spring_constants)),
        "mean_E": float(np.mean(state.E)),
        "min_E": float(np.min(state.E)),
        "max_E": float(np.max(state.E)),
        "mean_A": float(np.mean(state.A)),
        "min_A": float(np.min(state.A)),
        "max_A": float(np.max(state.A)),
        "mean_strain": float(np.mean(state.strain)),
        "mean_stress": float(np.mean(state.stress)),
        "mean_activation": float(np.mean(state.activation)),
        "mean_D": float(np.mean(state.D)),
        "softened_fraction_current": float(np.mean(state.softening_active)),
        "softened_fraction_ever": float(np.mean(state.ever_softened)),
        "softened_count_current": int(np.count_nonzero(state.softening_active)),
        "softened_count_ever": int(np.count_nonzero(state.ever_softened)),
        "total_cell_area": float(np.sum(areas)),
        "mean_cell_area": float(np.mean(areas)),
        "median_cell_area": float(np.median(areas)),
        "mean_normalized_cell_area": float(np.nanmean(normalized)),
        "median_normalized_cell_area": float(np.nanmedian(normalized)),
    }
    for group_name, mask in {
        "current_softened": state.softening_active,
        "ever_softened": state.ever_softened,
        "non_softened": ~state.ever_softened,
    }.items():
        for metric_name, values in {
            "stiffness": state.spring_constants,
            "E": state.E,
            "A": state.A,
            "strain": state.strain,
            "stress": state.stress,
            "activation": state.activation,
            "D": state.D,
        }.items():
            row[f"{group_name}_mean_{metric_name}"] = _safe_mean(values[mask])
    row.update({f"initial_fibrotic_{key}": value for key, value in component_metrics(state.edges, state.initial_fibrotic).items()})
    if equilibrium is not None:
        row.update(
            {
                "mechanical_iterations": equilibrium.iterations,
                "mechanical_stagnated": int(equilibrium.converged_by_stagnation),
                "mechanical_accepted_uphill": equilibrium.accepted_uphill,
                "mechanical_rejected_uphill": equilibrium.rejected_uphill,
                "mechanical_final_energy": equilibrium.final_energy,
            }
        )
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0])
    for row in rows[1:]:
        fieldnames.extend(key for key in row if key not in fieldnames)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, contents: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(contents, stream, indent=2, sort_keys=True)
        stream.write("\n")


def write_spatial_snapshot(path: Path, state: ModelState, phase: int) -> None:
    """Compact spatial raw data, including cell areas and all FF state variables."""

    cell_areas = polygon_areas(state.points, state.hex_cells)
    normalized = np.divide(
        cell_areas,
        state.initial_cell_areas,
        out=np.full_like(cell_areas, np.nan),
        where=state.initial_cell_areas != 0,
    )
    np.savez_compressed(
        path,
        phase=np.asarray(phase),
        points=state.points,
        edges=state.edges,
        hex_cells=state.hex_cells,
        boundary_node_ids=state.boundary_node_ids,
        boundary_edge_ids=state.boundary_edge_ids,
        spring_stiffness=state.spring_constants,
        E=state.E,
        A=state.A,
        strain=state.strain,
        stress=state.stress,
        activation=state.activation,
        D=state.D,
        initial_fibrotic=state.initial_fibrotic,
        softening_active=state.softening_active,
        ever_softened=state.ever_softened,
        recovery_active=state.recovery_active,
        soften_age=state.soften_age,
        E_pre_softening=state.E_pre_softening,
        E_target=state.E_target,
        cell_area=cell_areas,
        normalized_cell_area=normalized,
    )
