#!/usr/bin/env python3
"""Continue a saved full-geometry branch for a bounded phase interval."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import baseline_config, generate_initial_state
from pulmonary_fibrosis_model.core import ModelState, apply_softening, compute_strain_and_stress, mechanical_equilibrate, update_activation_and_density
from pulmonary_fibrosis_model.io import state_summary, write_csv, write_spatial_snapshot


def load_state(path: Path) -> ModelState:
    raw = np.load(path)
    return ModelState(
        points=raw["points"], edges=raw["edges"], hex_cells=raw["hex_cells"],
        boundary_node_ids=raw["boundary_node_ids"], boundary_edge_ids=raw["boundary_edge_ids"],
        spring_constants=raw["spring_stiffness"], A=raw["A"], E=raw["E"], strain=raw["strain"], stress=raw["stress"],
        activation=raw["activation"], D=raw["D"], initial_fibrotic=raw["initial_fibrotic"],
        softening_active=raw["softening_active"], ever_softened=raw["ever_softened"], recovery_active=raw["recovery_active"],
        soften_age=raw["soften_age"], E_pre_softening=raw["E_pre_softening"], E_target=raw["E_target"],
        initial_cell_areas=raw["cell_area"] / raw["normalized_cell_area"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=int)
    parser.add_argument("start_phase", type=int)
    parser.add_argument("end_phase", type=int)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--softening", action="store_true")
    args = parser.parse_args()
    config = baseline_config()
    state = load_state(args.input)
    _, _, _, seeds = generate_initial_state(args.seed, config)
    condition = "softening" if args.softening else "control"
    free = np.where(state.boundary_node_ids == 0)[0]
    rows = []
    for phase in range(args.start_phase + 1, args.end_phase + 1):
        update_activation_and_density(state, config, stage=condition, outer_phase=phase, condition=condition)
        apply_softening(state, config, args.softening, outer_phase=phase, condition=condition)
        rng = np.random.default_rng(np.random.SeedSequence([seeds["paired_mechanics_common_random_numbers"], phase]))
        state.points, info = mechanical_equilibrate(
            state.points, state.edges, config.l0, state.spring_constants, free, config, rng,
            config.agent_update_interval, stage=condition, outer_phase=phase, condition=condition, E=state.E, A=state.A,
        )
        state.strain, state.stress = compute_strain_and_stress(state.points, state.edges, config.l0, state.spring_constants)
        rows.append(state_summary(state, phase, info))
    args.output.mkdir(parents=True, exist_ok=False)
    write_spatial_snapshot(args.output / f"spatial_phase_{args.end_phase:04d}.npz", state, args.end_phase)
    write_csv(args.output / "timeseries.csv", rows)
    print(f"completed {condition} phases {args.start_phase + 1}-{args.end_phase}; min_A={min(float(r['min_A']) for r in rows):.17g}; min_k={min(float(r['min_stiffness']) for r in rows):.17g}")


if __name__ == "__main__":
    main()
