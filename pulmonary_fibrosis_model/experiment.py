"""Paired-run and batch drivers for controlled FF-softening experiments."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import time
from typing import Any, Iterable

import numpy as np

from .config import ModelConfig
from .core import (
    EquilibriumInfo,
    ModelState,
    apply_softening,
    make_initial_state,
    mechanical_equilibrate,
    state_fingerprint,
    update_activation_and_density,
)
from .io import state_summary, write_csv, write_json, write_spatial_snapshot


PAIR_METRICS = (
    "mean_stiffness",
    "mean_E",
    "mean_A",
    "mean_strain",
    "mean_stress",
    "mean_activation",
    "mean_D",
    "softened_fraction_current",
    "softened_fraction_ever",
    "total_cell_area",
    "mean_cell_area",
    "median_cell_area",
    "mean_normalized_cell_area",
    "median_normalized_cell_area",
)


@dataclass
class SimulationResult:
    condition: str
    state: ModelState
    rows: list[dict[str, Any]]
    start_fingerprint: str
    end_fingerprint: str
    runtime_seconds: float


def _git_value(arguments: list[str]) -> str | None:
    try:
        return subprocess.check_output(["git", *arguments], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _config_id(config: ModelConfig) -> str:
    payload = json.dumps(asdict(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _source_hash() -> str:
    """Hash canonical Python sources when work has not yet been committed."""

    package_root = Path(__file__).resolve().parent
    digest = hashlib.sha256()
    for path in sorted(package_root.glob("*.py")):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _stream_metadata(streams: list[np.random.SeedSequence]) -> dict[str, int]:
    labels = (
        "biology_initialization",
        "shared_preconvergence_mechanics",
        "paired_mechanics_common_random_numbers",
        "control_mechanics_independent",
        "softening_mechanics_independent",
    )
    return {
        label: int(sequence.generate_state(1, dtype=np.uint64)[0])
        for label, sequence in zip(labels, streams, strict=True)
    }


def generate_initial_state(
    seed: int, config: ModelConfig
) -> tuple[ModelState, EquilibriumInfo, dict[str, Any], dict[str, int]]:
    """Create the one common post-preconvergence state used by both pair branches."""

    master = np.random.SeedSequence(seed)
    biology_ss, shared_ss, common_ss, control_ss, softening_ss = master.spawn(5)
    state, pre_info = make_initial_state(
        config, np.random.default_rng(biology_ss), np.random.default_rng(shared_ss)
    )
    metadata = {
        "master_seed": int(seed),
        "derived_stream_seeds": _stream_metadata(
            [biology_ss, shared_ss, common_ss, control_ss, softening_ss]
        ),
        "initial_state_fingerprint": state_fingerprint(state),
        "preconvergence": asdict(pre_info),
    }
    return state, pre_info, metadata, metadata["derived_stream_seeds"]


def _phase_mechanics_rng(stream_seed: int, phase: int) -> np.random.Generator:
    """A reproducible Metropolis stream indexed independently for each phase."""

    return np.random.default_rng(np.random.SeedSequence([stream_seed, phase]))


def simulate(
    state: ModelState,
    config: ModelConfig,
    enable_softening: bool,
    mechanics_stream_seed: int,
    output_dir: Path,
    condition: str | None = None,
) -> SimulationResult:
    """Run one branch; state is owned exclusively by this condition."""

    condition = condition or ("softening" if enable_softening else "control")
    output_dir.mkdir(parents=True, exist_ok=False)
    start_hash = state_fingerprint(state)
    started = time.perf_counter()
    rows: list[dict[str, Any]] = [state_summary(state, phase=0)]
    write_spatial_snapshot(output_dir / "spatial_phase_0000.npz", state, phase=0)
    free_nodes = np.where(state.boundary_node_ids == 0)[0]

    for phase in range(1, config.agent_total_iterations + 1):
        # This retains the 0606 agent/area law.  Only E handling varies by condition.
        stage = "softening" if enable_softening else "control"
        update_activation_and_density(state, config, stage=stage, outer_phase=phase, condition=condition)
        apply_softening(state, config, enable_softening, outer_phase=phase, condition=condition)
        state.points, equilibrium = mechanical_equilibrate(
            state.points,
            state.edges,
            config.l0,
            state.spring_constants,
            free_nodes,
            config,
            _phase_mechanics_rng(mechanics_stream_seed, phase),
            config.agent_update_interval,
            stage=stage,
            outer_phase=phase,
            condition=condition,
            E=state.E,
            A=state.A,
        )
        # Phase-level QOIs are sampled after the mechanical solve.
        from .core import compute_strain_and_stress

        state.strain, state.stress = compute_strain_and_stress(
            state.points, state.edges, config.l0, state.spring_constants
        )
        if phase % config.output_interval == 0 or phase == config.agent_total_iterations:
            rows.append(state_summary(state, phase=phase, equilibrium=equilibrium))
        if phase % config.spatial_output_interval == 0 or phase == config.agent_total_iterations:
            write_spatial_snapshot(output_dir / f"spatial_phase_{phase:04d}.npz", state, phase=phase)

    write_csv(output_dir / "timeseries.csv", rows)
    write_json(
        output_dir / "run_metadata.json",
        {
            "condition": condition,
            "enable_softening": enable_softening,
            "initial_state_fingerprint": start_hash,
            "final_state_fingerprint": state_fingerprint(state),
            "output_sampling": "phase 0 is the common post-preconvergence state; later rows are post-equilibrium.",
        },
    )
    return SimulationResult(
        condition=condition,
        state=state,
        rows=rows,
        start_fingerprint=start_hash,
        end_fingerprint=state_fingerprint(state),
        runtime_seconds=time.perf_counter() - started,
    )


def _pair_tables(
    left: SimulationResult, right: SimulationResult
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build a continuous right-minus-left comparison table for two branches."""

    left_by_phase = {int(row["phase"]): row for row in left.rows}
    right_by_phase = {int(row["phase"]): row for row in right.rows}
    common_phases = sorted(left_by_phase.keys() & right_by_phase.keys())
    pair_rows: list[dict[str, Any]] = []
    for phase in common_phases:
        row: dict[str, Any] = {"phase": phase}
        for metric in PAIR_METRICS:
            left_value = float(left_by_phase[phase][metric])
            right_value = float(right_by_phase[phase][metric])
            row[f"{left.condition}_{metric}"] = left_value
            row[f"{right.condition}_{metric}"] = right_value
            row[f"delta_{metric}"] = right_value - left_value
        pair_rows.append(row)

    final = pair_rows[-1]
    summary: dict[str, Any] = {
        "final_phase": int(final["phase"]),
        f"{left.condition}_runtime_seconds": left.runtime_seconds,
        f"{right.condition}_runtime_seconds": right.runtime_seconds,
    }
    for metric in PAIR_METRICS:
        deltas = np.asarray([float(row[f"delta_{metric}"]) for row in pair_rows])
        summary[f"{left.condition}_final_{metric}"] = final[f"{left.condition}_{metric}"]
        summary[f"{right.condition}_final_{metric}"] = final[f"{right.condition}_{metric}"]
        summary[f"delta_final_{metric}"] = final[f"delta_{metric}"]
        summary[f"max_abs_delta_{metric}"] = float(np.nanmax(np.abs(deltas)))
    return pair_rows, summary


def _run_two_conditions(
    seed: int,
    config: ModelConfig,
    output_root: str | Path,
    *,
    right_enable_softening: bool,
    common_random_numbers: bool,
    left_condition: str,
    right_condition: str,
    comparison_kind: str,
) -> dict[str, Any]:
    """Run two deep-copied branches from precisely one common initial state."""

    root = Path(output_root) / f"seed_{seed}"
    if root.exists():
        raise FileExistsError(f"Refusing to overwrite existing pair output: {root}")
    root.mkdir(parents=True)

    initial, _, seed_metadata, stream_seeds = generate_initial_state(seed, config)
    if common_random_numbers:
        left_stream_seed = right_stream_seed = stream_seeds["paired_mechanics_common_random_numbers"]
    else:
        left_stream_seed = stream_seeds["control_mechanics_independent"]
        right_stream_seed = stream_seeds["softening_mechanics_independent"]
    config_id = _config_id(config)
    write_json(
        root / "pair_metadata.json",
        {
            **seed_metadata,
            "seed": int(seed),
            "configuration_id": config_id,
            "configuration": asdict(config),
            "code_git_commit": _git_value(["rev-parse", "HEAD"]),
            "code_git_dirty": bool(_git_value(["status", "--porcelain"])),
            "canonical_source_sha256": _source_hash(),
            "comparison_kind": comparison_kind,
            "pairing_design": "One biological/random-walk realization and one post-preconvergence state are deep-copied into both branches.",
            "mechanical_rng_design": (
                "Both conditions use identical phase/iteration-indexed Metropolis uniforms (common random numbers)."
                if common_random_numbers
                else "Conditions use independent phase-indexed Metropolis streams."
            ),
        },
    )
    write_spatial_snapshot(root / "common_initial_state.npz", initial, phase=0)
    left = simulate(
        initial.copy(), config, False, left_stream_seed, root / left_condition, left_condition
    )
    right = simulate(
        initial.copy(), config, right_enable_softening, right_stream_seed, root / right_condition, right_condition
    )
    if left.start_fingerprint != right.start_fingerprint:
        raise AssertionError("Paired branches did not start from identical state content.")

    pair_rows, summary = _pair_tables(left, right)
    summary.update(
        {
            "seed": int(seed),
            "configuration_id": config_id,
            "comparison_kind": comparison_kind,
            "status": "success",
            "common_initial_state_fingerprint": left.start_fingerprint,
            f"{left_condition}_final_state_fingerprint": left.end_fingerprint,
            f"{right_condition}_final_state_fingerprint": right.end_fingerprint,
            "same_initial_state": True,
        }
    )
    write_csv(root / "pair_timeseries.csv", pair_rows)
    write_csv(root / "pair_summary.csv", [summary])
    return summary


def run_pair(
    seed: int,
    config: ModelConfig,
    output_root: str | Path = "experiment_results",
) -> dict[str, Any]:
    """Generate one shared state, then run control and FF-softening branches."""
    return _run_two_conditions(
        seed,
        config,
        output_root,
        right_enable_softening=True,
        common_random_numbers=True,
        left_condition="control",
        right_condition="softening",
        comparison_kind="FF-softening ON minus OFF",
    )


def run_null_pair(
    seed: int,
    config: ModelConfig,
    output_root: str | Path = "experiment_results",
) -> dict[str, Any]:
    """Run independent-mechanics OFF/OFF branches to quantify solver RNG noise."""

    return _run_two_conditions(
        seed,
        config,
        output_root,
        right_enable_softening=False,
        common_random_numbers=False,
        left_condition="control_a",
        right_condition="control_b",
        comparison_kind="mechanics-RNG null: OFF versus OFF",
    )


def _batch_worker(seed: int, config: ModelConfig, output_root: str) -> dict[str, Any]:
    return run_pair(seed, config, output_root)


def run_batch(
    seeds: Iterable[int],
    config: ModelConfig,
    output_root: str | Path = "experiment_results",
    *,
    workers: int = 1,
) -> list[dict[str, Any]]:
    """Run independent pairs sequentially by default, or in conservative processes."""

    seeds = [int(seed) for seed in seeds]
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []

    def failed(seed: int, error: Exception) -> dict[str, Any]:
        return {
            "seed": seed,
            "configuration_id": _config_id(config),
            "status": "failed",
            "error_type": type(error).__name__,
            "error_message": str(error),
        }

    if workers <= 1:
        for seed in seeds:
            try:
                results.append(run_pair(seed, config, root))
            except Exception as error:  # Batch reporting must retain failed seeds.
                results.append(failed(seed, error))
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_batch_worker, seed, config, str(root)): seed for seed in seeds
            }
            for future in as_completed(futures):
                seed = futures[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    results.append(failed(seed, error))
    results.sort(key=lambda row: int(row["seed"]))
    write_csv(root / "batch_summary.csv", results)
    write_json(
        root / "batch_metadata.json",
        {
            "configuration_id": _config_id(config),
            "seeds": seeds,
            "workers": workers,
            "parallelism_note": "Independent seed pairs may run in separate processes; pairing is unaffected by completion order.",
        },
    )
    return results
