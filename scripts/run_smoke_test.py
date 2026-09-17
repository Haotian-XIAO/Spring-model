#!/usr/bin/env python3
"""Small end-to-end verification; not a biological simulation ensemble."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import generate_initial_state, run_batch, smoke_config
from pulmonary_fibrosis_model.core import state_fingerprint


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root", type=Path, default=Path("experiment_results/smoke_test_20260917")
    )
    arguments = parser.parse_args()
    config = smoke_config()
    seeds = [101, 202]

    # Same seed must reconstruct the common post-preconvergence state exactly.
    first, _, _, _ = generate_initial_state(seeds[0], config)
    second, _, _, _ = generate_initial_state(seeds[0], config)
    assert state_fingerprint(first) == state_fingerprint(second), "initial realization was not reproducible"

    results = run_batch(seeds, config, arguments.output_root, workers=1)
    assert all(result["status"] == "success" for result in results), results
    for seed in seeds:
        pair_dir = arguments.output_root / f"seed_{seed}"
        control = read_rows(pair_dir / "control" / "timeseries.csv")
        softening = read_rows(pair_dir / "softening" / "timeseries.csv")
        pair = read_rows(pair_dir / "pair_timeseries.csv")
        assert (pair_dir / "common_initial_state.npz").exists()
        assert (pair_dir / "pair_summary.csv").exists()
        assert len(control) == config.agent_total_iterations + 1
        assert len(softening) == config.agent_total_iterations + 1
        assert float(control[-1]["softened_fraction_ever"]) == 0.0
        assert any(float(row["softened_fraction_ever"]) > 0.0 for row in softening)
        assert min(float(row["mean_E"]) for row in softening) < 1.0
        assert min(float(row["mean_stiffness"]) for row in softening[1:]) < float(control[1]["mean_stiffness"])
        assert any(abs(float(row["delta_mean_E"])) > 0.0 for row in pair[1:])

    print("SMOKE TEST PASSED")
    print(f"output_root={arguments.output_root}")
    print("seeds=101,202; common initialization reproducible; ON/OFF outputs and pair summaries written")


if __name__ == "__main__":
    main()
