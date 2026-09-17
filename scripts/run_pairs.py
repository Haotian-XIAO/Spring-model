#!/usr/bin/env python3
"""Canonical command-line entry point for paired 0606 experiments."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import baseline_config, run_batch, smoke_config


def selected_config(arguments: argparse.Namespace):
    """Select a named configuration and apply only explicit FF overrides."""

    config = baseline_config() if arguments.config == "baseline" else smoke_config()
    softening_values = {
        "magnitude": arguments.softening_magnitude,
        "duration": arguments.softening_duration,
        "threshold": arguments.softening_threshold,
        "recovery_rate": arguments.recovery_rate,
        "minimum_E": arguments.minimum_e,
    }
    overrides = {key: value for key, value in softening_values.items() if value is not None}
    if overrides:
        config = replace(config, softening=replace(config.softening, **overrides))
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seeds", nargs="+", type=int, help="Independent master seeds, one pair each")
    parser.add_argument("--output-root", type=Path, default=Path("experiment_results"))
    parser.add_argument("--workers", type=int, default=1, help="Process workers; default preserves simple sequential execution")
    parser.add_argument(
        "--config", choices=("baseline", "smoke"), default="baseline",
        help="baseline is the 25 x 40 / 600-phase recovered 0606 setup; smoke is a tiny software test",
    )
    parser.add_argument("--softening-magnitude", type=float, help="Fractional immediate E reduction at trigger")
    parser.add_argument("--softening-duration", type=int, help="Number of outer phases held at target E")
    parser.add_argument("--softening-threshold", type=float, help="Trigger when current k exceeds this value")
    parser.add_argument("--recovery-rate", type=float, help="Per-phase geometric E recovery fraction")
    parser.add_argument("--minimum-e", type=float, help="Lower bound on softening target E")
    parser.add_argument("--smoke-config", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error("--workers must be at least 1")
    if arguments.softening_duration is not None and arguments.softening_duration < 1:
        parser.error("--softening-duration must be at least 1")
    if arguments.softening_magnitude is not None and not 0.0 <= arguments.softening_magnitude <= 1.0:
        parser.error("--softening-magnitude must be between 0 and 1")
    if arguments.recovery_rate is not None and not 0.0 < arguments.recovery_rate <= 1.0:
        parser.error("--recovery-rate must be in (0, 1]")
    if arguments.softening_threshold is not None and arguments.softening_threshold < 0.0:
        parser.error("--softening-threshold must be non-negative")
    if arguments.minimum_e is not None and arguments.minimum_e < 0.0:
        parser.error("--minimum-e must be non-negative")
    if arguments.smoke_config:
        arguments.config = "smoke"
    config = selected_config(arguments)
    results = run_batch(arguments.seeds, config, arguments.output_root, workers=arguments.workers)
    for result in results:
        print(f"seed={result['seed']} status={result['status']}")


if __name__ == "__main__":
    main()
