#!/usr/bin/env python3
"""Run one full recovered-0606 (25 x 40, 600-phase) paired experiment."""

from __future__ import annotations

import argparse
import platform
from pathlib import Path
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import baseline_config, run_pair


def peak_memory_mib() -> float:
    """Return process peak resident size on macOS/Linux where resource supports it."""

    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    return value / (1024.0 * 1024.0) if platform.system() == "Darwin" else value / 1024.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=int)
    parser.add_argument("--output-root", type=Path, default=Path("experiment_results/full_pair"))
    parser.add_argument("--dry-run", action="store_true", help="Validate the full configuration without starting mechanics")
    arguments = parser.parse_args()
    config = baseline_config()
    if arguments.dry_run:
        print(
            f"validated full configuration: nx={config.nx}, ny={config.ny}, "
            f"agent_phases={config.agent_total_iterations}, "
            f"preconvergence_iterations={config.preconvergence_iterations}, "
            f"mechanical_iterations_per_phase={config.agent_update_interval}"
        )
        print(f"would write {arguments.output_root / f'seed_{arguments.seed}'}")
        return

    started = time.perf_counter()
    summary = run_pair(arguments.seed, config, arguments.output_root)
    total = time.perf_counter() - started
    print(f"seed={arguments.seed} status={summary['status']}")
    print(f"total_wall_time_seconds={total:.3f}")
    print(f"control_runtime_seconds={summary['control_runtime_seconds']:.3f}")
    print(f"softening_runtime_seconds={summary['softening_runtime_seconds']:.3f}")
    print(f"peak_process_memory_MiB={peak_memory_mib():.1f}")


if __name__ == "__main__":
    main()
