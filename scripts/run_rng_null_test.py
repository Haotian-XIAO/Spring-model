#!/usr/bin/env python3
"""Measure independent-mechanics OFF/OFF variability against FF ON/OFF effects.

This is a small software/numerical check, not a biological experiment.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import run_null_pair, run_pair, smoke_config


METRICS = ("mean_stiffness", "mean_E", "mean_strain", "mean_stress", "mean_activation", "mean_D", "mean_cell_area")


def read_summary(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seeds", nargs="*", type=int, default=[701, 702])
    parser.add_argument("--output-root", type=Path, default=Path("experiment_results/rng_null_test"))
    arguments = parser.parse_args()
    config = smoke_config()

    rows: list[dict[str, float | int]] = []
    for seed in arguments.seeds:
        null_root = arguments.output_root / "off_vs_off"
        effect_root = arguments.output_root / "off_vs_on"
        null = run_null_pair(seed, config, null_root)
        effect = run_pair(seed, config, effect_root)
        row: dict[str, float | int] = {"seed": seed}
        for metric in METRICS:
            null_value = abs(float(null[f"delta_final_{metric}"]))
            effect_value = abs(float(effect[f"delta_final_{metric}"]))
            row[f"null_abs_delta_{metric}"] = null_value
            row[f"effect_abs_delta_{metric}"] = effect_value
            row[f"effect_to_null_ratio_{metric}"] = effect_value / null_value if null_value else float("inf")
        rows.append(row)
    arguments.output_root.mkdir(parents=True, exist_ok=True)
    output = arguments.output_root / "rng_null_vs_softening.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output}")
    for row in rows:
        print(
            f"seed={row['seed']} "
            f"null |Δmean_E|={row['null_abs_delta_mean_E']:.6g}; "
            f"ON/OFF |Δmean_E|={row['effect_abs_delta_mean_E']:.6g}"
        )


if __name__ == "__main__":
    main()
