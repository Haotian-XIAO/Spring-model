#!/usr/bin/env python3
"""Run one bounded full-geometry branch and print any first-instability report."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pulmonary_fibrosis_model import baseline_config, generate_initial_state
from pulmonary_fibrosis_model.core import NumericalInstabilityError
from pulmonary_fibrosis_model.experiment import simulate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=int)
    parser.add_argument("--phases", type=int, default=160)
    parser.add_argument("--softening", action="store_true")
    parser.add_argument("--legacy-remodeling", action="store_true", help="Use archival additive A update for forensic comparison")
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    config = replace(
        baseline_config(),
        agent_total_iterations=args.phases,
        remodeling_law="legacy_additive" if args.legacy_remodeling else "positive_degradation",
    )
    state, _, _, streams = generate_initial_state(args.seed, config)
    condition = "softening" if args.softening else "control"
    try:
        simulate(state, config, args.softening, streams["paired_mechanics_common_random_numbers"], args.output_root, condition)
    except NumericalInstabilityError as error:
        print(error)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
