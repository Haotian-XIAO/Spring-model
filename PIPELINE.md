# Reproducible paired 0606 experiments

The archival 0606 notebooks are preserved unchanged. The canonical new entry
point is `scripts/run_pairs.py`, backed by `pulmonary_fibrosis_model`.

```bash
python scripts/run_pairs.py 123 456 --config baseline --output-root experiment_results/run_001
```

`baseline` is the reconstructed 0606 configuration: 25 x 40 cells, alpha 0.6,
random walk length 300, k_normal 1, k_fibrosis 50, l0 ratio 0.6, 600 outer
agent phases, 2,000 preconvergence iterations, and 1,500 mechanical iterations
per phase. A deliberately reduced, non-biological development configuration is
available with `--config smoke`.

The current exploratory FF parameters can be changed without source edits:

```bash
python scripts/run_pairs.py 101 102 103 --config baseline \
  --softening-magnitude 0.05 --softening-duration 150 \
  --softening-threshold 5 --recovery-rate 0.003 --minimum-e 0.1 \
  --workers 4 --output-root experiment_results/pilot_005
```
Independent seed pairs can be placed in conservative process workers with
`--workers N`; this changes neither seed derivation nor pairing.

For one pair, `seed_<master-seed>/` contains:

- `common_initial_state.npz`: the shared post-preconvergence state;
- `pair_metadata.json`: configuration, master/derived seeds, code revision,
  canonical-source hash, and a common-state fingerprint;
- `control/` and `softening/`: time series, selected-phase raw spatial NPZ data,
  and branch metadata;
- `pair_timeseries.csv` and `pair_summary.csv`: continuous paired effects.

The RNG design uses `SeedSequence` children for biological initialization/random
walk, shared preconvergence mechanics, paired Metropolis noise, and independent
null-test mechanics. Both branches deep-copy exactly the same initial state.
For the ON/OFF comparison, each outer phase uses a fresh reproducible
phase-indexed Metropolis stream and each mechanical iteration consumes one
uniform, so both conditions receive the same phase/iteration-indexed uniform
field. This avoids downstream numerical RNG noise being misattributed to the
FF intervention even if the branches converge in different iteration counts.

`softening=False` retains the 0606 route `A <- A + P*a*D`, `k <- E*A` with
unchanged `E`. With softening enabled, a spring that first exceeds the
threshold records its pre-event E, transitions immediately to a configurable
target E, holds for the configured number of outer phases, then recovers toward
its recorded pre-event E under a tolerance-based termination rule. The FF
intervention never modifies A.

The full configuration carries forward the visible 0606 value 0.05 as the
default softening magnitude, but it was inert in the archival notebook and is
not a biological calibration. The smoke test intentionally uses 0.5 only to
make the software intervention unambiguous on a tiny network.

Cell areas are direct shoelace areas of the persisted `hex_cells` polygons and
deformed node coordinates. No preserved/remodeled cell threshold is imposed:
the pipeline exports raw areas and normalized areas only.

Compact time-series QoIs are saved at every outer agent phase. Spatial NPZ
snapshots are saved at the configurable selected-phase cadence (25 phases in
the full 0606 configuration; every phase in the smoke configuration).

Run the small verification only with:

```bash
python scripts/run_smoke_test.py
```

It uses two small networks and six agent phases. It is a software test, not a
biological result.

For a fresh Mac Studio checkout:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m unittest discover -s tests -v
```

One production-sized pair can be launched with
`python scripts/run_full_pair.py 101 --output-root experiment_results/full_pair`.
Use `--dry-run` to validate that command without starting the simulation. The
script reports total, control and softening runtimes plus peak process memory.
