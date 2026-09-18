# Model parameters and morphometry definitions

**Source of truth:** this document is derived only from the current code in
`pulmonary_fibrosis_model/` and `analysis/` at this commit
(`pulmonary_fibrosis_model/config.py`, `core.py`, `io.py`;
`analysis/standard_reference.py`, `morphometry.py`). It does not derive
anything from the historical notebooks in `archive/legacy/`, which are
provenance-only and may disagree with the values below.

## Network and geometry (`ModelConfig`, `pulmonary_fibrosis_model/config.py`)

| Parameter | Value | Meaning |
|---|---|---|
| `nx`, `ny` | 25, 40 | hexagonal grid dimensions |
| `spacing` | 1.0 | hex grid spacing |
| `k_normal` / `k_fibrosis` | 1.0 / 50.0 | initial spring stiffness classes |
| `alpha` | 0.6 | target fraction of springs seeded fibrotic |
| `random_walk_length` | 300 | steps per fibrotic cluster random walk |
| `l0_ratio` | 0.6 | `l0 = l0_ratio * spacing / sqrt(3)` (exact 0606 expression) |

Fibrotic seeding repeatedly starts a random walk from a healthy spring and
marks springs along the walk as `k_fibrosis` until the fibrotic fraction
reaches `alpha` (`core.py: make_initial_state`).

## Mechanical relaxation (`core.py: mechanical_equilibrate`)

| Parameter | Value | Meaning |
|---|---|---|
| `mu` | 0.01 | initial step size; `*1.1` on accepted downhill step, `*0.9` on rejection |
| `temperature` | 1.0 | Metropolis temperature; decays `*0.99` per iteration |
| `tolerance` | 1e-7 | relative energy-change threshold counted as stagnation |
| `stagnation_window` | 20 | consecutive stagnant iterations to declare convergence |
| `preconvergence_iterations` | 2000 | mechanical iterations for the shared initial state |
| `agent_update_interval` | 1500 | mechanical iterations run per outer remodeling phase |
| `agent_total_iterations` | 600 | number of outer remodeling phases per run |

One uniform random draw is always consumed per mechanical iteration
(accepted or not), so paired control/softening branches share the same
phase/iteration-indexed Metropolis noise field even if their energy
trajectories diverge (`core.py: mechanical_equilibrate`).

## Remodeling and biological activation (`core.py: update_activation_and_density`)

| Parameter | Value | Meaning |
|---|---|---|
| `At` | 10.0 | stiffness scale in the activation function `a_k` |
| `A0` | 1.0 | reference remodeling-multiplier scale |
| `epsilon_s` | 1.0 | strain scale in the activation function `a_epsilon` |
| `beta`, `gamma` | 3.0, 1.0 | Hill-function exponents for `a_k`, `a_epsilon` |
| `w1`, `w2` | 1.0, 1.0 | activation-function weights |
| `memory_factor` | 1.0 | activation update memory weight |
| `area_update_scale` | 0.2 | scales the remodeling increment `P * activation * D` |
| `p1`, `D0`, `Dmax` | 0.01, 0.3, 3.0 | agent-density update parameters |
| `E_initial` | 1.0 | initial Young's-modulus-like state variable |
| `remodeling_law` | `"positive_degradation"` | **production** update rule; see below |

`A` is a positive spring remodeling multiplier (`k = E * A`), not a
geometric area — it is initialized from spring stiffness after fibrosis
seeding. Two update rules exist:

- `legacy_additive` — the exact archival 0606 additive update
  `A <- A + increment`. It is **not used in production**; it is retained
  only for forensic comparison and can drive `A` non-positive, which makes
  `k = E*A` mechanically invalid.
- `positive_degradation` (**production default**) — growth (`increment >= 0`)
  keeps the archival additive step; degradation (`increment < 0`) instead
  applies a local positive multiplicative decay,
  `A <- A * exp(increment / A_ref)` with `A_ref = k_fibrosis` or `k_normal`
  per spring, so `A` stays strictly positive without a floor or new
  parameter.

## FF-like softening intervention (`SofteningConfig`, applied in `core.py: apply_softening`)

| Parameter | Value | Meaning |
|---|---|---|
| `threshold` | 5.0 | spring stiffness `k` above which a spring becomes eligible to soften, once |
| `magnitude` | 0.05 | fractional drop applied to `E` at the trigger instant |
| `duration` | 150 | outer phases the softened `E` is held before recovery begins |
| `recovery_rate` | 0.003 | per-phase geometric recovery rate of `E` toward its recorded pre-event value |
| `recovery_tolerance` | 1e-8 | absolute tolerance for declaring recovery complete |
| `minimum_E` | 0.1 | floor on `E` after a softening event |

Each spring can soften at most once. At trigger, `E` drops immediately from
its current value to `max(minimum_E, pre_E * (1 - magnitude))`, is held for
`duration` phases, then recovers geometrically toward the recorded pre-event
`E` under a tolerance-based (not floating-point-equality) termination rule.
Softening never modifies `A` directly; it only changes `k` via `k = E * A`.
`magnitude = 0.05` is the visible 0606 notebook value carried forward for
continuity, but it was inert in that notebook (the softening path was never
exercised) and is **not a biological calibration**.
`--config smoke` intentionally uses `magnitude=0.5` on a tiny network purely
to make the software intervention unambiguous in tests; it is not a
biological result (see `README.md`, `scripts/run_smoke_test.py`).

## Output cadence (`ModelConfig`)

- `output_interval = 1`: a compact per-phase timeseries row is written every
  outer phase (`pulmonary_fibrosis_model/io.py: state_summary`).
- `spatial_output_interval = 25` (baseline) / `1` (smoke): cadence, in outer
  phases, of full spatial NPZ snapshots (`io.py: write_spatial_snapshot`).

## Morphometry: raw stored quantity vs. current standard-reference definition

Two distinct area-normalization definitions exist in this codebase. They
must not be conflated.

### 1. Raw stored self-normalized quantity (still saved by `io.py`)

`pulmonary_fibrosis_model/io.py: write_spatial_snapshot` still saves, for
every polygon `i` at every recorded phase `t`:

```
normalized_cell_area[i] = cell_area[i] / initial_cell_areas[i]
```

i.e. **`A_i(t) / A_i(0)`** — each polygon normalized against *its own*
phase-0 area. The compact per-phase table (`io.py: state_summary`) saves the
mean/median of this same per-polygon self-normalized quantity
(`mean_normalized_cell_area`, `median_normalized_cell_area`). This is a raw
model output, not itself the collaborator-facing definition.

### 2. Current standard-reference collaborator/publication morphometry

`analysis/standard_reference.py: reference_areas` defines, once per seed:

```
A_ref(seed) = median(phase-0 control polygon areas for that seed)
```

using all phase-0 polygons in that seed's control initial geometry (no
robust polygon-level healthy/non-fibrotic label exists to select a subset).
Every polygon `i` at phase `t` in that seed is then normalized as

```
A_i(t) / A_ref(seed)
```

— one fixed scalar denominator shared by all polygons and phases within a
seed, rather than each polygon's own initial area. Paired phase-0 control
and softening geometries are verified bit-identical before `A_ref` is used
(`standard_reference.py: reference_areas`).

**Current standard-reference figures, videos, and the collaborator package
use definition 2 (`A_i(t) / A_ref(seed)`), not definition 1.** This is wired
through `analysis/morphometry.py: metrics(..., reference_area=...)`,
`analysis/aggregate.py: build(..., reference_by_seed=...)`, and
`analysis/collaborator_package.py` (`render_region_frame`, `render_seed`),
all invoked from `analysis/build_collaborator_package.py` and
`analysis/run_analysis.py`, both of which call
`standard_reference.reference_areas` and pass its result through explicitly.

**Legacy path warning:** `analysis/generate_simple_main_figures.py` reads
`snapshot["normalized_cell_area"]` directly, i.e. definition 1
(`A_i(t) / A_i(0)`, per-polygon self-normalization), and does **not** use
`A_ref(seed)`. Its output directory and figures are legacy and must not be
mistaken for the standard-reference workflow; see the warning at the top of
that file.

### Operational region thresholds

Both `analysis/morphometry.py: metrics()` and
`analysis/collaborator_package.py: classify()` use the same fixed operational
thresholds on whichever normalized quantity is supplied:

- Contracted / low-area: normalized area `< 0.75`
- Preserved-like: `0.75 <= normalized area <= 1.25`
- Enlarged / remodeled: normalized area `> 1.25`

These are explicit operational, model-based definitions for collaborator
review; they are not validated against histology.
