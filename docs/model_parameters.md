# Model parameters and morphometry definitions

**Source of truth:** this document is derived only from the current code in
`pulmonary_fibrosis_model/` and `analysis/` at this commit
(`pulmonary_fibrosis_model/config.py`, `core.py`, `experiment.py`, `io.py`;
`analysis/standard_reference.py`, `morphometry.py`). It does not derive
anything from the historical notebooks in `archive/legacy/`, which are
provenance-only and may disagree with the values below.

Each table row gives: code name, symbol (if the code/comments use one),
current nominal value, units (if known — this is a nondimensional model, so
most quantities are dimensionless "model units"), a short meaning, and the
source file that defines or consumes it.

## 1. Geometry (`ModelConfig`, `pulmonary_fibrosis_model/config.py`, `core.py: build_geometry`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `nx` | — | 25 | cells | hexagonal grid width | `config.py` |
| `ny` | — | 40 | cells | hexagonal grid height | `config.py` |
| `spacing` | — | 1.0 | model length | hex grid center-to-center spacing (`min_diam` of `create_hex_grid`) | `config.py`, `core.py: build_geometry` |
| `l0_ratio` | — | 0.6 | dimensionless | ratio defining the spring rest length `l0` | `config.py` |
| `l0` (derived property) | $l_0$ | `l0_ratio * spacing / sqrt(3)` ≈ 0.346 | model length | spring rest length used everywhere strain/stress/energy are computed; exact 0606 expression | `config.py: ModelConfig.l0` |

## 2. Initial fibrosis (`core.py: make_initial_state`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `k_normal` | $k_{healthy}$ | 1.0 | model stiffness | initial stiffness of a non-fibrotic spring | `config.py` |
| `k_fibrosis` | $k_{fib}$ | 50.0 | model stiffness | initial stiffness of a fibrotic spring | `config.py` |
| `alpha` | $\alpha$ | 0.6 | fraction (0–1) | target fraction of springs seeded fibrotic before preconvergence | `config.py` |
| `random_walk_length` | — | 300 | steps | number of connected springs marked fibrotic per random-walk cluster seed | `config.py` |

Fibrotic seeding repeatedly starts a random walk from a healthy spring (chosen
via `biology_rng`) and marks springs along the walk `k_fibrosis` until the
fibrotic fraction reaches `alpha`. `A` is subsequently initialized from
`spring_constants` (i.e. `A = k_normal` or `k_fibrosis` per spring) — this is
the archival 0606 initialization, preserved exactly (`core.py:
make_initial_state`).

## 3. Mechanics (`core.py: mechanical_equilibrate`, `compute_forces`, `compute_total_energy`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `mu` | $\mu$ | 0.01 (initial) | model step size | annealed gradient-descent step size; `*1.1` on an accepted downhill step, `*0.9` on a rejection | `config.py`, `core.py: mechanical_equilibrate` |
| `temperature` | $T$ | 1.0 (initial) | model energy | Metropolis temperature for uphill moves; decays `*0.99` every mechanical iteration | `config.py`, `core.py: mechanical_equilibrate` |
| spring energy | $E_{el} = \tfrac12\sum k (\ell-\ell_0)^2$ | — | model energy | total elastic energy driving relaxation | `core.py: compute_total_energy` |
| strain | $\varepsilon = (\ell-\ell_0)/\ell_0$ | — | dimensionless | per-spring strain from current length vs. `l0` | `core.py: compute_strain_and_stress` |
| stress | $\sigma = k\varepsilon$ | — | model stress | per-spring stress | `core.py: compute_strain_and_stress` |
| `k` (`spring_constants`) | $k$ | `k = E * A` | model stiffness | current spring stiffness, set from `E` and `A` at every phase (`core.py: apply_softening`) | `core.py` |

## 4. Remodeling and biological activation (`core.py: update_activation_and_density`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `At` | $A_t$ | 10.0 | model stiffness scale | stiffness scale in the activation term `a_k` | `config.py` |
| `A0` | $A_0$ | 1.0 | model scale | reference remodeling-multiplier scale (also used in `activation_offset`) | `config.py` |
| `epsilon_s` | $\varepsilon_s$ | 1.0 | dimensionless | strain scale in the activation term `a_ε` | `config.py` |
| `beta` | $\beta$ | 3.0 | dimensionless | Hill-function exponent for `a_k` | `config.py` |
| `gamma` | $\gamma$ | 1.0 | dimensionless | Hill-function exponent for `a_ε` | `config.py` |
| `w1` | $w_1$ | 1.0 | dimensionless | weight of `a_ε` in the activation update | `config.py` |
| `w2` | $w_2$ | 1.0 | dimensionless | weight of `a_k` in the activation update | `config.py` |
| `memory_factor` | — | 1.0 | dimensionless | activation update memory weight (1.0 = no memory of prior activation) | `config.py` |
| `area_update_scale` | $P$ | 0.2 | dimensionless | scales the remodeling increment `P * activation * D` | `config.py` |
| `remodeling_law` | — | `"positive_degradation"` (production) | — | selects the `A` update rule; see below | `config.py` |
| `p1` | — | 0.01 | dimensionless | agent-density update parameter | `config.py` |
| `D0` | — | 0.3 | model density | minimum/reference agent density | `config.py` |
| `Dmax` | — | 3.0 | model density | maximum agent density (`D` is clipped to `[D0, Dmax]`) | `config.py` |
| `E_initial` | $E_0$ | 1.0 | model stiffness | initial value of `E` for every spring | `config.py` |
| `A` | $A$ | initialized to `k_normal`/`k_fibrosis` | model remodeling multiplier | positive spring remodeling multiplier, **not** a geometric area (`k = E * A`) | `core.py` |

Two `A` update rules exist:

- `legacy_additive` — the exact archival 0606 additive update `A <- A +
  increment`. **Not used in production**; retained only for forensic
  comparison. It can drive `A` non-positive, which makes `k = E*A`
  mechanically invalid.
- `positive_degradation` (**production default**) — growth
  (`increment >= 0`) keeps the archival additive step; degradation
  (`increment < 0`) instead applies a local positive multiplicative decay,
  `A <- A * exp(increment / A_ref)` with `A_ref = k_fibrosis` or `k_normal`
  per spring, so `A` stays strictly positive without a floor or new
  parameter.

## 5. FF-like softening intervention (`SofteningConfig`, applied in `core.py: apply_softening`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `threshold` | — | 5.0 | model stiffness | spring stiffness `k` above which a spring becomes eligible to soften, once | `config.py` |
| `magnitude` | — | 0.05 | fraction (0–1) | fractional drop applied to `E` at the trigger instant | `config.py` |
| `duration` | — | 150 | outer phases | number of outer phases the softened `E` is held before recovery begins | `config.py` |
| `recovery_rate` | — | 0.003 | fraction / phase | per-phase geometric recovery rate of `E` toward its recorded pre-event value | `config.py` |
| `recovery_tolerance` | — | 1e-8 | model stiffness (absolute) | absolute tolerance for declaring recovery complete | `config.py` |
| `minimum_E` | — | 0.1 | model stiffness | floor on `E` after a softening event | `config.py` |

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

## 6. Stochastic / RNG design (`pulmonary_fibrosis_model/experiment.py`)

Every seed's RNG state derives from one integer master seed via
`numpy.random.SeedSequence`, so runs are exactly reproducible. No parameter
here has a "nominal value" in the same sense as above — this section
documents the stream derivation itself.

| Stream (code name in `_stream_metadata`) | Derivation | Consumed by | Source file |
|---|---|---|---|
| `biology_initialization` | `SeedSequence(master_seed).spawn(5)[0]` | fibrotic-cluster random walk (`core.py: make_initial_state`) | `experiment.py: generate_initial_state` |
| `shared_preconvergence_mechanics` | `spawn(5)[1]` | shared preconvergence `mechanical_equilibrate` call producing the common initial state | `experiment.py: generate_initial_state` |
| `paired_mechanics_common_random_numbers` | `spawn(5)[2]` | control and softening branches when `common_random_numbers=True` (used by `run_pair`, the production paired run) | `experiment.py: run_pair` |
| `control_mechanics_independent` | `spawn(5)[3]` | control branch of the RNG-null test (`run_null_pair`, independent OFF/OFF mechanics) | `experiment.py: run_null_pair` |
| `softening_mechanics_independent` | `spawn(5)[4]` | second branch of the RNG-null test | `experiment.py: run_null_pair` |
| per-phase Metropolis stream | `SeedSequence([stream_seed, phase])` | one reproducible Metropolis-uniform stream per outer phase, re-derived from the branch's stream seed and the phase index | `experiment.py: _phase_mechanics_rng` |

For the production paired run (`run_pair`), both the control and softening
branches use the **same** `paired_mechanics_common_random_numbers` stream
seed, so at every outer phase both conditions draw from an identical
phase-indexed Metropolis stream: one uniform draw is consumed per mechanical
iteration regardless of whether the trial step is accepted (`core.py:
mechanical_equilibrate`), so the two branches share the same
phase/iteration-indexed noise field even if their energy trajectories and
iteration counts diverge. This is deliberate common-random-numbers pairing,
so that downstream differences are attributable to the FF-like softening
intervention rather than to solver RNG noise. `run_null_pair` instead uses
two *independent* streams (OFF vs. OFF) to quantify how much of a paired
difference could arise from mechanics RNG noise alone, with no softening
involved.

## 7. Numerical solver (`core.py: mechanical_equilibrate`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `tolerance` | — | 1e-7 | relative | relative energy-change threshold counted as one stagnant iteration | `config.py` |
| `stagnation_window` | — | 20 | iterations | consecutive stagnant iterations required to declare convergence | `config.py` |
| `preconvergence_iterations` | — | 2000 | mechanical iterations | iteration budget for the one shared common initial state | `config.py` |
| `agent_update_interval` | — | 1500 | mechanical iterations | mechanical-iteration budget run per outer remodeling phase | `config.py` |

The solver is an annealed force-relaxation/Metropolis scheme
(`mechanical_equilibrate`): each iteration computes nodal forces, proposes a
`mu`-scaled displacement, accepts it outright if energy decreases, otherwise
accepts uphill moves with Metropolis probability `exp(-ΔE/temperature)`.
Every array and scalar involved is checked for finiteness at each step, and
`NumericalInstabilityError` is raised with a detailed diagnostic before any
non-finite state could otherwise be accepted or written (`core.py:
_raise_if_nonfinite`, `_raise_if_nonpositive_stiffness`).

## 8. Simulation duration and output cadence (`config.py`, `experiment.py: simulate`)

| Code name | Symbol | Value | Units | Meaning | Source file |
|---|---|---|---|---|---|
| `agent_total_iterations` | — | 600 | outer phases | number of outer remodeling phases per run (the full production run length) | `config.py` |
| `output_interval` | — | 1 | phases | cadence, in outer phases, of the compact per-phase timeseries row | `config.py`, `io.py: state_summary` |
| `spatial_output_interval` | — | 25 (baseline) / 1 (smoke) | phases | cadence, in outer phases, of full spatial NPZ snapshots | `config.py`, `io.py: write_spatial_snapshot` |

One production paired run therefore executes 600 outer phases, each running
up to `agent_update_interval` (1500) mechanical iterations, preceded by one
shared `preconvergence_iterations` (2000) mechanical-iteration run to build
the common initial state. `scripts/run_smoke_test.py` uses the much smaller
`smoke_config()` (see `config.py: smoke_config`) purely as a software
end-to-end check, not a biological result.

## 9. Morphometry: raw stored quantity vs. current standard-reference definition

Two distinct area-normalization definitions exist in this codebase. They
must not be conflated.

### 9.1 Raw stored self-normalized quantity (still saved by `io.py`)

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

### 9.2 Current standard-reference collaborator/publication morphometry

`analysis/standard_reference.py: reference_areas` defines, once per seed:

```
A_ref(seed) = median(phase-0 control polygon areas for that seed)
```

using all phase-0 polygons in that seed's control initial geometry (no
robust polygon-level healthy/non-fibrotic label exists to select a subset).
Every polygon `i` at phase `t` in that seed is then normalized as

```
r_i(t) = A_i(t) / A_ref(seed)
```

— one fixed scalar denominator shared by all polygons and phases within a
seed, rather than each polygon's own initial area. Paired phase-0 control
and softening geometries are verified bit-identical before `A_ref` is used
(`standard_reference.py: reference_areas`).

**Current standard-reference figures, videos, and the collaborator package
use definition 9.2 (`r_i(t) = A_i(t) / A_ref(seed)`), not definition 9.1.**
This is wired through `analysis/morphometry.py: metrics(...,
reference_area=...)`, `analysis/aggregate.py: build(...,
reference_by_seed=...)`, and `analysis/collaborator_package.py`
(`render_region_frame`, `render_seed`), all invoked from
`analysis/build_collaborator_package.py` and `analysis/run_analysis.py`, both
of which call `standard_reference.reference_areas` and pass its result
through explicitly.

**Legacy path warning:** `analysis/generate_simple_main_figures.py` reads
`snapshot["normalized_cell_area"]` directly, i.e. definition 9.1
(`A_i(t) / A_i(0)`, per-polygon self-normalization), and does **not** use
`A_ref(seed)`. Its output directory and figures are legacy and must not be
mistaken for the standard-reference workflow; see the warning at the top of
that file.

### 9.3 Operational region thresholds

Both `analysis/morphometry.py: metrics()` and
`analysis/collaborator_package.py: classify()` use the same fixed operational
thresholds on whichever normalized quantity is supplied (`r`, per section 9.1
or 9.2):

| Region | Threshold |
|---|---|
| Contracted | `r < 0.75` |
| Preserved | `0.75 <= r <= 1.25` |
| Enlarged | `r > 1.25` |

These are explicit operational, model-based definitions for collaborator
review; they are not validated against histology.

**The current N=100 standard-reference figures, videos, and the
`collaborator_package/FF_softening_pilot_N100_standard_ref` delivery package
use definition 9.2 (`A_ref(seed)`-normalized, `analysis/standard_reference.py`).**
Older unsuffixed collaborator-package outputs
(`FF_softening_pilot_N100`, `FF_softening_pilot_N100_TEST_3seeds`) predate the
2026-09-18 standard-reference fix and use definition 9.1
(`A_i(t) / A_i(0)`, per-polygon self-normalization) instead; they must not be
cited as the current standard-reference result.
