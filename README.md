# Spring-model

Spring-network model of pulmonary fibrosis for paired simulations with and without FF-like softening.

## Scientific entry points

- `docs/model_parameters.md` — Current model parameters and their meaning
- `pulmonary_fibrosis_model/core.py` — Mechanical model, remodeling, and softening
- `pulmonary_fibrosis_model/experiment.py` — Paired simulation protocol and randomization
- `analysis/standard_reference.py` — Reference alveolar area definition
- `analysis/morphometry.py` — Morphometric indicators and alveolar classification

## Repository map

| Path | Purpose |
|---|---|
| `pulmonary_fibrosis_model/config.py` | Numerical parameter values |
| `pulmonary_fibrosis_model/core.py` | Mechanics, remodeling, and softening |
| `pulmonary_fibrosis_model/experiment.py` | Paired simulation protocol |
| `pulmonary_fibrosis_model/io.py` | Simulation output |
| `scripts/run_pairs.py` | Production simulation entry point |
| `analysis/standard_reference.py` | Reference-area definition |
| `analysis/morphometry.py` | Morphometric analysis |
| `analysis/aggregate.py` | Ensemble aggregation |
| `analysis/run_analysis.py` | Analysis entry point |
| `analysis/plotting.py` | Figure generation |
| `analysis/collaborator_package.py` | Morphology figures and videos |
| `analysis/build_collaborator_package.py` | N=100 collaborator package |
| `docs/model_parameters.md` | Parameter reference |
| `archive/legacy/` | Historical notebooks |
| `tests/` | Lightweight tests |

## Current morphometry

The current N=100 analysis uses:

`A_ref = median alveolar area at phase 0`

Each alveolus is classified from its current area relative to `A_ref`.

The older self-normalized definition `A_i(t)/A_i(0)` is retained only for historical reference.

## Reproducibility

- `PIPELINE.md` — Simulation workflow
- `analysis/README.md` — Analysis workflow
