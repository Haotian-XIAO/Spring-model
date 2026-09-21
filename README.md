# Spring-model

Reproducible pulmonary-fibrosis spring-network model with a paired FF-like
softening intervention, plus the analysis pipeline that produces the current
N=100 standard-reference collaborator package.

## Where to start

- **`docs/model_parameters.md`** — Human-readable reference for the active
  model parameters and their meaning.
- **`pulmonary_fibrosis_model/config.py`** — Source of truth for active
  numerical parameter values.
- **`pulmonary_fibrosis_model/core.py`** — Core spring-network mechanics,
  remodeling law, and softening mechanism.
- **`pulmonary_fibrosis_model/experiment.py`** — Paired no-softening /
  softening experiment design and RNG handling.
- **`analysis/standard_reference.py`** — Defines the reference alveolar area
  used by the current standard-reference morphometry.
- **`analysis/morphometry.py`** — Defines morphometric quantities and
  Contracted / Preserved / Enlarged classification.
- **`analysis/build_collaborator_package.py`** — Entry point used to generate
  the current N=100 collaborator-facing figures and videos.

## Repository map

| Path | Purpose |
|---|---|
| `pulmonary_fibrosis_model/config.py` | Active numerical parameter values (`ModelConfig`, `SofteningConfig`, `smoke_config`) |
| `pulmonary_fibrosis_model/core.py` | Geometry, mechanical relaxation, remodeling law, and FF-like softening state machine |
| `pulmonary_fibrosis_model/experiment.py` | Paired control/softening run driver and `SeedSequence`-based RNG design |
| `pulmonary_fibrosis_model/io.py` | Raw per-phase timeseries and spatial NPZ snapshot output (includes the raw self-normalized area) |
| `scripts/run_pairs.py` | Canonical CLI entry point for running paired experiments |
| `analysis/standard_reference.py` | Defines `A_ref(seed)`, the current standard-reference area normalization |
| `analysis/morphometry.py` | Per-polygon morphometric quantities and the Contracted / Preserved / Enlarged classification |
| `analysis/aggregate.py` | Aggregates per-seed metrics into ensemble-level tables using the standard-reference definition |
| `analysis/run_analysis.py` | CLI driver that runs the standard-reference analysis over saved experiment outputs |
| `analysis/plotting.py` | Shared plotting style/helpers used by the standard-reference figures |
| `analysis/collaborator_package.py` | Renders region-classification frames and videos for the collaborator package |
| `analysis/build_collaborator_package.py` | Assembles the full collaborator-facing package (figures, videos, tables) |
| `analysis/reporting.py` | Writes the text/markdown summary report for an ensemble |
| `analysis/load_results.py` | Discovers and loads saved per-seed experiment outputs |
| `analysis/generate_simple_main_figures.py` | **Legacy/deprecated.** Uses the old self-normalized `A_i(t)/A_i(0)` definition; NOT used for current standard-reference results |
| `docs/model_parameters.md` | Human-readable reference for all active model parameters |
| `archive/legacy/` | Historical pre-refactor notebooks, preserved byte-for-byte; provenance only |
| `tests/` | Unit test suite validating the model and analysis code |

## For biological collaborators

- Model parameters -> `docs/model_parameters.md`
- Mechanical softening -> `pulmonary_fibrosis_model/core.py`
- Paired simulation protocol -> `pulmonary_fibrosis_model/experiment.py`
- Definition of preserved alveoli -> `analysis/standard_reference.py` and `analysis/morphometry.py`

The current collaborator/publication analysis uses the standard-reference
morphometry. Files under `archive/legacy/` and
`analysis/generate_simple_main_figures.py` are retained for historical
provenance and should not be used to reproduce the current results.

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m unittest discover -s tests -v
```

Run one small paired verification (a software test, not a biological
result):

```bash
python scripts/run_smoke_test.py
```

See `PIPELINE.md` for the full paired-experiment workflow and RNG design, and
`analysis/README.md` for turning saved paired outputs into figures, tables,
and the collaborator package.

## Model summary

The model is a hexagonal spring network with fibrotic seeding, annealed
mechanical relaxation, biological remodeling, and an optional FF-like
softening intervention on individual springs. Full current parameter values
live in `docs/model_parameters.md`, sourced only from the current code (never
from the legacy notebooks). For the two area-normalization definitions in
the codebase and which one is current, see `docs/model_parameters.md` and
the warning above.
