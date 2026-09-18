# Spring-model

Reproducible pulmonary-fibrosis spring-network model with a paired FF-like
softening intervention, plus the analysis pipeline that produces the current
N=100 standard-reference collaborator package.

## Where to look

| Topic | Location |
|---|---|
| Current model definition, equations, and code | `pulmonary_fibrosis_model/` (`config.py`, `core.py`, `experiment.py`, `io.py`) |
| Active parameters (network, mechanics, remodeling, softening) | `docs/model_parameters.md` |
| Softening (FF-like) mechanism | `docs/model_parameters.md#ff-like-softening-intervention-softeningconfig-applied-in-corepy-apply_softening`, implemented in `pulmonary_fibrosis_model/core.py: apply_softening` |
| Paired experiment / RNG design | `PIPELINE.md`, `scripts/run_pairs.py` |
| Standard-reference morphometry (current publication definition) | `docs/model_parameters.md`, `analysis/standard_reference.py` |
| Analysis and collaborator-package pipeline | `analysis/README.md`, `analysis/run_analysis.py`, `analysis/build_collaborator_package.py` |
| Historical/legacy notebooks (not part of the validated pipeline) | `archive/legacy/` |

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

## Model and morphometry summary

The model is a hexagonal spring network with fibrotic seeding, annealed
mechanical relaxation, biological remodeling, and an optional FF-like
softening intervention on individual springs. Full current parameter values
live in `docs/model_parameters.md`, sourced only from the current code (never
from the legacy notebooks).

Two area-normalization definitions exist and must not be conflated:

- the raw stored per-polygon self-normalized quantity `A_i(t) / A_i(0)`
  (still saved by `pulmonary_fibrosis_model/io.py`), and
- the current standard-reference definition
  `A_ref(seed) = median phase-0 control polygon area`, then
  `A_i(t) / A_ref(seed)` (`analysis/standard_reference.py`).

**Current standard-reference figures and videos use the second definition.**
`analysis/generate_simple_main_figures.py` is a legacy path that still uses
the first (self-normalized) definition; it is clearly marked as legacy in
that file and must not be mistaken for the standard-reference workflow. See
`docs/model_parameters.md` for the full explanation.

## Historical notebooks

`archive/legacy/` preserves, byte-for-byte, the pre-refactor notebooks that
this repository's tested Python pipeline reconstructs and extends. They are
provenance only — see `archive/legacy/README.md`.
