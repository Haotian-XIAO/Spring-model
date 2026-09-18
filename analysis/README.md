# Paired-result analysis

This package analyzes saved paired simulation outputs only. It never imports or
changes the simulation equations.

## Commands

```bash
python analysis/run_analysis.py \
  --input-root experiment_results/pilot_m005_n100/batch_01 \
  --output-root analysis_results/pilot_m005_n100/batch_01

python analysis/run_analysis.py \
  --input-root experiment_results/pilot_m005_n100 \
  --output-root analysis_results/pilot_m005_n100/ensemble
```

Run these in the `spring-model` environment, which provides NumPy and
Matplotlib. Seed discovery is recursive, so the same command handles a single
batch or a multi-batch ensemble root.

## Morphometry definitions

Continuous polygon-derived quantities are primary: cell area, normalized cell
area, quantiles, coefficient of variation, interquartile range, Gini
coefficient, and mean absolute normalized-area contrast between neighboring
cells.

The optional operational preserved-alveolus definition is
`lower <= final/reference cell area <= upper`, defaulting to 0.75–1.25. It is
not claimed to be a histologically validated threshold. Change it explicitly
with `--preserved-lower` and `--preserved-upper` and compare sensitivity.
Connectivity summaries count connected components of cells satisfying that
operational definition and report the fraction belonging to the largest
component. The fixed hex-cell adjacency is used; these are candidate spatial
fragmentation descriptors, not evidence of literal septal rupture.

Model variable `A` is always labelled **spring remodeling multiplier** and is
kept distinct from geometric cell/alveolar area.
