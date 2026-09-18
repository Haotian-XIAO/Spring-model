# Legacy notebooks (historical provenance only)

The notebooks in this directory are **not part of the current validated
pipeline**. They are preserved byte-for-byte as historical provenance for the
archival "0606" agent-based model that `pulmonary_fibrosis_model/` and
`scripts/run_pairs.py` reconstruct and extend in versioned, tested Python code.

Do not import, run, or cite these notebooks as the source of current
parameters, equations, or results. They predate:

- the paired control/softening (FF-like) experiment design;
- the reproducible `SeedSequence`-based RNG design;
- the numerical-instability guards in `pulmonary_fibrosis_model/core.py`;
- the `positive_degradation` remodeling law;
- the standard-reference morphometry definition in
  `analysis/standard_reference.py`.

Files:

- `Main_agent_young.ipynb`, `Main_agent_young copy.ipynb` — early single-branch
  agent-based model notebooks, pre-softening.
- `Main_agent_young_withsoftening_0606.ipynb`,
  `Main_agent_young_withsoftening_0606_with table.ipynb` — the "0606" notebooks
  that introduced the softening intervention structure later reconstructed
  (with corrections) in `pulmonary_fibrosis_model/`.
- `Random_walk.ipynb`, `Random_walk_progressive.ipynb` — exploratory notebooks
  for the fibrotic-cluster random-walk seeding used by the initial condition.

For current model definitions, parameters, and morphometry, see
`docs/model_parameters.md`, `README.md`, and `PIPELINE.md` at the repository
root.
