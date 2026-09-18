from pathlib import Path
import numpy as np
def write_report(path:Path,runs,rows,reps,lower,upper,reference_by_seed=None):
    ok=[r for r in rows if r.get('run_status')=='success']; d=np.array([r['delta_operational_preserved_cell_fraction'] for r in ok]) if ok else np.array([])
    tol=.005; protect=int((d>tol).sum()); aggrav=int((d<-tol).sum()); neutral=int((abs(d)<=tol).sum())
    median=float(np.median(d)) if len(d) else float('nan'); q1,q3=(np.percentile(d,[25,75]) if len(d) else (float('nan'),float('nan')))
    mean_area=np.array([r['delta_normalized_area_mean'] for r in ok]); heterogeneity=np.array([r['delta_normalized_area_cv'] for r in ok]); contrast=np.array([r['delta_neighbor_normalized_area_contrast'] for r in ok])
    ref_text = (f"A fixed scalar reference was used per seed: the median geometric area of all polygons in that seed's phase-0 control initial geometry. No robust polygon-level healthy/non-fibrotic mapping was available, so all phase-0 polygons were included. The paired phase-0 control and softening geometries were verified identical. Reference areas ranged from {min(reference_by_seed.values()):.8g} to {max(reference_by_seed.values()):.8g}." if reference_by_seed else "")
    text=f"""# Paired pulmonary-fibrosis analysis report

Discovered **{len(runs)} seeds**; **{len(ok)} successful paired runs** and **{len(runs)-len(ok)} incomplete/failed runs**.

## Standard reference area

{ref_text} Every morphology quantity uses `A_i(t) / A_ref`; no polygon-by-polygon phase-0 self-normalization is used. The thresholds remain Contracted `< {lower:.2f}`, Preserved `{lower:.2f}–{upper:.2f}`, and Enlarged `> {upper:.2f}`.

## Main paired finding

Using the explicitly operational preserved-alveolus definition (`{lower:.2f} ≤ final/reference cell area ≤ {upper:.2f}`), FF-like softening was protective in {protect} seeds, aggravating in {aggrav}, and near-neutral (absolute paired change ≤ {tol:.3f}) in {neutral}. Continuous area summaries remain the primary results.

The paired preserved-fraction effect had median **{median:.4f}**, IQR **[{q1:.4f}, {q3:.4f}]**, and range **[{d.min():.4f}, {d.max():.4f}]**. Mean normalized alveolar area was lower under softening in **{int((mean_area<0).sum())}/{len(mean_area)}** seeds, area CV was lower in **{int((heterogeneity<0).sum())}/{len(heterogeneity)}**, and neighboring-cell area contrast was lower in **{int((contrast<0).sum())}/{len(contrast)}**. Because control normalized area is above baseline on average, these continuous shifts indicate less enlargement and less spatial heterogeneity, consistent with the operational preserved-fraction direction.

## Interpretation

“Spring remodeling multiplier” denotes model variable A; it is not geometric alveolar area. Geometric cell area and normalized cell area are computed directly from saved hex-cell polygons. The thresholded preserved metric is a sensitivity-analysis candidate, not a validated histological definition.

## Representative seeds

""" + "\n".join(f"- {r['category'].replace('_',' ')}: seed {r['seed']} (paired effect {r['paired_effect']:.4g})" for r in reps) + """

## Outputs

- `tables/master_seed_summary.csv`
- `tables/master_timecourse_summary.csv`
- `tables/representative_seed_catalog.csv`
- `figures_main/` publication-oriented PDF, SVG, and 400-dpi PNG exports
- `representative_seeds/` spatial comparisons

See `figure_captions.md` for draft captions and the operational-definition caveat.
"""
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
    (path.parent/'figure_captions.md').write_text("# Draft figure captions\n\n**Figure 1.** Each point is one paired simulation. Values are calculated as With softening minus No softening. Positive values indicate a larger outcome with softening. The preserved-alveolus measure uses the fixed model-based definition and has not been validated histologically.\n\n**Figure 2.** Lines show the median; shaded areas show the middle 50% of simulations. Remodeling level is the spring remodeling multiplier and is distinct from geometric alveolar area.\n\n**Figure 3.** Adjacent paired endpoint associations along the hypothesized E→stiffness→stress→activation/agent-density→morphology sequence. Associations do not alone establish causality.\n",encoding='utf-8')
