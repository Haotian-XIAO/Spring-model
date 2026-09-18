from __future__ import annotations

import numpy as np

from .load_results import SeedRun, load_snapshot, snapshot_path


def reference_areas(runs: list[SeedRun]) -> dict[int, float]:
    """Return one fixed phase-0 standard alveolus area per seed.

    Existing outputs do not contain a robust polygon-level healthy/non-fibrotic
    label.  Therefore all phase-0 polygons in the control initial geometry are
    used, with the median as the representative unstretched normal area.
    """
    refs: dict[int, float] = {}
    for run in runs:
        if run.status != "success":
            continue
        control_path = snapshot_path(run, "control", 0)
        softening_path = snapshot_path(run, "softening", 0)
        if control_path is None or softening_path is None:
            raise FileNotFoundError(f"seed {run.seed}: missing phase-0 spatial snapshot")
        control = load_snapshot(control_path)
        softening = load_snapshot(softening_path)
        if not np.allclose(control["cell_area"], softening["cell_area"], rtol=0, atol=1e-12):
            raise ValueError(f"seed {run.seed}: paired phase-0 geometries are not identical")
        refs[run.seed] = float(np.median(np.asarray(control["cell_area"], float)))
    return refs
