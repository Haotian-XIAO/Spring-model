from __future__ import annotations
import csv, json, os, re, warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np

@dataclass(frozen=True)
class SeedRun:
    seed: int
    path: Path
    status: str
    problems: tuple[str, ...]

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def discover_seeds(root: Path) -> list[SeedRun]:
    if not root.exists():
        raise FileNotFoundError(f"Analysis input root does not exist: {root}")
    found: dict[int, Path] = {}
    visited: set[Path] = set()
    for directory, names, _ in os.walk(root, followlinks=True):
        path = Path(directory)
        resolved = path.resolve()
        if resolved in visited:
            names[:] = []
            continue
        visited.add(resolved)
        m = re.fullmatch(r"seed_(\d+)", path.name)
        if m and (path / "pair_metadata.json").exists():
            seed = int(m.group(1))
            if seed in found and found[seed].resolve() != resolved:
                warnings.warn(f"duplicate seed {seed}: using {found[seed]}, ignoring {path}")
            else:
                found[seed] = path
            names[:] = []
    runs=[]
    for seed,path in sorted(found.items()):
        missing=[]
        for name in ("pair_metadata.json","pair_timeseries.csv"):
            if not (path/name).exists(): missing.append(name)
        for condition in ("control","softening"):
            if not (path/condition/"timeseries.csv").exists(): missing.append(f"{condition}/timeseries.csv")
        status="success" if not missing else "incomplete"
        if missing: warnings.warn(f"seed {seed}: missing {', '.join(missing)}")
        runs.append(SeedRun(seed,path,status,tuple(missing)))
    return runs

def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))

def snapshot_path(run: SeedRun, condition: str, phase: int | None=None) -> Path | None:
    paths=sorted((run.path/condition).glob("spatial_phase_*.npz"))
    if not paths: return None
    if phase is None: return paths[-1]
    exact=run.path/condition/f"spatial_phase_{phase:04d}.npz"
    if exact.exists(): return exact
    return min(paths,key=lambda p:abs(int(p.stem.rsplit('_',1)[1])-phase))

def load_snapshot(path: Path) -> dict[str,np.ndarray]:
    with np.load(path) as raw: return {k:raw[k] for k in raw.files}
