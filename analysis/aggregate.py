from __future__ import annotations
import csv, math
import numpy as np
from pathlib import Path
from .load_results import SeedRun, read_csv, snapshot_path, load_snapshot
from .morphometry import metrics
from .standard_reference import reference_areas

def number(x):
    try: return float(x)
    except (TypeError,ValueError): return x

def build(runs:list[SeedRun], lower=.75, upper=1.25, reference_by_seed: dict[int, float] | None = None):
    reference_by_seed = reference_areas(runs) if reference_by_seed is None else reference_by_seed
    seeds=[]; time=[]
    for run in runs:
        row={"seed":run.seed,"run_status":run.status,"failure_or_guard_status":"; ".join(run.problems)}
        if run.status!="success": seeds.append(row); continue
        pair=read_csv(run.path/"pair_timeseries.csv")
        for r in pair:
            out={"seed":run.seed,**{k:number(v) for k,v in r.items()}}
            for condition in ("control", "softening"):
                # The saved time-course already contains mean geometric area;
                # only its denominator changes for this reanalysis.
                mean_area=float(r[f"{condition}_mean_cell_area"])
                out[f"{condition}_mean_normalized_cell_area"]=mean_area/reference_by_seed[run.seed]
            for key in ("mean_normalized_cell_area", "normalized_area_median", "normalized_area_cv", "neighbor_normalized_area_contrast", "operational_preserved_cell_fraction"):
                base=f"{key}"
                if f"control_{base}" in out and f"softening_{base}" in out:
                    out[f"delta_{base}"]=out[f"softening_{base}"]-out[f"control_{base}"]
            time.append(out)
        final=pair[-1]
        row.update({f"endpoint_{k}":number(v) for k,v in final.items() if k!="phase"})
        for k in final:
            if k.startswith("delta_"):
                vals=[abs(float(x[k])) for x in pair if x.get(k)]
                row[f"max_abs_{k}"]=max(vals) if vals else float("nan")
        for condition in ("control","softening"):
            p=snapshot_path(run,condition)
            if p:
                for k,v in metrics(load_snapshot(p),lower,upper,reference_by_seed[run.seed]).items(): row[f"{condition}_{k}"]=v
        for k in list(row):
            if k.startswith("control_") and f"softening_{k[8:]}" in row:
                row[f"delta_{k[8:]}"]=row[f"softening_{k[8:]}"]-row[k]
        seeds.append(row)
    return seeds,time

def representative(rows, metric="delta_operational_preserved_cell_fraction"):
    ok=[r for r in rows if r.get("run_status")=="success" and isinstance(r.get(metric),(int,float)) and math.isfinite(r[metric])]
    if not ok:return []
    ordered=sorted(ok,key=lambda r:(r[metric],r["seed"])); near=sorted(ok,key=lambda r:(abs(r[metric]),r["seed"]))[0]
    med=ordered[(len(ordered)-1)//2]
    picks=[("strongest_aggravating" if ordered[0][metric]<0 else "least_protective",ordered[0]),("near_neutral",near),("median_like",med),("strongest_protective" if ordered[-1][metric]>0 else "least_aggravating",ordered[-1])]
    unique=[];seen=set()
    for item in picks:
        if item[1]["seed"] not in seen: unique.append(item);seen.add(item[1]["seed"])
    picks=unique
    return [{"category":c,"seed":r["seed"],"selection_metric":metric,"paired_effect":r[metric]} for c,r in picks]

def write_csv(path:Path,rows:list[dict]):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    fields=[]
    for r in rows:
        fields.extend(k for k in r if k not in fields)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
