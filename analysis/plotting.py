from __future__ import annotations
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from .load_results import SeedRun, snapshot_path, load_snapshot

CONTROL="#4C78A8"; SOFT="#D55E00"; DELTA="#6A3D9A"
LABELS={"mean_normalized_cell_area":"Relative alveolar area","median_normalized_cell_area":"Median normalized alveolar area","mean_stiffness":"Mean spring stiffness","mean_A":"Remodeling level","mean_stress":"Average spring stress","mean_strain":"Average spring strain","mean_activation":"Biological activation","mean_D":"Agent density"}
def style():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.titlesize":8,"axes.labelsize":7,"legend.fontsize":6.5,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"savefig.dpi":400})
def save(fig,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","png","svg"): fig.savefig(path.with_suffix('.'+ext),bbox_inches="tight")
    plt.close(fig)
def endpoint(rows,out):
    metrics=[("delta_operational_preserved_cell_fraction","Preserved alveoli"),("delta_normalized_area_cv","Variation in alveolar area"),("endpoint_delta_mean_stress","Average spring stress")]
    fig,axs=plt.subplots(1,3,figsize=(7.48,2.35))
    for ax,(key,label) in zip(axs,metrics):
        vals=np.sort(np.array([r[key] for r in rows if isinstance(r.get(key),(int,float))])); x=np.arange(len(vals))
        ax.scatter(x,vals,c=np.where(vals>=0,SOFT,CONTROL),s=14);ax.axhline(0,color=".3",lw=.8);ax.set_title(label);ax.set_xlabel("Simulations ordered by effect");ax.set_ylabel("Difference (with − no softening)")
    fig.suptitle("Final effects of softening",y=1.02);fig.tight_layout();save(fig,out/"main_01_paired_endpoints")
def effect_distributions(rows,out):
    specs=[("delta_operational_preserved_cell_fraction","Operationally preserved alveoli"),("delta_normalized_area_cv","Alveolar-area heterogeneity"),("endpoint_delta_mean_stress","Mean spring stress")]
    fig,axs=plt.subplots(1,3,figsize=(7.48,2.25))
    for ax,(key,label) in zip(axs,specs):
        v=np.sort(np.array([r[key] for r in rows if isinstance(r.get(key),(int,float))])); y=np.arange(1,len(v)+1)/len(v)
        ax.step(v,y,where='post',color=DELTA);ax.axvline(0,color='.35',lw=.8);ax.set_title(label);ax.set_xlabel("Softening − control");ax.set_ylabel("Cumulative seed fraction")
    fig.suptitle("Distribution of paired endpoint effects");fig.tight_layout();save(fig,out/"supp_01_paired_effect_ecdf")
def timecourses(time,out):
    if not time:
        raise ValueError("Cannot plot ensemble time courses: no successful paired time-course rows were loaded. Check input discovery, symlink targets, and pair_timeseries.csv files.")
    fig,axs=plt.subplots(2,3,figsize=(7.48,4.6)); metrics=["mean_normalized_cell_area","mean_stress","mean_strain","mean_A","mean_activation","mean_D"]
    phases=sorted({int(r['phase']) for r in time})
    for ax,m in zip(axs.flat,metrics):
        for prefix,color,label in (("control",CONTROL,"No softening"),("softening",SOFT,"With softening")):
            values=[[float(r[f'{prefix}_{m}']) for r in time if int(r['phase'])==p] for p in phases]
            if not values or any(not row for row in values):
                raise ValueError(f"Cannot plot {m}: one or more phases have no {prefix} observations.")
            arr=np.array(values); med=np.median(arr,axis=1);lo=np.percentile(arr,25,axis=1);hi=np.percentile(arr,75,axis=1)
            ax.plot(phases,med,color=color,label=label);ax.fill_between(phases,lo,hi,color=color,alpha=.18,lw=0)
        ax.set_title(LABELS[m]);ax.set_xlabel("Iteration")
    axs[0,0].legend(frameon=False);fig.suptitle("Evolution across simulations");fig.tight_layout();save(fig,out/"main_02_ensemble_timecourses")
def spatial(run:SeedRun,category:str,out:Path,phases=(0,100,300,600), reference_area: float | None = None):
    fig,axs=plt.subplots(2,len(phases),figsize=(7.48,3.8),sharex=True,sharey=True)
    for col,phase in enumerate(phases):
      data=[]
      for row,condition in enumerate(("control","softening")):
        p=snapshot_path(run,condition,phase)
        if not p: axs[row,col].set_axis_off();continue
        s=load_snapshot(p); polys=s['points'][s['hex_cells']]; vals=(s['cell_area']/reference_area if reference_area is not None else s['normalized_cell_area']); data.append((polys,vals))
      if not data:continue
      vmin=min(np.percentile(v[1],2) for v in data);vmax=max(np.percentile(v[1],98) for v in data)
      for row,(polys,vals) in enumerate(data):
        pc=PolyCollection(polys,array=vals,cmap="viridis",clim=(vmin,vmax),edgecolor="white",linewidth=.08);axs[row,col].add_collection(pc);axs[row,col].autoscale();axs[row,col].set_aspect('equal');axs[row,col].set_title(f"Phase {phase}");
        if col==0:axs[row,col].set_ylabel("Control" if row==0 else "FF-like softening")
    fig.suptitle(f"{category.replace('_',' ').title()} seed {run.seed}: normalized alveolar area");fig.tight_layout();save(fig,out/f"seed_{run.seed}_{category}_spatial")
def mechanism(time,out):
    final={}
    for r in time:
      if int(r['phase'])>=final.get(r['seed'],(-1,None))[0]:final[r['seed']]=(int(r['phase']),r)
    xs=[]
    for _,r in final.values(): xs.append([float(r[f'delta_{m}']) for m in ('mean_E','mean_stiffness','mean_stress','mean_activation','mean_D','mean_normalized_cell_area')])
    a=np.array(xs);fig,axs=plt.subplots(1,5,figsize=(7.48,1.9)); names=["E","stiffness","stress","activation","agent density","alveolar area"]
    for i,ax in enumerate(axs):ax.scatter(a[:,i],a[:,i+1],s=13,c=DELTA,alpha=.75);ax.axhline(0,color='.7',lw=.5);ax.axvline(0,color='.7',lw=.5);ax.set_xlabel(f"Δ {names[i]}");ax.set_ylabel(f"Δ {names[i+1]}")
    fig.suptitle("Paired endpoint links in the proposed mechanical–biological sequence");fig.tight_layout();save(fig,out/"main_03_mechanistic_chain")
