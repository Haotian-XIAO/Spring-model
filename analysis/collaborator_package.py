from __future__ import annotations

import csv
import shutil
import subprocess
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

from .aggregate import build, representative, write_csv
from .load_results import SeedRun, load_snapshot, snapshot_path
from .standard_reference import reference_areas

REGION_NAMES = ("Contracted", "Preserved", "Enlarged")
REGION_COLORS = ("#0072B2", "#009E73", "#D55E00")


def classify(normalized_area: np.ndarray, lower: float = .75, upper: float = 1.25) -> np.ndarray:
    values = np.asarray(normalized_area, float)
    return np.where(values < lower, 0, np.where(values <= upper, 1, 2)).astype(np.int8)


def paired_checkpoints(run: SeedRun) -> list[tuple[int, Path, Path]]:
    def indexed(condition: str) -> dict[int, Path]:
        return {int(p.stem.rsplit("_", 1)[1]): p for p in (run.path / condition).glob("spatial_phase_*.npz")}
    control, softening = indexed("control"), indexed("softening")
    missing_control = sorted(set(softening) - set(control)); missing_softening = sorted(set(control) - set(softening))
    if missing_control or missing_softening:
        warnings.warn(f"seed {run.seed}: unmatched checkpoints; control missing {missing_control}, softening missing {missing_softening}")
    return [(phase, control[phase], softening[phase]) for phase in sorted(set(control) & set(softening))]


def geometry_limits(checkpoints: list[tuple[int, Path, Path]]) -> tuple[float, float, float, float]:
    xmin=ymin=float("inf"); xmax=ymax=float("-inf")
    for _, left, right in checkpoints:
        for path in (left, right):
            points=load_snapshot(path)["points"]
            xmin=min(xmin,float(points[:,0].min()));xmax=max(xmax,float(points[:,0].max()))
            ymin=min(ymin,float(points[:,1].min()));ymax=max(ymax,float(points[:,1].max()))
    padx=.025*(xmax-xmin);pady=.025*(ymax-ymin)
    return xmin-padx,xmax+padx,ymin-pady,ymax+pady


def render_region_frame(run: SeedRun, phase: int, control_path: Path, softening_path: Path, output: Path, limits, lower=.75, upper=1.25, reference_area: float | None = None) -> dict[str,float]:
    snapshots=[load_snapshot(control_path),load_snapshot(softening_path)]
    cmap=ListedColormap(REGION_COLORS); norm=BoundaryNorm([-.5,.5,1.5,2.5],cmap.N)
    fig,axes=plt.subplots(1,2,figsize=(9,4.5),dpi=200)
    fractions=[]
    for ax,snapshot,title in zip(axes,snapshots,("No softening","With softening")):
        normalized = snapshot["cell_area"] / reference_area if reference_area is not None else snapshot["normalized_cell_area"]
        classes=classify(normalized,lower,upper)
        fraction=float(np.mean(classes==1));fractions.append(fraction)
        polygons=snapshot["points"][snapshot["hex_cells"]]
        collection=PolyCollection(polygons,array=classes,cmap=cmap,norm=norm,edgecolor="white",linewidth=.08)
        ax.add_collection(collection);ax.set_xlim(limits[0],limits[1]);ax.set_ylim(limits[2],limits[3]);ax.set_aspect("equal");ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(f"{title}\nPreserved fraction: {fraction:.1%}",fontsize=11)
    handles=[plt.Line2D([0],[0],marker='s',linestyle='',markersize=8,color=color,label=name) for name,color in zip(REGION_NAMES,REGION_COLORS)]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.5,.055),ncol=3,frameon=False,fontsize=9)
    fig.suptitle(f"Seed {run.seed} · remodeling phase {phase}",fontsize=13)
    fig.subplots_adjust(left=.02,right=.98,top=.86,bottom=.18,wspace=.03)
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=200,facecolor="white");plt.close(fig)
    return {"phase":phase,"control_preserved_like_fraction":fractions[0],"softening_preserved_like_fraction":fractions[1],"paired_effect":fractions[1]-fractions[0]}


def find_ffmpeg() -> tuple[str|None,str|None]:
    return shutil.which("ffmpeg"), shutil.which("ffprobe")


def encode_video(frames: list[Path], output: Path, ffmpeg: str) -> None:
    manifest=output.parent/".frames.txt"
    lines=[]
    for i,path in enumerate(frames):
        duration=.75 if i in (0,len(frames)-1) else .25
        escaped=str(path.resolve()).replace("'", "'\\''")
        lines.extend((f"file '{escaped}'",f"duration {duration}"))
    lines.append(f"file '{str(frames[-1].resolve()).replace(chr(39), chr(39)+'\\'+chr(39)+chr(39))}'")
    manifest.write_text("\n".join(lines)+"\n",encoding="utf-8")
    command=[ffmpeg,"-hide_banner","-loglevel","error","-y","-f","concat","-safe","0","-i",str(manifest),"-vf","fps=4,format=yuv420p","-c:v","libx264","-crf","23","-an","-movflags","+faststart",str(output)]
    try: subprocess.run(command,check=True)
    finally: manifest.unlink(missing_ok=True)


def render_seed(run: SeedRun, output_root: Path, *, make_video=True, lower=.75, upper=1.25, reference_area: float | None = None) -> dict:
    checkpoints=paired_checkpoints(run)
    if not checkpoints: raise ValueError(f"seed {run.seed}: no paired spatial checkpoints")
    seed_root=output_root/f"seed_{run.seed}"; figures=seed_root/"figures";limits=geometry_limits(checkpoints)
    rows=[];frames=[]
    for phase,left,right in checkpoints:
        frame=figures/f"seed_{run.seed}_phase_{phase:04d}_regions.png"
        rows.append(render_region_frame(run,phase,left,right,frame,limits,lower,upper,reference_area));frames.append(frame)
    write_csv(seed_root/f"seed_{run.seed}_summary.csv",rows)
    video=seed_root/f"seed_{run.seed}_evolution.mp4"
    if make_video:
        ffmpeg,_=find_ffmpeg()
        if not ffmpeg: raise RuntimeError("ffmpeg is required for MP4 generation but was not found on PATH; no software was installed")
        encode_video(frames,video,ffmpeg)
    return {"seed":run.seed,"checkpoints":len(checkpoints),"video":video if video.exists() else None}


def copy_ensemble_outputs(analysis_root: Path, package_root: Path) -> None:
    target=package_root/"01_ensemble_results";target.mkdir(parents=True,exist_ok=True)
    for folder in ("figures_main","figures_supplementary","representative_seeds","report"):
        source=analysis_root/folder
        if source.exists(): shutil.copytree(source,target/folder,dirs_exist_ok=True)
    numerical=package_root/"03_numerical_summary";numerical.mkdir(parents=True,exist_ok=True)
    if (analysis_root/"tables").exists(): shutil.copytree(analysis_root/"tables",numerical/"tables",dirs_exist_ok=True)


def representative_morphology(package_root: Path, representative_rows: list[dict]) -> None:
    if not representative_rows:return
    fig,axes=plt.subplots(len(representative_rows),1,figsize=(7.48,3.55*len(representative_rows)))
    axes=np.atleast_1d(axes)
    for ax,row in zip(axes,representative_rows):
        seed=int(row["seed"]);frames=sorted((package_root/"02_seed_results"/f"seed_{seed}"/"figures").glob("*.png"))
        if not frames:ax.set_axis_off();continue
        ax.imshow(plt.imread(frames[-1]));ax.set_axis_off();ax.set_title(row["category"].replace('_',' ').title(),fontsize=9,pad=2)
    fig.suptitle("Representative final morphology classifications",fontsize=10);fig.tight_layout()
    out=package_root/"01_ensemble_results"/"main_04_representative_morphology"
    for ext in ("pdf","png","svg"):fig.savefig(out.with_suffix('.'+ext),dpi=300,bbox_inches="tight")
    plt.close(fig)


def convergence(seed_rows: list[dict], output_root: Path, checkpoints=(24,48,72,100)) -> list[dict]:
    ordered=sorted((r for r in seed_rows if r.get("run_status")=="success"),key=lambda r:int(r["seed"])); rows=[]
    for count in checkpoints:
        if len(ordered)<count: continue
        subset=ordered[:count]; preserved=np.array([float(r["delta_operational_preserved_cell_fraction"]) for r in subset]); cv=np.array([float(r["delta_normalized_area_cv"]) for r in subset]);tol=.005
        rows.append({"N":count,"preserved_effect_median":np.median(preserved),"preserved_effect_q25":np.percentile(preserved,25),"preserved_effect_q75":np.percentile(preserved,75),"protective_count":int((preserved>tol).sum()),"near_zero_count":int((np.abs(preserved)<=tol).sum()),"aggravating_count":int((preserved<-tol).sum()),"area_cv_effect_median":np.median(cv),"area_cv_effect_q25":np.percentile(cv,25),"area_cv_effect_q75":np.percentile(cv,75)})
    write_csv(output_root/"convergence_summary.csv",rows)
    if rows:
        fig,axes=plt.subplots(1,2,figsize=(7.48,2.5));n=np.array([r["N"] for r in rows])
        for ax,key,title in ((axes[0],"preserved_effect","Operational preserved-like fraction"),(axes[1],"area_cv_effect","Normalized-area CV")):
            med=np.array([r[f"{key}_median"] for r in rows]);lo=np.array([r[f"{key}_q25"] for r in rows]);hi=np.array([r[f"{key}_q75"] for r in rows]);ax.plot(n,med,"o-",color="#6A3D9A");ax.fill_between(n,lo,hi,color="#6A3D9A",alpha=.18);ax.axhline(0,color=".4",lw=.7);ax.set_xticks(n);ax.set_xlabel("Number of paired seeds");ax.set_ylabel("Softening − control");ax.set_title(title)
        fig.suptitle("Ensemble-size convergence of paired morphometric effects");fig.tight_layout()
        for ext in ("pdf","png","svg"):fig.savefig((output_root/"main_05_ensemble_convergence").with_suffix('.'+ext),dpi=400,bbox_inches="tight")
        plt.close(fig)
    return rows


def write_readmes(package_root: Path, reference_by_seed: dict[int, float] | None = None) -> None:
    readme=package_root/"00_README";readme.mkdir(parents=True,exist_ok=True)
    (readme/"README.md").write_text("""# FF-like softening paired pulmonary-fibrosis pilot

Each paired simulation contains an initial state followed by 600 remodeling phases. Spatial images and videos show only the saved spatial checkpoints; no intermediate states are interpolated or invented.

Control and FF-like softening start from the same initial biological realization. The region-identification workflow is automated directly from the saved alveolar-like polygon geometry. Control and softening are displayed side by side at the same phase and spatial scale.

The three region classes are Contracted, Preserved, and Enlarged. These are operational model classifications, not validated histological labels. For this standard-reference reanalysis, each seed uses one fixed scalar `A_ref`, the median geometric area of all polygons in that seed's phase-0 control initial geometry. No robust polygon-level healthy/non-fibrotic mapping was available, so all phase-0 polygons were included. Every class uses `geometric area_i(t) / A_ref`; no polygon-by-polygon self-normalization is used. Paired phase-0 geometries were verified identical.

**Question for Jean-François and Patrice:** Do the automatically identified preserved-like and remodeled regions correspond to the preserved/remodeled alveolar areas previously assessed semi-manually? Which additional shape, continuity, or size criteria would make the automated classification biologically faithful?

No raw NPZ simulation files are included in this delivery package.
""",encoding="utf-8")
    (readme/"region_definition.md").write_text("""# Operational region definition

For polygon *i* at phase *t*, `r_i(t) = geometric area_i(t) / A_ref`, where `A_ref` is one fixed per-seed scalar: the median geometric area of all phase-0 polygons in that seed's control initial geometry.

- Contracted / low-area: `r_i < 0.75`
- Preserved-like: `0.75 <= r_i <= 1.25`
- Enlarged / remodeled: `r_i > 1.25`

The thresholds are an explicit operational model definition for collaborator review. They have not been validated against histology or Patrice Callard's semi-automatic workflow.
""",encoding="utf-8")


def package_size(root: Path) -> tuple[int,list[int]]:
    files=[p for p in root.rglob('*') if p.is_file()];videos=[p.stat().st_size for p in files if p.suffix.lower()=='.mp4']
    return sum(p.stat().st_size for p in files),videos
