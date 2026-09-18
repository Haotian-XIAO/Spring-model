#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,subprocess,sys,warnings
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analysis.aggregate import build,representative
from analysis.collaborator_package import copy_ensemble_outputs,convergence,package_size,render_seed,representative_morphology,write_readmes,find_ffmpeg
from analysis.load_results import discover_seeds
from analysis.standard_reference import reference_areas

def main():
    p=argparse.ArgumentParser(description="Build the FileSender-ready collaborator package from saved paired outputs only.")
    p.add_argument('--input-root',type=Path,required=True);p.add_argument('--analysis-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--representative-only',action='store_true');p.add_argument('--seeds',type=int,nargs='*');p.add_argument('--skip-videos',action='store_true');p.add_argument('--finalize-n100',action='store_true')
    a=p.parse_args();runs=discover_seeds(a.input_root);complete=[r for r in runs if r.status=='success'];incomplete=[r for r in runs if r.status!='success']
    for run in incomplete:warnings.warn(f"ignoring incomplete seed {run.seed}: {', '.join(run.problems)}")
    if a.finalize_n100:
        expected=set(range(1001,1101));actual={r.seed for r in complete}
        if actual!=expected: p.error(f"N=100 finalization requires exactly complete seeds 1001-1100; missing={sorted(expected-actual)}, unexpected={sorted(actual-expected)}")
        ffmpeg,_=find_ffmpeg()
        if not ffmpeg:p.error('N=100 finalization requires ffmpeg for H.264 MP4 generation, but ffmpeg was not found on PATH')
        subprocess.run([sys.executable,str(Path(__file__).with_name('run_analysis.py')),'--input-root',str(a.input_root),'--output-root',str(a.analysis_root)],check=True)
    refs=reference_areas(complete);rows,time=build(complete,reference_by_seed=refs);reps=representative(rows);byseed={r.seed:r for r in complete}
    if a.seeds is not None:selected=a.seeds
    elif a.representative_only:selected=[int(r['seed']) for r in reps]
    else:selected=sorted(byseed)
    unknown=sorted(set(selected)-set(byseed))
    if unknown:p.error(f"requested seeds are not complete paired runs: {unknown}")
    if not a.finalize_n100 and not a.representative_only and a.seeds is None:p.error('refusing bulk rendering without --finalize-n100; use --representative-only or explicit --seeds for testing')
    write_readmes(a.output_root,refs);copy_ensemble_outputs(a.analysis_root,a.output_root)
    seed_root=a.output_root/'02_seed_results';results=[render_seed(byseed[s],seed_root,make_video=not a.skip_videos,reference_area=refs[s]) for s in selected]
    rendered=set(selected);representative_morphology(a.output_root,[r for r in reps if int(r['seed']) in rendered])
    convergence(rows,a.output_root/'03_numerical_summary')
    total,videos=package_size(a.output_root);print(f"complete paired seeds discovered: {len(complete)}; rendered seeds: {len(results)}; package bytes: {total}")
    if videos:print(f"videos: {len(videos)}; min/median/max bytes: {min(videos)}/{sorted(videos)[len(videos)//2]}/{max(videos)}")
    elif a.skip_videos:print('videos skipped by request')
if __name__=='__main__':main()
