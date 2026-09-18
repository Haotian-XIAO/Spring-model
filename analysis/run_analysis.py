#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analysis.load_results import discover_seeds
from analysis.aggregate import build,representative,write_csv
from analysis.plotting import style,endpoint,effect_distributions,timecourses,spatial,mechanism
from analysis.reporting import write_report
from analysis.standard_reference import reference_areas
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);p.add_argument('--preserved-lower',type=float,default=.75);p.add_argument('--preserved-upper',type=float,default=1.25);a=p.parse_args()
 if not 0<a.preserved_lower<1<a.preserved_upper: p.error('preserved bounds must bracket 1 and be positive')
 runs=discover_seeds(a.input_root);refs=reference_areas(runs);rows,time=build(runs,a.preserved_lower,a.preserved_upper,refs);reps=representative(rows)
 if not runs: p.error(f'no paired seed directories discovered under {a.input_root}; symlink targets and pair_metadata.json files should be checked')
 successful=sum(r.status=='success' for r in runs)
 if not successful: p.error(f'discovered {len(runs)} seed directories but none contains a complete paired result')
 if not time: p.error(f'discovered {successful} successful seeds but loaded zero paired time-course rows')
 tables=a.output_root/'tables';write_csv(tables/'master_seed_summary.csv',rows);write_csv(tables/'master_timecourse_summary.csv',time);write_csv(tables/'representative_seed_catalog.csv',reps)
 style();mainfig=a.output_root/'figures_main';endpoint(rows,mainfig);timecourses(time,mainfig);mechanism(time,mainfig);effect_distributions(rows,a.output_root/'figures_supplementary')
 byseed={r.seed:r for r in runs}
 for rep in reps:
  if rep['seed'] in byseed:spatial(byseed[rep['seed']],rep['category'],a.output_root/'representative_seeds',reference_area=refs[rep['seed']])
 write_report(a.output_root/'report'/'analysis_report.md',runs,rows,reps,a.preserved_lower,a.preserved_upper,refs)
 print(f"Analyzed {len(runs)} seeds; outputs: {a.output_root}")
if __name__=='__main__':main()
