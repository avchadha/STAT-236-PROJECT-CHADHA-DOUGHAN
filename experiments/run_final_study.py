"""Run and export the single active 1,500/600 study used by notebook, paper, and site."""
import json
import os
from pathlib import Path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import pandas as pd
from threadpoolctl import threadpool_limits

from notebook_studies import (broad_plan, threshold_plan, power_law_plan, run_study,
                              analytic_detection_check)
from verify_final_study import (threshold_consistency, verify_source_thresholds,
                               summarize_recovery,exact_gaussian_consistency)

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT/'results/final'


def main():
    definitions = [('broad',broad_plan()[0],True,(101,202)),
                   ('thresholds',threshold_plan(),True,(301,302)),
                   ('power_laws',power_law_plan(),False,(401,402))]
    source_checks = verify_source_thresholds()
    analytic_family_size=4*sum(len(c['mus']) for _,plan,_,_ in definitions for c in plan)
    manifest = dict(seed=2362026,alpha=.05,calibration_reps=1500,evaluation_reps=600,
                    recovery_reps=600,source_threshold_checks=source_checks,
                    exact_gaussian_family_size=analytic_family_size,studies=[],complete=False)
    for name,plan,recovery,keys in definitions:
        with threadpool_limits(limits=1):
            summary,trials,calibration,folder = run_study(name,plan,ROOT/'outputs_notebook',
                cal_reps=1500,eval_reps=600,recovery=recovery,stream_keys=keys)
        destination = FINAL/name
        destination.mkdir(parents=True,exist_ok=True)
        for filename in ['summary.csv','calibration.csv','configurations.csv','design.json','environment.json']:
            (destination/filename).write_bytes((folder/filename).read_bytes())
        trials.to_csv(destination/'trials.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        analytic=exact_gaussian_consistency(analytic_detection_check(summary,calibration),
                                            family_size=analytic_family_size)
        analytic.to_csv(destination/'analytic_sum_max_check.csv',index=False)
        item=dict(name=name,configurations=len(plan),signal_conditions=sum(len(c['mus']) for c in plan),
                  evaluation_trials=len(trials),summary_rows=len(summary),analytic_checks=len(analytic),
                  recovery=recovery,folder='results/final/'+name)
        if name=='thresholds':
            errors,checks=threshold_consistency(trials)
            errors.to_csv(destination/'theory_threshold_errors.csv',index=False)
            checks.to_csv(destination/'theory_threshold_checks.csv',index=False)
            summarize_recovery(summary,plan).to_csv(destination/'recovery_reference_checks.csv',index=False)
            item['theory_threshold_checks']=len(checks)
            item['theory_threshold_checks_passed']=int(checks.consistent.sum())
        manifest['studies'].append(item)
    manifest['evaluation_trials']=sum(s['evaluation_trials'] for s in manifest['studies'])
    assert manifest['evaluation_trials']==577800
    manifest['calibration_matrices']=1500*sum(s['configurations'] for s in manifest['studies'])
    manifest['complete']=True
    (FINAL/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':
    main()
