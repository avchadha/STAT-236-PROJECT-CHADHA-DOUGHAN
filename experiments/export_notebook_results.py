"""Export website data by running the canonical notebook code, without changing it."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
path = ROOT/'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb'
notebook = nbformat.read(path, as_version=4)
namespace = {}
plt.show = lambda: plt.close('all')
for i, cell in enumerate(notebook.cells):
    if cell.cell_type == 'code':
        print(f'Running notebook cell {i}', flush=True)
        exec(compile(cell.source, f'notebook-cell-{i}', 'exec'), namespace)

out = ROOT/'results/notebook'
out.mkdir(parents=True, exist_ok=True)
detection = namespace['detection']
records = []
for (n,k), group in detection.groupby(['n','k']):
    fpr = group.loc[group.signal == 0, namespace['detectors']].mean()
    for mu, planted in group[group.signal == 1].groupby('mu'):
        for method in namespace['detectors']:
            fnr = 1-planted[method].mean()
            records.append(dict(study='finite',task='detection',n=n,k=k,mu=mu,
                                method=method,error=(fpr[method]+fnr)/2,
                                false_positive=fpr[method],false_negative=fnr))
for r in namespace['recovery_error_rates'].to_dict('records'):
    for method in ['exhaustive','spectral','soft spectral']:
        if pd.notna(r[method]):
            records.append(dict(study='finite',task='recovery',n=r['n'],k=r['k'],
                                mu=r['mu'],method=method,error=r[method]))
for r in namespace['growing'].to_dict('records'):
    for task, method in [('detection','sum'),('recovery','spectral')]:
        records.append(dict(study='growing',task=task,n=r['n'],k=r['k'],mu=r['mu'],
                            beta=r['beta'],method=method,error=r[task+' error']))
summary = pd.DataFrame(records)
summary.to_csv(out/'summary.csv',index=False)
manifest = dict(notebook='notebooks/'+path.name,
                notebook_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                seed=236,evaluation_reps=600,scan_calibration_reps=1500,
                finite_matrices=len(namespace['data']),
                summary_records=len(summary),
                summary_sha256=hashlib.sha256((out/'summary.csv').read_bytes()).hexdigest())
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Exported {len(summary)} summary points from the canonical notebook.',flush=True)
