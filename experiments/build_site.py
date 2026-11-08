#!/usr/bin/env python3
"""Publish the canonical notebook results and the current report."""
import base64
import hashlib
import csv
import json
import re
from pathlib import Path
from site_layout import reorganize_methods

ROOT=Path(__file__).resolve().parents[1]


def build():
    notebook_path=ROOT/'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb'
    folder=ROOT/'results/notebook'
    manifest=json.loads((folder/'manifest.json').read_text())
    assert manifest['notebook_sha256']==hashlib.sha256(notebook_path.read_bytes()).hexdigest()
    assert manifest['summary_sha256']==hashlib.sha256((folder/'summary.csv').read_bytes()).hexdigest()
    assert manifest['seed']==236 and manifest['evaluation_reps']==600
    with (folder/'summary.csv').open(newline='') as stream:
        data=[]
        for row in csv.DictReader(stream):
            data.append({k:(v if k in ['study','task','method'] else float(v))
                         for k,v in row.items() if v!=''})
    empirical=(ROOT/'site/empirical.html').read_text(encoding='utf-8')
    source=(ROOT/'site/base-v4.html').read_text(encoding='utf-8')
    insertion='<section aria-labelledby="mc-title" class="montecarlo section-pad" id="montecarlo">'
    assert source.count(insertion)==1
    source=source.replace(insertion,empirical+'\n'+insertion)
    source=source.replace('href="#montecarlo">Simulations','href="#empirical">Simulations')
    source=source.replace('href="#rectangles">Geometry</a>','')
    source=source.replace('<h2 id="mc-title">Monte Carlo comparison</h2>','<h2 id="mc-title">Run a simulation</h2>')
    source=source.replace('The current Monte Carlo benchmark compares global sum, exhaustive scan, and column pooling under the same simulated Gaussian matrices. The additional detection and recovery procedures above are not yet included in this benchmark.',
        'Try the notebook’s four detectors on fresh browser-generated matrices. Error is half the false-positive plus false-negative rate; lower is better. This independent browser run uses seed 236, but a different random generator from Python, so it does not reproduce the saved observations above.')
    source=source.replace('The sum uses its exact Gaussian 5% cutoff. Pooling uses its exact null maximum cutoff. For the exhaustive scan, a separate set of 250 null matrices estimates a 5% rejection threshold; fresh null and signal matrices evaluate performance. This is an empirical comparison, not an asymptotic theorem.',
        'Sum and maximum use exact Gaussian 5% cutoffs. Exact scan uses 1,500 independent null matrices for calibration. Sum+scan splits the 5% level between its components. All four detectors share the same evaluation matrices. Recovery is evaluated in the notebook.')
    source=source.replace('Run simulations to compare empirical power.','Run simulations to compare detection error.')
    start=source.index('  function statsForMonteCarlo(sample)')
    end=source.index('  // Debuggable, opt-in public entry point',start)
    source=source[:start]+(ROOT/'site/live-simulation.js').read_text(encoding='utf-8')+'\n'+source[end:]
    # Add the notebook's dense block sizes to the existing browser controls.
    source=source.replace('<option value="20">20 × 20</option>','')
    source=source.replace('<option value="4">4 × 4</option>',
        '<option value="4">4 × 4</option><option value="6">6 × 6</option><option value="8">8 × 8</option>')
    source=source.replace('id="mc-reps" max="300"','id="mc-reps" max="600"')
    start=source.index('<section aria-labelledby="rectangles-title" class="rectangles section-pad" id="rectangles">')
    end=source.index('</section>',start)+len('</section>')
    compact='''<section id="rectangles" class="section-pad geometry-note" aria-labelledby="rectangles-title"><div class="container"><div class="section-head"><h2 id="rectangles-title">Beyond square blocks</h2></div><p>The report and simulations consider <strong>square <em>n × n</em> Gaussian matrices containing square <em>k × k</em> planted submatrices</strong>, with unknown, not necessarily consecutive row and column sets. Rectangular affected blocks can behave differently: pooling across a long dimension may then be helpful. They appear in the literature's broader lower-bound formulations, but are not a separate focus of our present analysis.</p></div></section>'''
    source=source[:start]+compact+source[end:]
    source=reorganize_methods(source)
    report=ROOT/'STAT236_Project_With_Simulations.pdf'
    assert report.is_file(),'Compile the current report before building the website.'
    encoded='data:application/pdf;base64,'+base64.b64encode(report.read_bytes()).decode('ascii')
    source,count=re.subn(r'data:application/pdf;base64,[A-Za-z0-9+/=]+',encoded,source)
    assert count==2,count
    css=(ROOT/'site/empirical.css').read_text(encoding='utf-8')
    js=(ROOT/'site/empirical.js').read_text(encoding='utf-8')
    css+='\n'+(ROOT/'site/procedure-tree.css').read_text(encoding='utf-8')
    js+='\n'+(ROOT/'site/procedure-tree.js').read_text(encoding='utf-8')
    source=source.replace('</head>','<style>\n'+css+'\n</style>\n</head>',1)
    payload='<script id="notebook-provenance-json" type="application/json">'+json.dumps(manifest,separators=(',',':'))+'</script>'
    payload+='<script id="empirical-results-json" type="application/json">'+json.dumps(data,separators=(',',':'),allow_nan=False)+'</script>'
    source=source.replace('</body>',payload+'\n<script>\n'+js+'\n</script>\n</body>',1)
    source='\n'.join(line.rstrip() for line in source.splitlines())+'\n'
    for path in [ROOT/'index.html',ROOT/'site/index.html']:
        path.write_text(source,encoding='utf-8')
    print('Built root and site/index.html from the canonical notebook: canonical numerical results, seed 236.')


if __name__=='__main__':
    build()
