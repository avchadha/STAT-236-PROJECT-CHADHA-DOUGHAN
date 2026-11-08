#!/usr/bin/env python3
"""Embed the verified final study and current report in the project's static website."""
import base64
import csv
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FINAL=ROOT/'results/final'


def build():
    manifest=json.loads((FINAL/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['complete'] and manifest['calibration_reps']==1500 and manifest['evaluation_reps']==600
    data=[];designs={}
    for study in manifest['studies']:
        name=study['name'];folder=FINAL/name
        design=json.loads((folder/'design.json').read_text(encoding='utf-8'))
        assert design['cal_reps']==1500 and design['eval_reps']==600
        designs[name]=design['plan']
        with (folder/'summary.csv').open(encoding='utf-8',newline='') as stream:
            for r in csv.DictReader(stream):
                assert int(r['replicates'])==600
                data.append(dict(study=name,task=r['task'],s=r['scenario'],n=int(r['n']),k=int(r['k']),
                    mu=float(r['mu']),m=r['method'],metric=r['metric'],y=float(r['estimate']),
                    lo=float(r['ci_low']) if r['ci_low'] else None,hi=float(r['ci_high']) if r['ci_high'] else None,
                    reps=int(r['replicates']),cal=int(float(r['calibration_reps'])) if r['calibration_reps'] else None))
    source=(ROOT/'site/base-v4.html').read_text(encoding='utf-8')
    insertion='<section aria-labelledby="mc-title" class="montecarlo section-pad" id="montecarlo">'
    assert source.count(insertion)==1
    source=source.replace(insertion,(ROOT/'site/empirical.html').read_text(encoding='utf-8')+'\n'+insertion)
    source=source.replace('href="#montecarlo">Simulations','href="#empirical">Simulations')
    source=source.replace('href="#rectangles">Geometry</a>','')
    source=source.replace('<h2 id="mc-title">Monte Carlo comparison</h2>','<h2 id="mc-title">Run a simulation</h2>')
    source=source.replace('The current Monte Carlo benchmark compares global sum, exhaustive scan, and column pooling under the same simulated Gaussian matrices. The additional detection and recovery procedures above are not yet included in this benchmark.',
        'Explore a smaller, independent Monte Carlo study in your browser. The completed 600-replication studies above remain fixed; these live trial results are regenerated when you adjust the controls. The live benchmark evaluates global sum, exhaustive scan, and column pooling.')
    start=source.index('<section aria-labelledby="rectangles-title" class="rectangles section-pad" id="rectangles">')
    end=source.index('</section>',start)+len('</section>')
    compact='''<section id="rectangles" class="section-pad geometry-note" aria-labelledby="rectangles-title"><div class="container"><div class="section-head"><h2 id="rectangles-title">Beyond square blocks</h2></div><p>The report and simulations consider <strong>square <em>n × n</em> Gaussian matrices containing square <em>k × k</em> planted submatrices</strong>, with unknown, not necessarily consecutive row and column sets. Rectangular affected blocks can behave differently: pooling across a long dimension may then be helpful. They appear in the literature's broader lower-bound formulations, but are not a separate focus of our present analysis.</p></div></section>'''
    source=source[:start]+compact+source[end:]
    report=ROOT/'STAT236_Project_With_Simulations.pdf'
    assert report.is_file(),'Compile the current report before building the website.'
    encoded='data:application/pdf;base64,'+base64.b64encode(report.read_bytes()).decode('ascii')
    source,count=re.subn(r'data:application/pdf;base64,[A-Za-z0-9+/=]+',encoded,source)
    assert count==2,count
    css=(ROOT/'site/empirical.css').read_text(encoding='utf-8')
    js=(ROOT/'site/empirical.js').read_text(encoding='utf-8')
    source=source.replace('</head>','<style>\n'+css+'\n</style>\n</head>',1)
    payload='<script id="empirical-results-json" type="application/json">'+json.dumps(data,separators=(',',':'),allow_nan=False)+'</script>'
    payload+='<script id="empirical-metadata-json" type="application/json">'+json.dumps(dict(manifest=manifest,designs=designs),separators=(',',':'),allow_nan=False)+'</script>'
    source=source.replace('</body>',payload+'\n<script>\n'+js+'\n</script>\n</body>',1)
    source='\n'.join(line.rstrip() for line in source.splitlines())+'\n'
    for path in [ROOT/'index.html',ROOT/'site/index.html']:
        path.write_text(source,encoding='utf-8')
    print('Built root and site/index.html:',len(data),'final summary records; all evaluations use 600 replications.')


if __name__=='__main__':
    build()
