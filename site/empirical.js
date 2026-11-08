(() => {
 'use strict';
 const $=id=>document.getElementById(id);
 if(!$('empirical-results-json'))return;
 const rows=JSON.parse($('empirical-results-json').textContent),meta=JSON.parse($('empirical-metadata-json').textContent);
 const names={sum:'Global sum',scan:'Exhaustive scan',maximum:'Maximum entry',column_pool:'Column pooling',sum_or_scan:'Sum or scan',exhaustive:'Block maximization',spectral:'Spectral projection',soft_spectral:'Thresholded spectral'};
 const colors={sum:'#a51c30',scan:'#26272e',maximum:'#747b86',column_pool:'#bb8790',sum_or_scan:'#7c2332',exhaustive:'#a51c30',spectral:'#34373f',soft_spectral:'#ad7380'};
 const dash={sum:'',scan:'',maximum:'2,6',column_pool:'9,5,2,5',sum_or_scan:'9,5',exhaustive:'',spectral:'9,5',soft_spectral:'9,5,2,5'};
 const order=['sum','scan','maximum','column_pool','sum_or_scan','exhaustive','spectral','soft_spectral'];
 const notes={broad:'31 configurations; five signals around the scan reference at each matrix/block size.',thresholds:'Four configurations; signals below, at, and above every applicable reference scale.',power_laws:'12 configurations; beta values on both sides of the exponent boundaries. Detection only.'};
 let study='broad',task='detection';
 const pct=x=>(100*x).toFixed(1)+'%',number=x=>Number(x.toPrecision(4)).toString();
 const svgNS='http://www.w3.org/2000/svg';
 function el(tag,attributes,parent){const a=document.createElementNS(svgNS,tag);for(const [k,v] of Object.entries(attributes))a.setAttribute(k,v);parent?.append(a);return a;}
 function text(svg,x,y,value,attributes={}){const t=el('text',{x,y,...attributes},svg);t.textContent=value;return t;}
 function option(select,value,label){const o=document.createElement('option');o.value=value;o.textContent=label;select.append(o);}
 function configuration(){return meta.designs[study].find(c=>`${c.n}_${c.k}`===$('emp-scenario').value);}
 function signals(cfg){return (study==='power_laws'?cfg.signals.filter(s=>s.reference==='beta'):study==='broad'?cfg.signals:cfg.mus.map(mu=>({mu,parameter:mu}))).map(s=>({...s,x:study==='thresholds'?s.mu:s.parameter})).sort((a,b)=>a.x-b.x);}
 function refreshStudy(){
   $('emp-study-note').textContent=notes[study];$('emp-task-recovery').disabled=study==='power_laws';
   if(study==='power_laws')task='detection';
   for(const b of document.querySelectorAll('.emp-tab')){const active=b.dataset.task===task;b.classList.toggle('active',active);b.setAttribute('aria-pressed',String(active));}
   const select=$('emp-scenario');select.replaceChildren();let lastN=null,group=null;
   for(const c of meta.designs[study]){if(c.n!==lastN){group=document.createElement('optgroup');group.label=`n = ${c.n}`;select.append(group);lastN=c.n;}option(group,`${c.n}_${c.k}`,`n = ${c.n}, k = ${c.k} · gamma = ${(Math.log(c.k)/Math.log(c.n)).toFixed(3)}`);}
   const preferred={broad:'16_4',thresholds:'64_4',power_laws:'64_27'}[study];if([...select.options].some(o=>o.value===preferred))select.value=preferred;refreshSignals();
 }
 function refreshSignals(){const values=signals(configuration()),select=$('emp-mu');select.replaceChildren();for(const s of values)option(select,String(s.mu),study==='broad'?`c = ${number(s.x)} · μ = ${number(s.mu)}`:study==='power_laws'?`β = ${number(s.x)} · μ = ${number(s.mu)}`:`μ = ${number(s.mu)}`);select.value=String(values[Math.floor(values.length/2)].mu);draw();}
 function appendRow(body,values){const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value;tr.append(td);}body.append(tr);}
 function draw(){
   const cfg=configuration(),values=signals(cfg),mu=+$('emp-mu').value,metric=task==='detection'?'power':'exact_recovery';
   const all=rows.filter(r=>r.study===study&&r.task===task&&r.n===cfg.n&&r.k===cfg.k&&r.metric===metric),methods=order.filter(m=>all.some(r=>r.m===m));
   const svg=$('emp-chart');svg.replaceChildren();const compact=window.innerWidth<=540,W=compact?440:820,H=compact?365:412,L=55,R=16,T=18,B=58;
   svg.setAttribute('viewBox',`0 0 ${W} ${H}`);svg.setAttribute('aria-label',`${task==='detection'?'Detection power':'Exact recovery'} for n=${cfg.n}, k=${cfg.k}; 600 replications`);
   const points=new Map(values.map(s=>[s.mu,s.x])),xmin=values[0].x,xmax=values.at(-1).x,px=x=>L+(x-xmin)/(xmax-xmin||1)*(W-L-R),py=y=>T+(1-y)*(H-T-B);
   for(let i=0;i<=4;i++){const y=i/4;el('line',{x1:L,x2:W-R,y1:py(y),y2:py(y),class:'emp-grid'},svg);text(svg,L-10,py(y)+4,pct(y),{'text-anchor':'end'});}
   let ticks=values.map(v=>v.x);if(ticks.length>6)ticks=Array.from({length:5},(_,i)=>xmin+(xmax-xmin)*i/4);
   for(const x of ticks)text(svg,px(x),H-B+24,number(x),{'text-anchor':'middle'});
   el('line',{x1:L,x2:L,y1:T,y2:H-B,class:'emp-axis'},svg);el('line',{x1:L,x2:W-R,y1:H-B,y2:H-B,class:'emp-axis'},svg);
   text(svg,L+(W-L-R)/2,H-8,study==='broad'?'Scan-reference multiplier c':study==='power_laws'?'β (larger means weaker signal)':'Mean shift μ',{'text-anchor':'middle'});
   text(svg,16,(H-B+T)/2,task==='detection'?'Power':'Exact recovery',{'text-anchor':'middle',transform:`rotate(-90 16 ${(H-B+T)/2})`});
   const xref=points.get(mu);el('line',{x1:px(xref),x2:px(xref),y1:T,y2:H-B,class:'emp-selected'},svg);
   if(study==='power_laws'){const g=Math.log(cfg.k)/Math.log(cfg.n);for(const [boundary,color,label] of [[Math.max(g/2,2*g-1),'#a51c30','Statistical exponent'],[Math.max(0,2*g-1),'#999','Conditional computational exponent']]){const line=el('line',{x1:px(boundary),x2:px(boundary),y1:T,y2:H-B,stroke:color,'stroke-dasharray':'3 5','stroke-width':1.3},svg);el('title',{},line).textContent=label+': '+number(boundary);}}
   for(const m of methods){
     const group=all.filter(r=>r.m===m&&points.has(r.mu)).sort((a,b)=>points.get(a.mu)-points.get(b.mu)),lo=group.map(r=>`${px(points.get(r.mu))},${py(r.lo)}`),hi=group.map(r=>`${px(points.get(r.mu))},${py(r.hi)}`);
     el('polygon',{points:lo.concat(hi.reverse()).join(' '),fill:colors[m],opacity:.075},svg);
     el('polyline',{points:group.map(r=>`${px(points.get(r.mu))},${py(r.y)}`).join(' '),fill:'none',stroke:colors[m],'stroke-width':2.6,'stroke-dasharray':dash[m]},svg);
     for(const r of group){const c=el('circle',{cx:px(points.get(r.mu)),cy:py(r.y),r:r.mu===mu?5:3.2,fill:'#fff',stroke:colors[m],'stroke-width':1.8},svg);el('title',{},c).textContent=`${names[m]} · μ=${number(r.mu)}: ${pct(r.y)} [${pct(r.lo)}, ${pct(r.hi)}]`;c.style.cursor='pointer';c.addEventListener('click',()=>{$('emp-mu').value=String(r.mu);draw();});}
   }
   const legend=$('emp-legend');legend.replaceChildren();for(const m of methods){const entry=document.createElement('span');entry.className='emp-legend-entry';const key=document.createElement('i');key.className='emp-line-key';key.style.background=colors[m];entry.append(key,document.createTextNode(names[m]));legend.append(entry);}
   $('emp-reading-title').textContent=`At μ = ${number(mu)} · n = ${cfg.n}, k = ${cfg.k}`;
   const head=$('emp-table-heading');head.replaceChildren();for(const label of ['Procedure',task==='detection'?'Power':'Exact recovery','95% CI']){const th=document.createElement('th');th.textContent=label;head.append(th);}
   const body=$('emp-table-body');body.replaceChildren();for(const m of methods){const r=all.find(r=>r.m===m&&r.mu===mu);if(r)appendRow(body,[names[m],pct(r.y),`${pct(r.lo)}–${pct(r.hi)}`]);}
   const nullRows=rows.filter(r=>r.study===study&&r.n===cfg.n&&r.k===cfg.k&&r.metric==='false_positive_rate'&&r.mu===mu);
   $('emp-table-note').textContent=task==='detection'?'Held-out false positives: '+nullRows.map(r=>names[r.m]+' '+pct(r.y)).join(' · ')+'.':'Both index sets must match. Zero means no exact successes in 600 independent trials.';
   const g=Math.log(cfg.k)/Math.log(cfg.n),detect=Math.abs(g-2/3)<1e-10?'detection crossover':g<2/3?'highly sparse detection regime':'moderately sparse detection regime',recover=cfg.k*cfg.k===cfg.n?'recovery crossover':cfg.k*cfg.k<cfg.n?'small-block recovery regime':'large-block recovery regime';
   $('emp-regime').textContent=`Actual gamma = ${g.toFixed(3)} · ${detect} · ${recover}.`;
   let caption=`All ${values.length} sampled signal values are shown. `;if(!methods.includes(task==='detection'?'scan':'exhaustive'))caption+='Exact block search was not run for this configuration; no approximate scan is substituted. ';
   caption+=study==='power_laws'?'Dotted lines mark the asymptotic exponent references, using the actual rounded block size. A separate zero-signal check is saved in the raw results.':'Bands show the evaluation uncertainty. These finite matrices do not establish asymptotic boundaries.';$('emp-caption').textContent=caption;
   const n=cfg.n,k=cfg.k,scan=2*Math.sqrt(Math.log(n/k)/k),scales=[['Sum',n/(k*k)],['Scan',scan],['Maximum entry',2*Math.sqrt(Math.log(n))],['Exact recovery',Math.max((Math.sqrt(2*Math.log(k))+Math.sqrt(2*Math.log(n-k)))/Math.sqrt(k),scan)],['Spectral',Math.sqrt(n)/k+Math.sqrt(Math.log(n)/k)]];
   if(n>k*k)scales.push(['Thresholded spectral',Math.sqrt(Math.log(n/(k*k)))]);const scaleBody=$('emp-scale-body');scaleBody.replaceChildren();for(const [label,a] of scales)appendRow(scaleBody,[label,number(a),number(mu/a)]);
 }
 $('emp-study').addEventListener('change',()=>{study=$('emp-study').value;refreshStudy();});
 for(const b of document.querySelectorAll('.emp-tab'))b.addEventListener('click',()=>{if(b.disabled)return;task=b.dataset.task;for(const x of document.querySelectorAll('.emp-tab')){const active=x===b;x.classList.toggle('active',active);x.setAttribute('aria-pressed',String(active));}draw();});
 $('emp-scenario').addEventListener('change',refreshSignals);$('emp-mu').addEventListener('change',draw);window.addEventListener('resize',draw);
 $('emp-download').addEventListener('click',e=>{e.preventDefault();const header=['study','task','n','k','mu','method','metric','estimate','ci_low','ci_high','replicates','calibration_reps'],lines=[header.join(',')];for(const r of rows.filter(r=>r.study===study))lines.push([r.study,r.task,r.n,r.k,r.mu,r.m,r.metric,r.y,r.lo??'',r.hi??'',r.reps,r.cal??''].join(','));const url=URL.createObjectURL(new Blob([lines.join('\n')+'\n'],{type:'text/csv'})),a=document.createElement('a');a.href=url;a.download=`stat236-${study}-1500-600.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
 const analytic=meta.manifest.studies.reduce((n,s)=>n+s.analytic_checks,0),bounds=meta.manifest.studies.find(s=>s.name==='thresholds').theory_threshold_checks;
 $('emp-validation').textContent=`Verified ${meta.manifest.evaluation_trials.toLocaleString()} trial records and ${analytic.toLocaleString()} sum/maximum rejection rates against exact Gaussian probabilities. The theoretical-threshold study passed ${bounds} checks against the paper's finite-sample error bounds (1% familywise level). Recovery reference curves are reported directly: finite-sample experiments cannot verify asymptotic impossibility or computational hardness.`;
 refreshStudy();
})();
