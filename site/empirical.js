(() => {
 'use strict';
 const $=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
 const rows=JSON.parse($('empirical-results-json').textContent);
 const names={sum:'Global sum',scan:'Exhaustive scan',max:'Maximum entry','sum+scan':'Sum+scan',exhaustive:'Block maximization',spectral:'Spectral projection','soft spectral':'Thresholded spectral'};
 const palette=['#a51c30','#26272e','#747b86','#7c2332'],dashes=['','','2,6','9,5'];
 const configs=[[12,2,.4],[12,6,.75],[16,3,.4],[16,8,.75]];
 let study='finite',task='detection';
 const pct=x=>(100*x).toFixed(1)+'%',num=x=>Number(x.toPrecision(4)).toString();
 function element(tag,attrs,parent){const e=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);parent.append(e);return e;}
 function text(svg,x,y,value,attrs={}){const e=element('text',{x,y,...attrs},svg);e.textContent=value;}
 function option(select,value,label){const o=document.createElement('option');o.value=value;o.textContent=label;select.append(o);}
 function tableRow(body,values){const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value;tr.append(td);}body.append(tr);}
 function refresh(){
   $('emp-study-note').textContent=study==='finite'?'Four configurations, eight positive signal strengths; 600 evaluations per condition.':'k = floor(n^.75), μ = n^(−β); β = 0, .25, .75. Sum detection and spectral recovery only.';
   const cfg=$('emp-scenario');cfg.parentElement.hidden=study==='growing';
   if(!cfg.options.length)configs.forEach(([n,k,g])=>option(cfg,`${n}_${k}`,`n = ${n}, k = ${k} · γ = ${g}`));
   const select=$('emp-mu');select.replaceChildren();
   const values=study==='finite'?[.25,.5,.75,1,1.5,2,3,4]:[16,32,64,128,256];
   values.forEach(v=>option(select,v,(study==='finite'?'μ = ':'n = ')+v));select.value=values[3];draw();
 }
 function draw(){
   const [n,k]=study==='finite'?$('emp-scenario').value.split('_').map(Number):[+$('emp-mu').value,Math.floor((+$('emp-mu').value)**.75)];
   const selected=+$('emp-mu').value;
   const all=rows.filter(r=>r.study===study&&r.task===task&&(study==='growing'||r.n===n&&r.k===k));
   const keys=study==='growing'?[0,.25,.75]:[...new Set(all.map(r=>r.method))];
   const label=key=>study==='growing'?`β = ${key}`:names[key];
   const points=key=>all.filter(r=>study==='growing'?r.beta===key:r.method===key).sort((a,b)=>xvalue(a)-xvalue(b));
   const xvalue=r=>study==='growing'?r.n:r.mu,values=[...new Set(all.map(xvalue))].sort((a,b)=>a-b);
   const svg=$('emp-chart');svg.replaceChildren();const compact=innerWidth<=540,W=compact?440:820,H=compact?365:412,L=55,R=16,T=18,B=58,ymax=task==='detection'?.55:1;
   svg.setAttribute('viewBox',`0 0 ${W} ${H}`);svg.setAttribute('aria-label',`${task} error, ${study==='finite'?`n=${n}, k=${k}`:'growing matrix size'}, canonical seed-236 notebook`);
   const px=x=>L+(x-values[0])/(values.at(-1)-values[0])*(W-L-R),py=y=>T+(1-y/ymax)*(H-T-B);
   for(let i=0;i<=4;i++){const y=ymax*i/4;element('line',{x1:L,x2:W-R,y1:py(y),y2:py(y),class:'emp-grid'},svg);text(svg,L-10,py(y)+4,pct(y),{'text-anchor':'end'});}
   for(const x of values)text(svg,px(x),H-B+24,num(x),{'text-anchor':'middle'});
   element('line',{x1:L,x2:L,y1:T,y2:H-B,class:'emp-axis'},svg);element('line',{x1:L,x2:W-R,y1:H-B,y2:H-B,class:'emp-axis'},svg);
   text(svg,(L+W-R)/2,H-8,study==='finite'?'Signal strength μ':'Matrix size n',{'text-anchor':'middle'});
   text(svg,16,(H-B+T)/2,task==='detection'?'Detection error':'Exact recovery error',{'text-anchor':'middle',transform:`rotate(-90 16 ${(H-B+T)/2})`});
   element('line',{x1:px(selected),x2:px(selected),y1:T,y2:H-B,class:'emp-selected'},svg);
   const legend=$('emp-legend');legend.replaceChildren();
   keys.forEach((key,i)=>{const group=points(key),color=palette[i];
     element('polyline',{points:group.map(r=>`${px(xvalue(r))},${py(r.error)}`).join(' '),fill:'none',stroke:color,'stroke-width':2.6,'stroke-dasharray':dashes[i]},svg);
     for(const r of group){const dot=element('circle',{cx:px(xvalue(r)),cy:py(r.error),r:xvalue(r)===selected?5:3.2,fill:'#fff',stroke:color,'stroke-width':1.8},svg);element('title',{},dot).textContent=`${label(key)}: ${pct(r.error)} error`;dot.style.cursor='pointer';dot.addEventListener('click',()=>{$('emp-mu').value=xvalue(r);draw();});}
     const entry=document.createElement('span');entry.className='emp-legend-entry';const line=document.createElement('i');line.className='emp-line-key';line.style.background=color;entry.append(line,document.createTextNode(label(key)));legend.append(entry);
   });
   $('emp-reading-title').textContent=study==='finite'?`At μ = ${selected} · n = ${n}, k = ${k}`:`At n = ${n} · k = ${k}`;
   const body=$('emp-table-body');body.replaceChildren();
   keys.forEach(key=>{const r=points(key).find(r=>xvalue(r)===selected);tableRow(body,[label(key),pct(r.error),study==='growing'?`μ = ${num(r.mu)}`:task==='detection'?`FP ${pct(r.false_positive)} · FN ${pct(r.false_negative)}`:'600 planted matrices']);});
   $('emp-table-note').textContent=task==='detection'?'Error = ½(false positive + false negative). Chance gives 50%; a 5% false-positive rate leaves 2.5% error at perfect power.':'Error means at least one index set is wrong. Partial overlap does not count as exact recovery.';
   $('emp-caption').textContent='Every plotted value is exported from the canonical notebook. Click a point to inspect its error in the table. Small differences may reflect sampling variation.';
   $('emp-finding').textContent=study==='growing'?'The sum mean shift scales as n^(.5−β): detection improves for β = 0 and .25, while β = .75 approaches chance. Spectral recovery improves for fixed μ = 1, but fails for both decaying signals at these sizes. These sequences illustrate trends, not asymptotic impossibility.':task==='detection'?'Small blocks favor scan because sum dilutes their signal across the matrix. Larger blocks make sum effective at much smaller μ. Maximum entry needs stronger individual entries; sum+scan tracks the stronger component, with stricter cutoffs. Finite-sample rankings need not coincide with asymptotic boundaries.':'Recovery needs stronger signals than detection. Exhaustive maximization recovers weaker signals than spectral projection in these trials. Larger blocks are easier at the same μ. Soft thresholding modestly improves the small-block spectral results.';
   const scale=$('emp-scale-body');scale.replaceChildren();tableRow(scale,['n, k',`${n}, ${k}`]);
   if(study==='finite'){tableRow(scale,['Requested γ',configs.find(c=>c[0]===n&&c[1]===k)[2]]);tableRow(scale,['Sum mean shift μk²/n',num(selected*k*k/n)]);tableRow(scale,['Scan reference 2√(log(n/k)/k)',num(2*Math.sqrt(Math.log(n/k)/k))]);}
   else {tableRow(scale,['Detection level 1/n',pct(1/n)]);tableRow(scale,['Spectral sufficient scale √n/k + √(log(n)/k)',num(Math.sqrt(n)/k+Math.sqrt(Math.log(n)/k))]);}
   $('emp-regime').textContent=k*k<n?'Small-block recovery regime: k² < n.':'Large-block recovery regime: k² ≥ n.';
 }
 $('emp-study').addEventListener('change',()=>{study=$('emp-study').value;refresh();});
 for(const b of document.querySelectorAll('.emp-tab'))b.addEventListener('click',()=>{task=b.dataset.task;for(const button of document.querySelectorAll('.emp-tab')){const active=button===b;button.classList.toggle('active',active);button.setAttribute('aria-pressed',active);}draw();});
 $('emp-scenario').addEventListener('change',draw);$('emp-mu').addEventListener('change',draw);addEventListener('resize',draw);
 $('emp-download').href='data:text/csv;charset=utf-8,'+encodeURIComponent(['study,task,n,k,mu,beta,method,error,false_positive,false_negative',...rows.map(r=>[r.study,r.task,r.n,r.k,r.mu,r.beta??'',r.method,r.error,r.false_positive??'',r.false_negative??''].join(','))].join('\n'));
 $('emp-download').download='stat236-notebook-seed236.csv';
 $('emp-validation').textContent='Notebook of record: seed 236; 600 evaluations per condition. The older results/final study remains with the report as historical context and does not supply this section.';
 refresh();
})();
