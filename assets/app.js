(function(){
  'use strict';
  const S=window.MatrixStatistics, $=id=>document.getElementById(id);
  const canvas=$('matrix'), context=canvas.getContext('2d');
  const state={n:16,k:3,mu:1.73,seed:236,planted:true,reveal:false,group:false};
  let draw=null, worker=null, generation=0, currentStats=null, exact=null, scanPending=false, displayRows=[],displayCols=[];
  const number=x=>Number.isFinite(x)?x.toFixed(3):'—';
  function sync(){
    $('dimension').value=state.n;
    for(const id of ['block-size','block-number']){$(id).max=Math.floor(state.n/2);$(id).value=state.k;}
    $('strength').value=state.mu; $('strength-number').value=state.mu;
    $('planted').checked=state.planted; $('reveal').checked=state.reveal; $('group').checked=state.group;
    $('reveal').disabled=!state.planted; $('group').disabled=!state.planted;
    $('matrix-summary').textContent=state.n+' × '+state.n+' · '+(state.planted?state.k+' × '+state.k+' block':'noise only');
    $('seed-label').textContent='Draw '+state.seed;
    const lower=Math.sqrt(state.n), upper=state.n**(2/3);
    $('regime').textContent=state.k<lower?'k < √n · small block':state.k<upper?'√n ≤ k < n²ᐟ³ · intermediate block':'k ≥ n²ᐟ³ · larger block';
    document.querySelectorAll('[data-gamma]').forEach(button=>button.setAttribute('aria-pressed',state.k===Math.floor(state.n**Number(button.dataset.gamma))?'true':'false'));
    canvas.setAttribute('aria-label',`${state.n} by ${state.n} Gaussian matrix, ${state.planted?'with a hidden '+state.k+' by '+state.k+' block':'noise only'}. The mean shift is ${state.mu}.`);
  }
  function paint(){
    const {n}=state, rowSet=new Set(draw.rows),colSet=new Set(draw.cols);
    const natural=Array.from({length:n},(_,i)=>i);
    displayRows=state.group&&state.planted?[...draw.rows,...natural.filter(i=>!rowSet.has(i))]:natural;
    displayCols=state.group&&state.planted?[...draw.cols,...natural.filter(i=>!colSet.has(i))]:natural;
    const raster=document.createElement('canvas');raster.width=n;raster.height=n;
    const rc=raster.getContext('2d'), pixels=rc.createImageData(n,n);
    const center=[250,248,244],blue=[51,94,133],red=[180,81,49];
    for(let y=0;y<n;y++)for(let x=0;x<n;x++){
      const value=draw.values[displayRows[y]*n+displayCols[x]], mix=Math.min(1,Math.abs(value)/3),end=value<0?blue:red,base=(y*n+x)*4;
      for(let c=0;c<3;c++)pixels.data[base+c]=Math.round(center[c]+mix*(end[c]-center[c]));pixels.data[base+3]=255;
    }
    rc.putImageData(pixels,0,0);
    context.clearRect(0,0,640,640);context.fillStyle='#fffefa';context.fillRect(0,0,640,640);
    context.imageSmoothingEnabled=false;const offset=state.reveal&&state.planted?12:0, size=640-offset;
    context.drawImage(raster,offset,offset,size,size);
    if(n<=32){context.strokeStyle='rgba(255,255,255,.2)';context.lineWidth=.5;for(let i=1;i<n;i++){const p=offset+i*size/n;context.beginPath();context.moveTo(p,offset);context.lineTo(p,640);context.moveTo(offset,p);context.lineTo(640,p);context.stroke();}}
    if(offset){context.fillStyle='#a1482c';displayRows.forEach((r,i)=>{if(rowSet.has(r))context.fillRect(0,offset+i*size/n,8,size/n);});displayCols.forEach((c,j)=>{if(colSet.has(c))context.fillRect(offset+j*size/n,0,size/n,8);});}
    $('matrix-caption').textContent=state.group&&state.planted?'Rows and columns are reordered using the planted membership for this display only. The statistics use the original matrix.':'The same fixed color scale is used for every sample. Hover over a cell to read its value.';
  }
  function statistic(name,value,cutoff){
    $(name+'-value').textContent=number(value);$(name+'-value').classList.remove('unavailable');
    $(name+'-detail').textContent='Cutoff '+number(cutoff)+(name==='scan'?' · conservative bound':'');
    const decision=$(name+'-decision');decision.textContent=value>cutoff?'Above cutoff':'Below cutoff';decision.className=value>cutoff?'above':'';
  }
  function updateCombined(){
    const cut=S.cutoffs(state.n,state.k,.025), out=$('combined-decision');out.className='';
    if(currentStats.sum>cut.sum||(exact&&exact.value>cut.scan)){out.textContent='Above a component cutoff';out.className='above';}
    else if(exact){out.textContent='Below both cutoffs';}
    else out.textContent=scanPending?'Waiting for scan':'Needs the scan';
  }
  function sample(){
    generation++;if(worker){worker.terminate();worker=null;}
    exact=null;scanPending=false;draw=S.matrix(state.n,state.k,state.mu,state.planted,state.seed);
    sync();paint();currentStats=S.basic(draw.values,state.n);
    const cut=S.cutoffs(state.n,state.k);
    statistic('sum',currentStats.sum,cut.sum);statistic('max',currentStats.max,cut.max);
    const possible=S.canScan(state.n,state.k);
    $('scan-value').className='unavailable';$('scan-decision').className='';
    if(!possible){
      $('scan-value').textContent='Not computed';$('scan-decision').textContent='Search too large';
      const logCount=S.logChoose(state.n,state.k),count=logCount<Math.log(1e9)?Math.round(Math.exp(logCount)).toLocaleString():'about 10^'+(logCount/Math.LN10).toFixed(1);
      $('scan-detail').textContent=count+' row sets. No substitute statistic is used.';
      document.querySelector('.statistics').setAttribute('aria-busy','false');updateCombined();return;
    }
    scanPending=true;$('scan-value').textContent='Computing…';$('scan-decision').textContent='Exact row search';$('scan-detail').textContent='The planted location is hidden from this search.';
    document.querySelector('.statistics').setAttribute('aria-busy','true');updateCombined();
    const id=generation;
    function done(result,error){if(id!==generation)return;scanPending=false;document.querySelector('.statistics').setAttribute('aria-busy','false');if(error||!result){$('scan-value').textContent='Unavailable';$('scan-decision').textContent='Search did not finish';$('scan-detail').textContent=error||'Try a smaller configuration.';}else{exact=result;statistic('scan',result.value,cut.scan);$('scan-detail').textContent+=' · '+result.examined.toLocaleString()+' row sets';}updateCombined();}
    try{
      worker=new Worker('assets/scan-worker.js');
      worker.onmessage=event=>{if(event.data.id===id&&id===generation){done(event.data.result,event.data.error);worker?.terminate();worker=null;}};
      worker.onerror=()=>{if(id!==generation)return;worker?.terminate();worker=null;done(null,'The scan worker could not load. Serve the site over HTTP.');};
      worker.postMessage({id,values:draw.values.slice().buffer,n:state.n,k:state.k});
    }catch(error){done(null,error.message);}
  }
  let changeTimer;
  function queue(){clearTimeout(changeTimer);changeTimer=setTimeout(sample,80);}
  function parameter(inputId,numberId,key,integer){
    for(const id of [inputId,numberId])$(id).addEventListener('input',event=>{
      const v=Number(event.target.value);if(event.target.value===''||!Number.isFinite(v))return;
      state[key]=Math.max(Number(event.target.min),Math.min(Number(event.target.max),integer?Math.round(v):v));
      $(inputId).value=state[key];$(numberId).value=state[key];queue();
    });
  }
  parameter('block-size','block-number','k',true);parameter('strength','strength-number','mu',false);
  $('dimension').addEventListener('change',event=>{state.n=Number(event.target.value);state.k=Math.min(state.k,Math.floor(state.n/2));sample();});
  document.querySelectorAll('[data-gamma]').forEach(button=>button.addEventListener('click',()=>{state.k=Math.min(Math.floor(state.n/2),Math.max(1,Math.floor(state.n**Number(button.dataset.gamma))));sample();}));
  $('planted').addEventListener('change',event=>{state.planted=event.target.checked;sample();});
  for(const id of ['reveal','group'])$(id).addEventListener('change',event=>{state[id]=event.target.checked;sync();paint();});
  $('sample').addEventListener('click',()=>{clearTimeout(changeTimer);state.seed++;sample();});
  $('reset').addEventListener('click',()=>{clearTimeout(changeTimer);Object.assign(state,{n:16,k:3,mu:1.73,seed:236,planted:true,reveal:false,group:false});sample();});
  $('controls').addEventListener('submit',event=>event.preventDefault());
  canvas.addEventListener('pointermove',event=>{
    const bounds=canvas.getBoundingClientRect(),scale=640/bounds.width,offset=state.reveal&&state.planted?12:0;
    const x=(event.clientX-bounds.left)*scale,y=(event.clientY-bounds.top)*scale;
    const col=Math.floor((x-offset)/(640-offset)*state.n),row=Math.floor((y-offset)/(640-offset)*state.n);
    if(row<0||col<0||row>=state.n||col>=state.n){$('cell-tooltip').hidden=true;return;}
    const r=displayRows[row],c=displayCols[col],t=$('cell-tooltip');
    t.textContent=`X[${r+1}, ${c+1}] = ${draw.values[r*state.n+c].toFixed(3)}`;t.hidden=false;
    t.style.left=Math.max(4,Math.min(bounds.width-175,event.clientX-bounds.left+12))+'px';
    t.style.top=Math.max(4,event.clientY-bounds.top-36)+'px';
  });
  canvas.addEventListener('pointerleave',()=>{$('cell-tooltip').hidden=true;});
  sample();
})();
