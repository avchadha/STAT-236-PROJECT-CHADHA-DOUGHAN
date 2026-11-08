  // Independent browser experiment; the saved Python results appear above.

  const MONTE_CARLO_TESTS=[{id:'sum',label:'Sum'},{id:'scan',label:'Scan'},

    {id:'max',label:'Maximum entry'},{id:'combined',label:'Sum+scan'}];

  function statsForMonteCarlo({n,k,data}) {

    let total=0,maximum=-Infinity;

    for(const x of data){total+=x;maximum=Math.max(maximum,x);}

    return {sum:total/n,scan:exactScan(data,n,k),max:maximum};

  }

  function updateSimControls(){

    const n=+$('mc-n').value,k=+$('mc-k').value,mu=+$('mc-mu').value,reps=+$('mc-reps').value;

    setText('mc-mu-value',mu.toFixed(2));setText('mc-reps-value',reps);

    setText('mc-candidates',binomial(n,k).toLocaleString()+' row subsets; exact top-k column selection');

    return {n,k,mu,reps};

  }

  let running=false,lastRun=null;

  const pause=()=>new Promise(resolve=>setTimeout(resolve,0));

  function mcProgress(message,fraction){setText('mc-progress-label',message);$('mc-progress').style.width=(100*fraction)+'%';}

  function decisions(s,c){return [s.sum>c.sum,s.scan>c.scan,s.max>c.max,

    s.sum>c.sumHalf||s.scan>c.scanHalf];}

  async function runMonteCarlo(){

    if(running)return;

    const {n,k,mu,reps}=updateSimControls();

    if(binomial(n,k)>20000)throw Error('Choose a smaller exact search.');

    running=true;$('mc-run').disabled=true;$('mc-export').disabled=true;

    const rng=mulberry32(236),calibration=1500,nullScan=[];

    for(let i=0;i<calibration;i++){

      nullScan.push(exactScan(makeMatrix(n,k,0,false,rng).data,n,k));

      if(i%20===0){mcProgress('Calibrating the exact scan…',i/(calibration+2*reps));await pause();}

    }

    nullScan.sort((a,b)=>a-b);

    const quantile=alpha=>nullScan[Math.ceil((calibration+1)*(1-alpha))-1];

    const c={sum:inversePhi(.95),max:inversePhi(Math.pow(.95,1/(n*n))),

      scan:quantile(.05),sumHalf:inversePhi(.975),scanHalf:quantile(.025)};

    const counts=MONTE_CARLO_TESTS.map(()=>[0,0]),rows=[];

    for(const signal of [0,1])for(let i=0;i<reps;i++){

      const s=statsForMonteCarlo(makeMatrix(n,k,mu,Boolean(signal),rng)),pred=decisions(s,c);

      pred.forEach((value,j)=>counts[j][signal]+=Number(value));

      rows.push([signal,signal?mu:0,i+1,s.sum,s.scan,s.max,...pred.map(Number)]);

      if(i%20===0){mcProgress('Evaluating fresh labeled matrices…',(calibration+signal*reps+i)/(calibration+2*reps));await pause();}

    }

    const host=$('mc-results');host.replaceChildren();

    MONTE_CARLO_TESTS.forEach((test,j)=>{

      const fpr=counts[j][0]/reps,fnr=1-counts[j][1]/reps,error=(fpr+fnr)/2;

      const row=document.createElement('div');row.className='mc-result-row';

      const top=document.createElement('div');top.className='mc-result-top';

      const label=document.createElement('strong');label.textContent=test.label;

      const value=document.createElement('strong');value.textContent=(100*error).toFixed(1)+'% error';top.append(label,value);

      const bar=document.createElement('div');bar.className='mc-bar';

      const fill=document.createElement('div');fill.style.width=(100*error)+'%';bar.append(fill);

      const info=document.createElement('p');info.className='mc-result-info';

      info.textContent=`False positive ${(100*fpr).toFixed(1)}% · false negative ${(100*fnr).toFixed(1)}% · ${reps} trials per hypothesis`;

      row.append(top,bar,info);host.append(row);

    });lastRun={n,k,mu,reps,rows};

    mcProgress(`Complete · n=${n}, k=${k}, μ=${mu}: 1,500 calibration + ${reps} null + ${reps} planted matrices.`,1);

    $('mc-run').disabled=false;$('mc-export').disabled=false;running=false;

  }

  function exportMonteCarlo(){

    if(!lastRun)return;

    const header='n,k,seed,signal,mu,rep,sum_stat,scan_stat,max_stat,sum_reject,scan_reject,max_reject,combined_reject';

    const blob=new Blob([header+'\n'+lastRun.rows.map(r=>[lastRun.n,lastRun.k,236,...r].join(',')).join('\n')+'\n'],{type:'text/csv'}),url=URL.createObjectURL(blob),a=document.createElement('a');

    a.href=url;a.download=`stat236-browser-n${lastRun.n}-k${lastRun.k}-seed236.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);

  }

  if($('mc-run')){

    ['mc-n','mc-k','mc-mu','mc-reps'].forEach(id=>$(id).addEventListener('input',updateSimControls));

    $('mc-run').addEventListener('click',()=>runMonteCarlo().catch(err=>{mcProgress(err.message,0);running=false;$('mc-run').disabled=false;}));

    $('mc-export').addEventListener('click',exportMonteCarlo);updateSimControls();

  }
