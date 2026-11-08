/* Gaussian matrix helpers. These functions never receive the planted location. */
(function(root) {
  'use strict';
  function random(seed) {
    let x = seed >>> 0;
    return function() { x += 0x6D2B79F5; let t = x; t = Math.imul(t ^ t >>> 15, t | 1); t ^= t + Math.imul(t ^ t >>> 7, t | 61); return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  }
  function matrix(n, k, mu, planted, seed) {
    const rng = random(seed), values = new Float64Array(n*n);
    for (let i=0; i<values.length; i+=2) {
      const r = Math.sqrt(-2*Math.log(Math.max(rng(), Number.MIN_VALUE))), angle = 2*Math.PI*rng();
      values[i] = r*Math.cos(angle);
      if (i+1<values.length) values[i+1] = r*Math.sin(angle);
    }
    function permutation() {
      const a = Array.from({length:n}, (_,i)=>i);
      for(let i=n-1; i>0; --i) { const j=Math.floor(rng()*(i+1)); [a[i],a[j]]=[a[j],a[i]]; }
      return a;
    }
    const rows = planted ? permutation().slice(0,k).sort((a,b)=>a-b) : [];
    const cols = planted ? permutation().slice(0,k).sort((a,b)=>a-b) : [];
    for(const i of rows) for(const j of cols) values[i*n+j] += mu;
    return {values, rows, cols};
  }
  function quantile(p) {
    if(!(p>0 && p<1)) throw new RangeError('Probability must be between zero and one.');
    const a=[-39.6968302866538,220.946098424521,-275.928510446969,138.357751867269,-30.6647980661472,2.50662827745924];
    const b=[-54.4760987982241,161.585836858041,-155.698979859887,66.8013118877197,-13.2806815528857];
    const c=[-.00778489400243029,-.322396458041136,-2.40075827716184,-2.54973253934373,4.37466414146497,2.93816398269878];
    const d=[.00778469570904146,.32246712907004,2.445134137143,3.75440866190742];
    function tail(q) { return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
    if(p<.02425) return tail(Math.sqrt(-2*Math.log(p)));
    if(p>1-.02425) return -tail(Math.sqrt(-2*Math.log1p(-p)));
    const q=p-.5, r=q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
  function logChoose(n,k) { k=Math.min(k,n-k); let s=0; for(let j=1;j<=k;j++) s+=Math.log(n-k+j)-Math.log(j); return s; }
  function basic(values,n) { let sum=0,max=-Infinity; for(const x of values){sum+=x; if(x>max)max=x;} return {sum:sum/n,max}; }
  function cutoffs(n,k,alpha=.05) {
    return {sum:quantile(1-alpha), max:quantile(Math.exp(Math.log1p(-alpha)/(n*n))), scan:Math.sqrt(2*(2*logChoose(n,k)-Math.log(alpha)))};
  }
  function canScan(n,k) { const logCount=logChoose(n,k); return logCount<=Math.log(20000)+1e-10 && logCount+Math.log(n*k)<=Math.log(8000000)+1e-10; }
  function scan(values,n,k) {
    if(!canScan(n,k)) return null;
    const sums=new Float64Array(n), chosen=[], order=Array.from({length:n},(_,i)=>i);
    let best=-Infinity, rows=null, cols=null, examined=0;
    function visit(start,depth) {
      if(depth===k) {
        const sorted=order.slice().sort((a,b)=>sums[b]-sums[a]);
        let total=0; for(let j=0;j<k;j++) total+=sums[sorted[j]];
        if(total>best){best=total;rows=chosen.slice();cols=sorted.slice(0,k).sort((a,b)=>a-b);}
        examined++; return;
      }
      for(let i=start;i<=n-(k-depth);i++) {
        chosen.push(i); for(let j=0;j<n;j++)sums[j]+=values[i*n+j];
        visit(i+1,depth+1);
        for(let j=0;j<n;j++)sums[j]-=values[i*n+j]; chosen.pop();
      }
    }
    visit(0,0); return {value:best/k, rows, cols, examined};
  }
  const api={random,matrix,quantile,logChoose,basic,cutoffs,canScan,scan};
  root.MatrixStatistics=api;
  if(typeof module!=='undefined' && module.exports) module.exports=api;
})(typeof self!=='undefined'?self:globalThis);
