'use strict';
const assert = require('node:assert/strict');
const S = require('../assets/statistics.js');
function close(a,b,tolerance=1e-8){assert.ok(Math.abs(a-b)<tolerance,`${a} differs from ${b}`);}
function subsets(n,k){const out=[];function walk(start,a){if(a.length===k){out.push(a.slice());return;}for(let i=start;i<n;i++){a.push(i);walk(i+1,a);a.pop();}}walk(0,[]);return out;}
function brute(values,n,k){let best=-Infinity;for(const rows of subsets(n,k))for(const cols of subsets(n,k)){let total=0;for(const i of rows)for(const j of cols)total+=values[i*n+j];best=Math.max(best,total/k);}return best;}
close(S.quantile(.5),0);close(S.quantile(.975),1.95996398454);close(S.quantile(.05),-1.64485362695);
close(Math.exp(S.logChoose(16,8)),12870,1e-7);
for(const seed of [1,236,1001])for(const n of [3,4,5])for(let k=1;k<n;k++){
  const {values}=S.matrix(n,k,.8,true,seed),result=S.scan(values,n,k);
  close(result.value,brute(values,n,k));assert.equal(result.examined,subsets(n,k).length);
  close(result.rows.reduce((v,i)=>v+result.cols.reduce((r,j)=>r+values[i*n+j],0),0)/k,result.value);
  if(k===1)close(result.value,S.basic(values,n).max);
}
const noise=S.matrix(16,3,0,true,236),shifted=S.matrix(16,3,2,true,236);
assert.deepEqual(noise.rows,shifted.rows);assert.deepEqual(noise.cols,shifted.cols);
const rows=new Set(shifted.rows),cols=new Set(shifted.cols);
for(let i=0;i<16;i++)for(let j=0;j<16;j++)close(shifted.values[i*16+j]-noise.values[i*16+j],rows.has(i)&&cols.has(j)?2:0);
close(S.basic(shifted.values,16).sum-S.basic(noise.values,16).sum,2*9/16);
assert.notDeepEqual(S.matrix(16,3,1,true,236).values,S.matrix(16,3,1,true,237).values);
assert.ok(S.canScan(16,8));assert.ok(!S.canScan(512,12));assert.equal(S.scan(new Float64Array(512*512),512,12),null);
for(const alpha of [.05,.025]){const c=S.cutoffs(16,3,alpha);close(c.scan**2/2,2*S.logChoose(16,3)-Math.log(alpha));assert.ok(c.max>c.sum);}
assert.ok(S.cutoffs(16,3,.025).scan>S.cutoffs(16,3).scan);
console.log('Passed quantiles, exact scan against full enumeration, planted shifts, search limits, and cutoff checks.');

// Exercise the UI controller with a small DOM and worker harness.
// This checks behavior without requiring a browser installation.
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const html=fs.readFileSync(path.join(__dirname,'../index.html'),'utf8');
class Element{
  constructor(attrs=''){this.attrs={};for(const m of attrs.matchAll(/([\w-]+)="([^"]*)"/g))this.attrs[m[1]]=m[2];this.value=this.attrs.value||'';this.min=this.attrs.min||'';this.max=this.attrs.max||'';this.checked=attrs.includes('checked');this.listeners={};this.textContent='';this.classList={remove(){}};this.style={};this.dataset={};}
  addEventListener(name,fn){this.listeners[name]=fn;}
  setAttribute(k,v){this.attrs[k]=v;}
  fire(type){this.listeners[type]?.({target:this,preventDefault(){}});}
}
const elements={};for(const m of html.matchAll(/<[^>]*\bid="([^"]+)"[^>]*>/g))elements[m[1]]=new Element(m[0]);
const drawing={createImageData:(n,m)=>({data:new Uint8ClampedArray(n*m*4)}),putImageData(){},clearRect(){},fillRect(){},drawImage(){},beginPath(){},moveTo(){},lineTo(){},stroke(){}};
elements.matrix.getContext=()=>drawing;
const buttons=[.4,.6,.75].map(g=>{const el=new Element();el.dataset.gamma=String(g);return el;});
const statistics=new Element();
const document={getElementById:id=>elements[id],querySelectorAll:()=>buttons,querySelector:()=>statistics,createElement:()=>({getContext:()=>drawing})};
class WorkerHarness{
 constructor(){this.stopped=false;}
 postMessage(message){const self=this;setTimeout(()=>{if(self.stopped)return;const scope={MatrixStatistics:S,Float64Array,importScripts(){},postMessage:data=>self.onmessage?.({data})};scope.self=scope;vm.createContext(scope);vm.runInContext(fs.readFileSync(path.join(__dirname,'../assets/scan-worker.js'),'utf8'),scope);scope.onmessage({data:message});},0);}
 terminate(){this.stopped=true;}
}
async function controllerChecks(){
 const scope={window:{MatrixStatistics:S},document,Worker:WorkerHarness,Float64Array,setTimeout,clearTimeout,Set,Math,Number};
 vm.createContext(scope);vm.runInContext(fs.readFileSync(path.join(__dirname,'../assets/app.js'),'utf8'),scope);
 const settle=()=>new Promise(resolve=>setTimeout(resolve,160));await settle();
 assert.match(elements['scan-value'].textContent,/^\d/);assert.equal(statistics.attrs['aria-busy'],'false');
 const initial=elements['sum-value'].textContent;
 elements.reveal.checked=true;elements.reveal.fire('change');elements.group.checked=true;elements.group.fire('change');assert.equal(elements['sum-value'].textContent,initial);
 elements.sample.fire('click');await settle();assert.notEqual(elements['sum-value'].textContent,initial);
 elements.reset.fire('click');await settle();assert.equal(elements['sum-value'].textContent,initial);
 elements.dimension.value='512';elements.dimension.fire('change');elements['block-number'].value='12';elements['block-number'].fire('input');elements['strength-number'].value='0';elements['strength-number'].fire('input');await settle();
 assert.equal(elements['scan-value'].textContent,'Not computed');
 assert.equal(elements['combined-decision'].textContent,'Not evaluated');
 elements['strength-number'].value='15';elements['strength-number'].fire('input');await settle();assert.equal(elements['sum-decision'].textContent,'Reject');
 elements.planted.checked=false;elements.planted.fire('change');await settle();assert.equal(elements.reveal.disabled,true);assert.equal(elements.group.disabled,true);
 elements.reset.fire('click');await settle();assert.equal(elements['sum-value'].textContent,initial);
 console.log('Passed controller checks for samples, reset, reveal, grouping, parameter changes, worker results, and skipped scans.');
}
controllerChecks().catch(error=>{console.error(error);process.exitCode=1;});
