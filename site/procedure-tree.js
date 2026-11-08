(() => {
 'use strict';
 const tree=document.querySelector('.procedure-tree');if(!tree)return;
 const reduce=window.matchMedia('(prefers-reduced-motion: reduce)');
 const nodes=[...tree.querySelectorAll('.procedure-node')],lanes=[...tree.querySelectorAll('.procedure-lane')];
 if(reduce.matches||!('IntersectionObserver' in window))return;
 const observer=new IntersectionObserver(entries=>{for(const e of entries)if(e.isIntersecting){e.target.classList.add('is-seen');observer.unobserve(e.target);}},{rootMargin:'0px 0px -6% 0px',threshold:.03});
 tree.classList.add('is-animated');for(const node of nodes)observer.observe(node);
 let queued=false;
 function update(){queued=false;if(reduce.matches)return;const height=innerHeight;
   for(const lane of lanes){const box=lane.getBoundingClientRect();lane.style.setProperty('--lane-progress',String(Math.max(0,Math.min(1,(height*.75-box.top)/Math.max(1,box.height)))));}
   for(const node of nodes){const box=node.getBoundingClientRect();if(box.top>height+100||box.bottom<-100)continue;const shift=Math.max(-8,Math.min(8,(height/2-box.top-box.height/2)*.018));node.style.setProperty('--depth-y',shift.toFixed(2)+'px');}
 }
 const queue=()=>{if(!queued){queued=true;requestAnimationFrame(update);}};
 addEventListener('scroll',queue,{passive:true});addEventListener('resize',queue,{passive:true});
 tree.addEventListener('focusin',e=>e.target.closest('.procedure-node')?.classList.add('is-seen'));
 reduce.addEventListener?.('change',()=>{if(reduce.matches){tree.classList.remove('is-animated');observer.disconnect();}queue();});
 update();
})();
