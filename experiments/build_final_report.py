"""Generate final-study figures and the report's simulations section from verified CSVs."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from notebook_studies import read_results, signal_view

ROOT=Path(__file__).resolve().parents[1]
FINAL=ROOT/'results/final'
FIGURES=ROOT/'figures'
LABELS={'sum':'Global sum','scan':'Exact scan','maximum':'Maximum entry','column_pool':'Column pooling',
        'sum_or_scan':'Sum or scan','exhaustive':'Block maximization','spectral':'Spectral',
        'soft_spectral':'Thresholded spectral'}
COLORS={'sum':'#a51c30','scan':'#20232a','maximum':'#7b6baa','column_pool':'#ad8551',
        'sum_or_scan':'#287b8e','exhaustive':'#20232a','spectral':'#a51c30','soft_spectral':'#287b8e'}


def finish(fig,name,methods):
    handles=[Line2D([0],[0],color=COLORS[m],marker='o',markersize=3,label=LABELS[m]) for m in methods]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,-.012),ncol=len(methods),frameon=False)
    fig.tight_layout(rect=(0,.055,1,1))
    for ext in ['pdf','svg','png']:
        path=FIGURES/(name+'.'+ext)
        fig.savefig(path,dpi=180,bbox_inches='tight')
        if ext=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    plt.close(fig)


def main():
    manifest=json.loads((FINAL/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['complete'] and manifest['evaluation_reps']==600
    FIGURES.mkdir(exist_ok=True)
    summary=read_results(FINAL/'thresholds/summary.csv')
    plan=json.loads((FINAL/'thresholds/design.json').read_text())['plan']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
                         'axes.spines.top':False,'axes.spines.right':False})
    for metric,name,methods in [('power','detection_power',['sum','scan','maximum','column_pool','sum_or_scan']),
                                ('exact_recovery','exact_recovery',['exhaustive','spectral','soft_spectral'])]:
        fig,axes=plt.subplots(2,2,figsize=(11.5,7.5))
        for ax,cfg in zip(axes.flat,plan):
            n,k=cfg['n'],cfg['k']
            subset=summary[(summary.n==n)&(summary.k==k)&(summary.metric==metric)]
            for method,group in subset.groupby('method'):
                group=group.sort_values('mu')
                ax.plot(group.mu,group.estimate,'o-',color=COLORS[method],ms=3,lw=1.6)
                ax.fill_between(group.mu.to_numpy(),group.ci_low.to_numpy(),group.ci_high.to_numpy(),color=COLORS[method],alpha=.075)
            ax.set(title=f'n={n}, k={k}; effective gamma={math.log(k)/math.log(n):.3f}',
                   xlabel='Mean shift $\\mu$',ylabel='Power' if metric=='power' else 'Exact recovery probability',
                   ylim=(-.025,1.025))
            ax.grid(axis='y',alpha=.16)
        finish(fig,name,methods)

    beta=read_results(FINAL/'power_laws/summary.csv')
    beta_plan=json.loads((FINAL/'power_laws/design.json').read_text())['plan']
    view=signal_view(beta,beta_plan,'beta')
    fig,axes=plt.subplots(3,4,figsize=(12,8.3),sharey=True)
    for ax,cfg in zip(axes.flat,beta_plan):
        n,k=cfg['n'],cfg['k'];g=cfg['gamma_effective']
        ss=view[(view.n==n)&(view.k==k)&(view.metric=='power')]
        for method in ['sum','scan','maximum','sum_or_scan']:
            group=ss[ss.method==method].sort_values('parameter')
            if group.empty:continue
            ax.plot(group.parameter,group.estimate,'o-',color=COLORS[method],ms=2,lw=1.2)
        ax.axvline(max(g/2,2*g-1),color='#a51c30',ls=':',lw=1)
        ax.axvline(max(0,2*g-1),color='#999',ls='--',lw=1)
        ax.set(title=f'n={n}, k={k}; gamma={g:.3f}',xlabel='$\\beta$',ylim=(-.025,1.025))
        ax.grid(alpha=.12)
    finish(fig,'power_law_regimes',['sum','scan','maximum','sum_or_scan'])

    # Selected points are joined by actual mu, never rounded chart labels.
    selections=[(64,16,'sum',1.,'sum','power'),(64,16,'sum',2.,'sum','power'),
                (64,4,'maximum',1.,'maximum','power'),(12,2,'exact_recovery',1.,'exhaustive','exact_recovery'),
                (12,2,'exact_recovery',2.,'exhaustive','exact_recovery'),
                (64,4,'spectral',1.,'spectral','exact_recovery'),
                (64,4,'spectral',1.,'soft_spectral','exact_recovery'),
                (64,16,'spectral',2.,'spectral','exact_recovery')]
    selected=[]
    for n,k,reference,factor,method,metric in selections:
        cfg=next(c for c in plan if c['n']==n and c['k']==k)
        signal=next(s for s in cfg['signals'] if s['reference']==reference and s['parameter']==factor)
        r=summary[(summary.n==n)&(summary.k==k)&(summary.method==method)&(summary.metric==metric)&
                  np.isclose(summary.mu,signal['mu'],atol=1e-12,rtol=0)].iloc[0]
        selected.append(dict(n=n,k=k,mu=r.mu,method=method,metric=metric,estimate=r.estimate,
                             ci_low=r.ci_low,ci_high=r.ci_high,reference=reference,multiplier=factor))
    import pandas as pd
    pd.DataFrame(selected).to_csv(FINAL/'paper_selected_results.csv',index=False)
    total=manifest['evaluation_trials'];analytic=sum(s['analytic_checks'] for s in manifest['studies'])
    checks=next(s for s in manifest['studies'] if s['name']=='thresholds')['theory_threshold_checks']
    text=r'''% Generated by experiments/build_final_report.py from results/final/.
\section{Simulations}\label{sec:simulations}
We used three fixed simulation designs, summarized in Table~\ref{tab:study-designs}.
Under the null, the entries of each $n\times n$ matrix were independent standard
Gaussians. Under the alternative, we sampled $k$ rows and $k$ columns uniformly
without replacement and added $\mu$ to their intersections. These indices were
unknown to the procedures. All results in this section use $1{,}500$ independent
null calibration matrices per configuration and $600$ fresh null and $600$ fresh
planted detection matrices per signal. Studies with recovery use another $600$
planted matrices per signal. Calibration data are separate from evaluation data;
procedures share the same matrix within each evaluation trial.

\begin{table}[H]\centering\small
\caption{Final simulation designs. Evaluation counts exclude calibration matrices.}
\label{tab:study-designs}
\begin{tabularx}{\textwidth}{@{}lXrr@{}}\toprule
Study & Parameter grid & Pairs & Evaluations \\
\midrule
Broad sweep & Seven matrix sizes; five requested sparsity exponents; five scan-scaled signals & 31 & 279,000 \\
Signal references & Four matrix/block pairs; half, once, and twice each applicable signal reference, plus zero & 4 & 126,000 \\
Power laws & Three matrix sizes; four requested exponents; beta values at and around both detection boundaries & 12 & 172,800 \\
\bottomrule\end{tabularx}
\end{table}

For the broad sweep, $n\in\{6,8,12,16,24,32,64\}$ and requested
$\gamma\in\{.25,.40,.55,.70,.85\}$. We set $k=\lfloor n^\gamma\rfloor$,
clipped to $[2,n-1]$, and $\mu=c\,2\sqrt{\log(n/k)/k}$ for
$c\in\{0,.65,1,1.35,1.8\}$. Repeated integer block sizes are simulated once;
the effective exponent is $\log k/\log n$. The reference study uses
$(n,k)=(12,2),(12,5),(64,4),(64,16)$ and samples both sides of every applicable
sum, scan, maximum-entry, exact-recovery, spectral and thresholded-spectral
signal reference from the theory sections. A reference multiplier of one is
not a finite-sample cutoff: the sum condition and the small-block spectral
condition require diverging ratios, and the large-block spectral guarantee
has an unspecified sufficiently large constant.

In the power-law study, $n\in\{12,24,64\}$, requested
$\gamma\in\{.4,.5,2/3,.8\}$, and $\mu=n^{-\beta}$.
We include $\beta\in\{0,.1,.2,.3,.4,.6,.8\}$, both exponent boundaries
computed from the effective exponent, and offsets $\pm.04$ where nonnegative.
A separate $\mu=0$ evaluation is included; $\beta=0$ instead means $\mu=1$.

\subsection{Calibration and finite-sample error checks}
Detection uses nominal level $.05$ and an upper order-statistic cutoff
$\lceil(m+1)(1-\alpha)\rceil$ from $m=1{,}500$ calibration draws.
Sum-or-scan uses two level-$.025$ cutoffs. Wilson $95\%$ intervals describe
evaluation uncertainty conditional on the cutoff and do not include
calibration variability. Exhaustive scanning and block maximization are
computed only when $\binom nk\leq2{,}500$; uncomputed methods are omitted.
Spectral methods use a full SVD and the known $k$. Small-block recovery also
uses the fixed soft threshold $t=\sqrt{2\log(n/k^2)}$.

'''
    text+=f'All {total:,} evaluation records were checked against their summary counts. '
    text+=f'We compared {analytic:,} sum and maximum-entry rejection rates with their exact Gaussian probabilities. '
    text+=r'''A separate check reused the reference study's held-out matrices with
$c_n=\mu k^2/(2n)$, $h_n=\sqrt{2\log\binom nk^2}$ and
$u_n=\sqrt{4.5\log n}$ from the error calculations. We checked sum errors
against their exact Gaussian probabilities, scan errors against the union
and true-block bounds, and maximum-entry errors against the stated bounds.
'''
    text+=f'All {checks} bound checks were consistent. '
    text+=r'''The exact-probability and error-bound comparisons each use
exact binomial checks with a Bonferroni-adjusted $1\%$ familywise level.
These diagnostics assess finite-sample implementation consistency; they do
not verify asymptotic impossibility or computational-hardness claims.

\subsection{Results across regimes}
Table~\ref{tab:simulation-results} gives selected points from the reference
study. Figures~\ref{fig:power-curves} and~\ref{fig:recovery-curves} show all
sampled signals for its four configurations. The small-matrix cases allow
direct comparison with exact search. The larger cases compare efficient
methods in small- and large-block recovery regimes. Stronger signals yield
high exact-recovery probabilities, while the signal needed for exact recovery
can exceed that needed for detection. The full broad sweep, null checks,
timings and raw trials are available in the executed notebook and repository.

\begin{table}[H]\centering\small
\caption{Selected final-study estimates, with Wilson $95\%$ intervals from $600$ trials.}
\label{tab:simulation-results}
\begin{tabular}{@{}lcll@{}}\toprule
$(n,k)$ & $\mu$ & Procedure / outcome & Estimate [95\% CI] \\
\midrule
'''
    for r in selected:
        outcome='power' if r['metric']=='power' else 'recovery'
        text+=f"$({r['n']},{r['k']})$ & ${r['mu']:.3f}$ & {LABELS[r['method']]} / {outcome} & ${r['estimate']:.3f}\;[{r['ci_low']:.3f},{r['ci_high']:.3f}]$"+r' \\'+'\n'
    text+=r'''\bottomrule\end{tabular}
\end{table}

Figure~\ref{fig:power-law} displays the power-law study. The exponent lines
are asymptotic references, with logarithmic factors omitted; they do not
predict a sharp transition at these finite dimensions. In particular, the
smallest matrices are diagnostics rather than evidence of a limiting theorem.
The original D1--D4 and R1--R4 results have been superseded in this section;
their historical source and data remain available in repository history.

\clearpage
\begin{figure}[p]\centering
\includegraphics[width=\textwidth]{figures/detection_power.pdf}
\caption{Detection power in the signal-reference study. All sampled mean shifts
are shown; shaded bands are marginal Wilson $95\%$ intervals. Exact scanning
and the combined test are present only in feasible small-matrix cases.}
\label{fig:power-curves}\end{figure}
\begin{figure}[p]\centering
\includegraphics[width=\textwidth]{figures/exact_recovery.pdf}
\caption{Exact recovery in the signal-reference study. Both planted index sets
must be correct. The thresholded spectral method is included only for $k^2<n$.}
\label{fig:recovery-curves}\end{figure}
\clearpage
\begin{figure}[p]\centering
\includegraphics[width=\textwidth]{figures/power_law_regimes.pdf}
\caption{Detection across the power-law study. Larger $\beta$ gives weaker
signals. Dotted red and dashed gray lines are the statistical and conditional
computational exponent references for the actual rounded block size.
Exact scanning is omitted when the row-subset cap is exceeded.}
\label{fig:power-law}\end{figure}
\clearpage
'''
    (ROOT/'paper_simulation_section.tex').write_text(text,encoding='utf-8')
    print('Generated three final-study figures and the current simulations section.')


if __name__=='__main__':
    main()
