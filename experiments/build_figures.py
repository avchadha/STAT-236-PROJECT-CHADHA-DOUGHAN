#!/usr/bin/env python3
"""Build publication-ready figures and a LaTeX-ready numerical results section."""
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np
from run_experiments import ROOT, RESULTS, DETECTION, RECOVERY

FIGURES=ROOT/'figures';FIGURES.mkdir(exist_ok=True)
rows=list(csv.DictReader((RESULTS/'summary.csv').open()))
for r in rows:
  for k in ('mu','estimate','ci_low','ci_high'):
    if r[k]: r[k]=float(r[k])

palette={'sum':'#a51c30','scan':'#26272e','maximum':'#747b86','column_pool':'#bb8790','sum_or_scan':'#7c2332',
         'exhaustive':'#a51c30','spectral':'#34373f','soft_spectral':'#ad7380'}
linestyles={'sum':'-','scan':'-','maximum':':','column_pool':'-.','sum_or_scan':'--',
            'exhaustive':'-','spectral':'--','soft_spectral':'-.'}
labels={'sum':'Global sum','scan':'Exhaustive scan','maximum':'Maximum entry','column_pool':'Column pooling',
        'sum_or_scan':'Sum or scan','exhaustive':'Exhaustive maximizer','spectral':'Spectral projection','soft_spectral':'Thresholded spectral'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.titlesize':12,'axes.labelsize':10,
                     'axes.edgecolor':'#9a9a9a','axes.linewidth':.8,'xtick.color':'#454545','ytick.color':'#454545',
                     'text.color':'#1f2024','axes.labelcolor':'#303035','figure.facecolor':'white','axes.facecolor':'white',
                     'legend.fontsize':8,'savefig.facecolor':'white'})


def select(sc,metric):
    return [r for r in rows if r['scenario']==sc and r['metric']==metric]


def panel(ax,scenario,metric):
    ss=select(scenario['id'],metric)
    names=list(dict.fromkeys(r['method'] for r in ss))
    for name in names:
        rr=sorted((r for r in ss if r['method']==name),key=lambda r:r['mu'])
        xs=np.array([x['mu'] for x in rr]);ys=np.array([x['estimate'] for x in rr]);
        lows=np.array([x['ci_low'] for x in rr]); highs=np.array([x['ci_high'] for x in rr]);
        clr=palette[name]
        ax.plot(xs,ys,color=clr,label=labels[name],linestyle=linestyles[name],linewidth=2,marker='o',markersize=3.5)
        ax.fill_between(xs,lows,highs,color=clr,alpha=.075,linewidth=0)
    ax.set_ylim(-.035,1.035)
    ax.yaxis.set_major_locator(MultipleLocator(.25))
    ax.set_xlabel('Mean shift $\\mu$')
    ax.set_ylabel('Power' if metric=='power' else 'Exact recovery probability')
    ax.grid(axis='y',linewidth=.7,color='#e8e8e8')
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.set_title(f"$n={scenario['n']}$, $k={scenario['k']}$",loc='left',pad=9)
    # A single figure legend below the panel grid avoids obscuring empirical curves.

fig,axs=plt.subplots(2,2,figsize=(12.0,8.0))
for ax,scenario in zip(axs.flat,DETECTION): panel(ax,scenario,'power')
fig.suptitle('Detection power across square Gaussian submatrix models',fontsize=14,weight='medium')
fig.tight_layout(rect=(0,.105,1,.955),w_pad=3,h_pad=3)
keys=['sum','scan','maximum','column_pool','sum_or_scan']
handles=[Line2D([0],[0],color=palette[x],linestyle=linestyles[x],linewidth=2,marker='o',markersize=3.5,label=labels[x]) for x in keys]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=5,frameon=False)
for ext in ('pdf','svg','png'):
    fig.savefig(FIGURES/f'detection_power.{ext}',dpi=210,bbox_inches='tight')
plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(12.0,8.0))
for ax,scenario in zip(axs.flat,RECOVERY): panel(ax,scenario,'exact_recovery')
fig.suptitle('Exact recovery probability across square Gaussian submatrix models',fontsize=14,weight='medium')
fig.tight_layout(rect=(0,.105,1,.955),w_pad=3,h_pad=3)
keys=['exhaustive','spectral','soft_spectral']
handles=[Line2D([0],[0],color=palette[x],linestyle=linestyles[x],linewidth=2,marker='o',markersize=3.5,label=labels[x]) for x in keys]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=3,frameon=False)
for ext in ('pdf','svg','png'):
    fig.savefig(FIGURES/f'exact_recovery.{ext}',dpi=210,bbox_inches='tight')
plt.close(fig)

# Small compact table for the report: selected points, with 95% Wilson intervals.
selections=[('D1',1.2,['sum','scan','maximum']),('D3',.6,['sum','column_pool','maximum']),
            ('D4',2.5,['sum','maximum','column_pool']),('R1',1.75,['exhaustive','spectral']),
            ('R4',3.0,['spectral','soft_spectral'])]
report_lines=[]
for sc,mu,methods in selections:
    for name in methods:
        matches=[r for r in rows if r['scenario']==sc and abs(r['mu']-mu)<1e-9 and r['method']==name and r['metric'] in ('power','exact_recovery')]
        assert len(matches)==1,(sc,mu,name)
        r=matches[0]
        report_lines.append((sc,r['n'],r['k'],mu,labels[name],r['estimate'],r['ci_low'],r['ci_high'],r['metric']))
tex=r'''% Automatically generated from results/summary.csv. Do not hand-edit values.
% IMPORTANT: estimates are finite-sample Monte Carlo results, not claims of asymptotic optimality.
\section{Simulations}\label{sec:simulations}
We simulated independent $n\times n$ matrices with $X_{ij}\sim N(0,1)$ under the null
and added a positive mean shift $\mu$ to a uniformly sampled, arbitrarily located
$k\times k$ block under the alternative. The block was unknown to each procedure.
All tests used a one-sided nominal level $\alpha=0.05$, with $1{,}500$
independent null samples per configuration to calibrate their rejection
thresholds; the combined sum-or-scan procedure used a Bonferroni split of
$\alpha/2$ between its two component thresholds. We evaluated each configuration
on $600$ independently generated null matrices and $600$ independently generated
alternatives at each signal strength. The tests used paired matrices within
each replicate. We recorded empirical Type I error and detection power, and
reported Wilson $95\%$ confidence intervals for each conditional binomial
proportion. Calibration variability is not included in those intervals.

For exact recovery, we evaluated $600$ independent alternative matrices per
configuration and recorded success only when \emph{both} affected index sets
were recovered exactly. The exhaustive scan and maximizing-block estimator
were computed exactly by enumerating row subsets and optimizing the column
subset for each. Consequently, exhaustive procedures were evaluated only
on small matrices. We used the leading singular vectors for spectral projection;
thresholded spectral projection applied fixed soft thresholding before
computing those vectors. The mean shift grid and threshold choices were
fixed before evaluating test outcomes. All random draws used the seeded
NumPy generator; the exact code, per-trial records and figures are archived
with the project repository.

\begin{table}[t]
\centering\small
\caption{Selected finite-sample simulation results. Each entry is an estimated detection power or exact-recovery probability and its Wilson $95\%$ interval from $600$ replications.}
\begin{tabular}{@{}llcll@{}}\toprule
Configuration & $(n,k)$ & $\mu$ & Procedure & Estimate (95\% CI) \\
\midrule
'''
for sc,n,k,mu,method,est,lo,hi,metric in report_lines:
    tex+=f"{sc} & $({n},{k})$ & ${mu:g}$ & {method} & ${est*100:.1f}\\%\;[{lo*100:.1f},\,{hi*100:.1f}]\\%$ \\\n"
tex+=r'''\bottomrule\end{tabular}\label{tab:simulation-results}
\end{table}

Figure~\ref{fig:power-curves} compares detection power, while
Figure~\ref{fig:recovery-curves} shows exact-recovery probabilities.
In the $(n,k)=(16,4)$ setting, scanning improves power over the global sum
at moderate shifts, reflecting its ability to concentrate on a sparse block.
For $(64,16)$, aggregation becomes more effective at weak shifts, whereas
for $(64,4)$ the maximum-entry test can outperform global aggregation when
individual signals are large. Exact recovery can require substantially more
signal than simple detection; in the $(64,4)$ example soft-thresholded
spectral projection improves over the unthresholded spectral baseline at
intermediate shifts. These finite-sample comparisons do not establish
computational lower bounds or prove asymptotic phase transitions.

\begin{figure}[t]\centering
\includegraphics[width=\textwidth]{figures/detection_power.pdf}
\caption{Empirical detection power at nominal level $\alpha=0.05$.
Shaded bands denote marginal Wilson $95\%$ confidence intervals; all
procedures within an experiment are run on paired matrices.}
\label{fig:power-curves}\end{figure}
\begin{figure}[t]\centering
\includegraphics[width=\textwidth]{figures/exact_recovery.pdf}
\caption{Exact recovery probability. An estimate counts as correct only if
both affected row and column index sets match their planted values.}
\label{fig:recovery-curves}\end{figure}
'''
tex='\n'.join(line + ('\\' if line.startswith(('D1 &','D3 &','D4 &','R1 &','R4 &')) else '') for line in tex.split('\n'))
(ROOT/'paper_simulation_section.tex').write_text(tex)
print('Wrote paper section, detection/recovery plots, and table entries:',len(report_lines))