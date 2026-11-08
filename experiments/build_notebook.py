"""Build the self-contained notebook from the reviewed simulation source."""
import ast
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
engine = (ROOT / 'experiments/run_experiments.py').read_text(encoding='utf-8')
tree = ast.parse(engine)


def definitions(*names):
    return '\n\n'.join(ast.get_source_segment(engine, node) for node in tree.body
                       if isinstance(node, ast.FunctionDef) and node.name in names)


cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text, hidden=False, title='Supporting code'):
    cell = nbf.v4.new_code_cell(text.strip())
    if hidden:
        cell.metadata.update(jupyter={'source_hidden': True},
                             colab={'cellView': 'form'}, cellView='form',
                             tags=['hide-input'])
        cell.source = '#@title ' + title + '\n' + cell.source
    cells.append(cell)


md(r'''
# Sparse submatrix detection and recovery
**Amrit Chadha and Faris Doughan · STAT 236 · Fall 2026**

We observe an $n\times n$ matrix of independent unit-variance Gaussian entries.
Under $H_0$ all means are zero. Under $H_1$, a randomly chosen set of $k$ rows
and $k$ columns has mean $\mu>0$; the indices need not be consecutive.
Detection asks whether a block exists. Exact recovery requires both index sets to be correct.

Run the notebook from top to bottom. All three studies below run as written,
and the saved outputs show their results without changing settings.
Setup, file handling, and plotting definitions are collapsed; their outputs remain visible.

| Study | Varies | Purpose |
|---|---|---|
| 1. Broad sweep | Seven $n$ values, five $\gamma$ values, five signals | Compare small and large blocks across matrix sizes |
| 2. Signal thresholds | Four $(n,k)$ pairs; signals below, at, and above each reference scale | Cover every detection and recovery procedure in the writeup |
| 3. Power laws | Three $n$ values, four $\gamma$ values, a $\beta$ grid | Sample both sides of the detection exponent boundaries |

Each detector uses 1,500 independent calibration matrices at level $\alpha=.05$.
Each signal has 600 fresh null and 600 fresh planted detection matrices;
Studies 1–2 add 600 recovery matrices. Methods share the same matrix within a trial.
Wilson 95% intervals describe evaluation uncertainty, conditional on the calibration cutoff.
''')
code('''
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import csv, hashlib, io, itertools, json, math, platform, time, zipfile
from pathlib import Path
from contextlib import redirect_stdout
from urllib.request import urlopen
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
try:
    from threadpoolctl import threadpool_limits
    _blas_limit = threadpool_limits(limits=1)
except ImportError:
    pass

SEED = 2362026
ALPHA = .05
CALIBRATION_REPS = 1500
EVALUATION_REPS = 600
EXACT_ROW_CAP = 2500
OUTPUT = Path('outputs_notebook')
OUTPUT.mkdir(exist_ok=True)
BASE_RAW = 'https://raw.githubusercontent.com/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/f87e43d1be9c002737c4068c1e84de4408026fba/'
ROW_COMBOS = {}
print('Seed:', SEED, '| calibration:', CALIBRATION_REPS, '| evaluations per signal:', EVALUATION_REPS)
''', hidden=True, title='Imports and fixed settings')
md(r'''
## Model and procedures
The generator selects the hidden indices independently of the noise. The procedures receive $X$ and $k$,
never the planted indices. Setting $\mu=0$ gives a no-signal check.
''')
code(definitions('planted_matrix'))
md(r'''
The sum is normalized by $n$ and the scan by $k$, giving unit null variance for a
fixed sum. Exact scanning enumerates every row set and takes its $k$ best column sums;
this is the same maximization as enumerating both sets. The scan and combined test
run only when $\binom nk\leq2500$. Column pooling is an additional comparator.
''')
code(definitions('exact_scan'), hidden=True, title='Exact scan implementation')
code(definitions('test_statistics', 'detection_decisions'))
md(r'''
Spectral recovery projects onto the leading singular vectors and takes the top $k$
row and column scores. For $k<\sqrt n$, we also soft-threshold entries at the fixed
$t=\sqrt{2\log(n/k^2)}$ before projection. We use a full SVD and known $k$;
the writeup's faster power-iteration cost does not describe this implementation.
''')
code(definitions('recovery_spectral'))
md(r'''
The supporting functions below calibrate each detector once per $(n,k)$,
then test every signal on fresh matrices. The sum-or-scan test uses two
level-.025 cutoffs. Recovery records exact success and average row/column Jaccard overlap.
''')
support = (ROOT / 'experiments/notebook_studies.py').read_text(encoding='utf-8')
support = support.replace('from run_experiments import run_detection, run_recovery\n', '')
code(definitions('wilson', 'empirical_cutoff', 'jaccard', 'run_detection', 'run_recovery')
     + '\n\n' + support, hidden=True, title='Calibration, evaluation, and checkpoint functions')
code('''
# These checks use independent brute-force enumeration and a known strong signal.
check_rng = np.random.default_rng(42)
for n, k in [(5, 2), (6, 3)]:
    for _ in range(6):
        X = check_rng.standard_normal((n, n))
        value, rows, cols = exact_scan(X, k, return_indices=True)
        brute = max(X[np.ix_(a, b)].sum()/k
                    for a in itertools.combinations(range(n), k)
                    for b in itertools.combinations(range(n), k))
        assert abs(value-brute) < 1e-12
        assert abs(X[np.ix_(rows, cols)].sum()/k-brute) < 1e-12
for n, k in [(16, 4), (64, 16)]:
    X, S, T = planted_matrix(n, k, 20., check_rng)
    for threshold in [None, .8]:
        rows, cols = recovery_spectral(X, k, threshold)
        assert np.array_equal(S, rows) and np.array_equal(T, cols)
assert empirical_cutoff(list(range(99)), .05) == 94
assert empirical_cutoff(list(range(99)), .025) == 97
assert k_from_gamma(64, 2/3) == 16
print('PASS: exact scan and maximizing sets; both spectral variants; calibration; integer crossover.')
''', hidden=True, title='Verify the procedures')

# Plotting definitions are embedded, so Colab needs no project-specific imports.
code('''
LABELS = {'sum':'Sum', 'maximum':'Maximum entry', 'column_pool':'Column pooling',
          'scan':'Exact scan', 'sum_or_scan':'Sum or scan', 'exhaustive':'Maximizing block',
          'spectral':'Spectral', 'soft_spectral':'Thresholded spectral'}
COLORS = {'sum':'#a51c30', 'maximum':'#7b6baa', 'column_pool':'#ad8551',
          'scan':'#20232a', 'sum_or_scan':'#287b8e', 'exhaustive':'#20232a',
          'spectral':'#a51c30', 'soft_spectral':'#287b8e'}

def finish(fig, name):
    fig.tight_layout()
    fig.savefig(OUTPUT / (name + '.png'), dpi=140, bbox_inches='tight')
    plt.show()

def plot_curves(table, metric, pairs, name, x='mu', columns=4):
    rows = math.ceil(len(pairs)/columns)
    fig, axes = plt.subplots(rows, columns, figsize=(3.6*columns, 3.0*rows), squeeze=False)
    for ax, (n, k) in zip(axes.flat, pairs):
        subset = table[(table.n == n) & (table.k == k) & (table.metric == metric)]
        for method, group in subset.groupby('method'):
            group = group.sort_values(x)
            ax.plot(group[x], group.estimate, 'o-', ms=3, lw=1.4,
                    label=LABELS[method], color=COLORS[method])
            ax.fill_between(group[x].to_numpy(), group.ci_low.to_numpy(),
                            group.ci_high.to_numpy(), color=COLORS[method], alpha=.08)
        ax.set(title=f'n={n}, k={k}; gamma={math.log(k)/math.log(n):.2f}',
               xlabel=x, ylabel=metric.replace('_', ' '), ylim=(-.03, 1.03))
        ax.grid(alpha=.12)
    for ax in list(axes.flat)[len(pairs):]: ax.set_visible(False)
    handles = {}
    for ax in axes.flat:
        handles.update(zip(*reversed(ax.get_legend_handles_labels())))
    fig.legend(handles.values(), handles.keys(), loc='lower center',
               bbox_to_anchor=(.5,-.1), ncol=5, frameon=False)
    finish(fig, name)

def broad_heatmaps(table, plan, mapping, metric, methods, name):
    factors = [0., .65, 1., 1.35, 1.8]
    sizes = sorted(mapping.n.unique())
    gammas = sorted(mapping.gamma_requested.unique())
    view = signal_view(table, plan, 'scan')
    fig, axes = plt.subplots(len(methods), len(factors),
                             figsize=(14, 2.5*len(methods)), squeeze=False, layout='constrained')
    cmap = plt.get_cmap('Reds').copy(); cmap.set_bad('#dddddd')
    for i, method in enumerate(methods):
        for j, factor in enumerate(factors):
            ss = view[(view.metric == metric) & (view.method == method) & (view.parameter == factor)]
            lookup = {(r.n, r.k): r.estimate for r in ss.itertuples()}
            values = [[lookup.get((n, k_from_gamma(n, g)), np.nan) for g in gammas] for n in sizes]
            axes[i,j].imshow(values, vmin=0, vmax=1, aspect='auto', cmap=cmap)
            axes[i,j].set(xticks=range(len(gammas)),
                          xticklabels=[f'{g:g}' for g in gammas] if i == len(methods)-1 else [],
                          yticks=range(len(sizes)), yticklabels=sizes if j == 0 else [])
            if i == 0: axes[i,j].set_title(f'c = {factor:g}')
            if j == 0: axes[i,j].set_ylabel(LABELS[method] + '\\nn')
            if i == len(methods)-1: axes[i,j].set_xlabel('requested gamma')
    fig.suptitle(metric.replace('_',' ').capitalize() + ': every signal multiplier')
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1), cmap=cmap),
                 ax=axes.ravel().tolist(), fraction=.012, pad=.02)
    # Reserve a separate right margin for the shared colorbar.
    fig.savefig(OUTPUT / (name + '.png'), dpi=140, bbox_inches='tight')
    plt.show()

def threshold_curves(table, plan):
    refs = [('sum','power','sum'), ('scan','power','scan'), ('maximum','power','maximum'),
            ('exact_recovery','exact_recovery','exhaustive'),
            ('spectral','exact_recovery','spectral'), ('soft_spectral','exact_recovery','soft_spectral')]
    fig, axes = plt.subplots(2, 3, figsize=(14,7.5))
    for ax, (reference, metric, method) in zip(axes.flat, refs):
        ss = signal_view(table, plan, reference)
        ss = ss[(ss.metric == metric) & (ss.method == method)]
        for (n,k), group in ss.groupby(['n','k']):
            group = group.sort_values('parameter')
            ax.errorbar(group.parameter, group.estimate,
                        yerr=[np.maximum(0., group.estimate-group.ci_low),
                              np.maximum(0., group.ci_high-group.estimate)],
                        marker='o', capsize=2, label=f'n={n}, k={k}')
        ax.axvline(1, color='#777', ls=':', lw=1)
        ax.set(title=reference.replace('_',' '), xlabel='mu / reference scale',
               ylabel=metric.replace('_',' '), ylim=(-.03,1.03), xticks=[.5,1,2])
        ax.grid(alpha=.12); ax.legend(fontsize=8, frameon=False)
    finish(fig, 'threshold_comparisons')

def beta_curves(table, plan):
    view = signal_view(table, plan, 'beta')
    fig, axes = plt.subplots(3,4,figsize=(15,9),sharey=True)
    for ax, cfg in zip(axes.flat, plan):
        n,k = cfg['n'],cfg['k']; g = cfg['gamma_effective']
        ss = view[(view.n == n) & (view.k == k) & (view.metric == 'power')]
        for method, group in ss.groupby('method'):
            if method == 'column_pool': continue
            group = group.sort_values('parameter')
            ax.plot(group.parameter,group.estimate,'o-',ms=2,lw=1.2,
                    color=COLORS[method], label=LABELS[method])
        ax.axvline(max(g/2,2*g-1), color='#a51c30',ls=':',lw=1.3)
        ax.axvline(max(0,2*g-1), color='#777',ls='--',lw=1.1)
        ax.set(title=f'n={n}, k={k}; gamma={g:.3f}',xlabel='beta (larger = weaker signal)',ylim=(-.03,1.03))
        ax.grid(alpha=.12)
    axes[0,0].set_ylabel('Detection power')
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.035),ncol=4,frameon=False)
    finish(fig, 'power_law_regimes')
''', hidden=True, title='Plotting functions')

md(r'''
## 1. Sweep matrix size and block size
$n\in\{6,8,12,16,24,32,64\}$ and $\gamma\in\{.25,.40,.55,.70,.85\}$.
Set $k=\lfloor n^\gamma\rfloor$, clipped to $[2,n-1]$, and
$\mu=c\,2\sqrt{\log(n/k)/k}$ for $c\in\{0,.65,1,1.35,1.8\}$.
Rounding gives **31 distinct $(n,k)$ pairs** from 35 requested pairs; each runs once.

The table shows the actual exponent $\log k/\log n$, both regime classifications,
and which methods run. $\gamma=2/3$ separates the detection exponents;
$k=\sqrt n$ separates the spectral recovery regimes.
''')
code('''
broad_design, gamma_mapping = broad_plan()
display(gamma_mapping.pivot(index='n', columns='gamma_requested', values='k'))
display(plan_table(broad_design).round(3))
''')
code('''
# Every distinct (n,k) runs all five signals, for detection and recovery.
broad_summary, broad_trials, broad_calibration, broad_folder = run_study(
    'broad', broad_design, OUTPUT, seed=SEED,
    cal_reps=CALIBRATION_REPS, eval_reps=EVALUATION_REPS)
gamma_mapping.to_csv(broad_folder / 'gamma_mapping.csv', index=False)
''')
md(r'''
Each column below is a different signal multiplier; all five are shown.
Red means higher success (0–1); gray means the method was not run.
Repeated $k$ values in the mapping share a result. Column pooling remains in the saved tables.
''')
code('''
broad_heatmaps(broad_summary, broad_design, gamma_mapping, 'power',
               ['sum','maximum','scan','sum_or_scan'], 'broad_detection')
broad_heatmaps(broad_summary, broad_design, gamma_mapping, 'exact_recovery',
               ['exhaustive','spectral','soft_spectral'], 'broad_recovery')
''', hidden=True, title='Plot every broad-sweep signal level')
md(r'''
## 2. Cross the writeup's signal reference levels
The scan-scaled grid alone does not reach every method's useful signal range.
For $(n,k)=(12,2),(12,5),(64,4),(64,16)$, this study tests
$\mu\in\{.5a,a,2a\}$ for **each** applicable reference $a$, plus $\mu=0$.
Shared signal values are evaluated once.

| Reference $a$ | Writeup condition |
|---|---|
| $n/k^2$ | Sum: $\mu k^2/n\to\infty$ |
| $2\sqrt{\log(n/k)/k}$ | Scan: fixed separation above this scale |
| $2\sqrt{\log n}$ | Maximum: a sufficient signal scale |
| $\max\{(\sqrt{2\log k}+\sqrt{2\log(n-k)})/\sqrt k,\ 2\sqrt{\log(n/k)/k}\}$ | Exact block recovery boundary |
| $\sqrt n/k+\sqrt{\log n/k}$ | Spectral: sufficiently large constant; large-block regime |
| $\sqrt{\log(n/k^2)}$ | Thresholded spectral: $\mu/a\to\infty$; small-block regime |

These are asymptotic references, not finite-sample cutoffs or guarantees at multiplier 1.
Exact scanning runs for the two $n=12$ cases. Spectral recovery runs throughout;
thresholded spectral runs for $(12,2)$ and $(64,4)$.
''')
code('''
threshold_design = threshold_plan()
display(plan_table(threshold_design).round(3))
# This checks below/at/above coverage for every listed reference in every pair.
for cfg in threshold_design:
    for reference in reference_scales(cfg['n'], cfg['k']):
        assert {s['parameter'] for s in cfg['signals'] if s['reference'] == reference} == {.5,1.,2.}
threshold_summary, threshold_trials, threshold_calibration, threshold_folder = run_study(
    'thresholds', threshold_design, OUTPUT, seed=SEED,
    cal_reps=CALIBRATION_REPS, eval_reps=EVALUATION_REPS, stream_keys=(301,302))
''')
code('threshold_curves(threshold_summary, threshold_design)', hidden=True,
     title='Plot the signal-threshold comparisons')
md(r'''
### Check the theoretical rejection thresholds
The plots above use calibrated 5% tests. Here we reuse Study 2's held-out matrices
with the writeup's thresholds: $c_n=\mu k^2/(2n)$ for the sum,
$h_n=\sqrt{2\log\binom nk^2}$ for the scan, and
$u_n=\sqrt{(4+.5)\log n}$ for the maximum.
We plot false-positive rate plus miss rate, the writeup's total testing error.
The sum's exact error is $2\Phi(-\mu k^2/(2n))$.
''')
code('''
risk = theory_errors(threshold_trials)
risk.to_csv(threshold_folder / 'theory_threshold_errors.csv', index=False)
fig, axes = plt.subplots(1,4,figsize=(15,3.3))
for ax, cfg in zip(axes, threshold_design):
    n,k = cfg['n'],cfg['k']
    for method, group in risk[(risk.n == n) & (risk.k == k)].groupby('method'):
        group = group.sort_values('mu')
        ax.plot(group.mu,group.total_error,'o-',ms=2,color=COLORS[method],label=LABELS[method])
        if method == 'sum':
            ax.plot(group.mu,group.exact_sum_risk,ls=':',color=COLORS[method],label='Sum: exact error')
    ax.set(title=f'n={n}, k={k}',xlabel='mu',ylabel='Total testing error',ylim=(-.03,1.15))
    ax.grid(alpha=.12)
handles,labels = axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.1),ncol=4,frameon=False)
finish(fig, 'theoretical_threshold_errors')
''', hidden=True, title='Evaluate the writeup rejection thresholds')
md(r'''
## 3. Sweep the power-law regimes
Here $\mu=n^{-\beta}$, with $n\in\{12,24,64\}$ and requested
$\gamma\in\{.4,.5,2/3,.8\}$. Every pair tests
$\beta\in\{0,.1,.2,.3,.4,.6,.8\}$, the two exponent boundaries,
and points .04 to either side when nonnegative. A separate $\mu=0$ check is included;
$\beta=0$ means $\mu=1$.

The plots mark $\beta_{\rm stat}=\max\{\gamma/2,2\gamma-1\}$ (dotted red)
and $\beta_{\rm comp}=\max\{0,2\gamma-1\}$ (dashed gray), using the **actual** exponent
after rounding $k$. This study runs detection only. Logarithmic factors and small $n$
can obscure exponent predictions, especially at equality.
''')
code('''
beta_design = power_law_plan()
display(plan_table(beta_design).round(3))
beta_summary, beta_trials, beta_calibration, beta_folder = run_study(
    'power_laws', beta_design, OUTPUT, seed=SEED,
    cal_reps=CALIBRATION_REPS, eval_reps=EVALUATION_REPS,
    recovery=False, stream_keys=(401,402))
''')
code('beta_curves(beta_summary, beta_design)', hidden=True, title='Plot power-law regimes')
md(r'''
## Coverage and calibration checks
The runner checks every expected $(n,k,\mu,\text{method})$, trial count, and reported
success count against the raw trials. Checkpoints are reused only for the same design;
the notebook reports each verified checkpoint. Sum and maximum-entry results are also
checked against their exact Gaussian rejection probabilities.
''')
code('''
study_frames = {'broad': broad_summary, 'thresholds': threshold_summary, 'power_laws': beta_summary}
for name, summary, calibration, folder in [
    ('broad', broad_summary, broad_calibration, broad_folder),
    ('thresholds', threshold_summary, threshold_calibration, threshold_folder),
    ('power_laws', beta_summary, beta_calibration, beta_folder)]:
    analytic = analytic_detection_check(summary, calibration)
    analytic.to_csv(folder / 'analytic_sum_max_check.csv', index=False)
    print(f'PASS: {name}: {len(analytic)} sum/max rates agree with exact Gaussian probabilities.')
fpr = pd.concat([frame[frame.metric == 'false_positive_rate'].assign(study=name)
                 for name,frame in study_frames.items()])
display(fpr.groupby(['study','method']).agg(
    signal_conditions=('estimate','size'), mean_fpr=('estimate','mean'),
    min_fpr=('estimate','min'), max_fpr=('estimate','max')).round(3))
# Each entry is a separately evaluated signal condition, sharing its pair's calibration cutoff.
fig, ax = plt.subplots(figsize=(10,3))
for i, (method, group) in enumerate(fpr.groupby('method')):
    ax.scatter([i]*len(group), group.estimate, s=9, alpha=.3, color=COLORS[method])
ax.set(xticks=range(fpr.method.nunique()), xticklabels=sorted(fpr.method.unique()),
       ylabel='Held-out false-positive rate')
ax.axhline(.05,color='#a51c30',ls='--',lw=1)
finish(fig,'false_positive_rates')
''', hidden=True, title='Check held-out false-positive rates')
md(r'''
| Writeup topic | Evidence here | Limit |
|---|---|---|
| Sum, scan, maximum, sum-or-scan | Studies 1–3; held-out false positives and power | Scan omitted when exact search exceeds the cap |
| Sum/scan/max error calculations | Study 2 theoretical-threshold errors; exact sum curve | Finite sample error, not an asymptotic proof |
| Exact block recovery | Studies 1–2, maximizing indices | Small feasible searches only |
| Spectral recovery for large blocks | Studies 1–2, including $(64,16)$ | Full SVD and known $k$ |
| Thresholded spectral for small blocks | Studies 1–2, including $(64,4)$ | Fixed threshold; no outcome-based tuning |
| Detection regimes and crossovers | Studies 1 and 3; requested $\gamma=.5,2/3$ in Study 3 | Integer rounding is reported |
| Statistical impossibility and computational lower bounds | Exponent references in Study 3 | Simulations cannot verify impossibility or the clique hypotheses |

Higher criticism is mentioned in the writeup but is not one of its implemented procedures.
The computational lower bounds concern broader model classes; algorithm failures here
do not establish those lower bounds.
''')
md(r'''
## Saved data
Each study saves its complete design, raw trials, calibrated cutoffs, configuration table,
and summaries. The ZIP also contains the figures. Timing columns measure this implementation
on this machine; combined-test timing is not recorded. Colab storage is temporary, so download
the ZIP to retain a run.
''')
code('''
bundle = Path('STAT236_notebook_results.zip')
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(OUTPUT.rglob('*')):
        if path.is_file():
            archive.write(path,path.relative_to(OUTPUT.parent).as_posix())
display(pd.DataFrame([
    {'study':'broad', 'configurations':len(broad_design), 'trials':len(broad_trials)},
    {'study':'thresholds', 'configurations':len(threshold_design), 'trials':len(threshold_trials)},
    {'study':'power_laws', 'configurations':len(beta_design), 'trials':len(beta_trials)}]))
print('Saved:', bundle.resolve())
''')
md(r'''
**Sources:** [Butucea & Ingster (2013)](https://doi.org/10.3150/12-BEJ470),
[Ma & Wu (2015)](https://doi.org/10.1214/14-AOS1300),
[Butucea, Ingster & Suslina (2015)](https://doi.org/10.1051/ps/2014017),
[Cai, Liang & Rakhlin (2017)](https://doi.org/10.1214/16-AOS1488).
''')

notebook = nbf.v4.new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name':'Python 3','language':'python','name':'python3'},
    'language_info': {'name':'python','version':'3.11'},
    'colab': {'name':'STAT236_Sparse_Submatrix_Experiments.ipynb','toc_visible':True},
})
nbf.write(notebook, ROOT / 'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb')
print(f'Built {len(cells)} cells; no project imports are needed when running the notebook.')
