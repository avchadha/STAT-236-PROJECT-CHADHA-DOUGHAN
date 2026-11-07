#!/usr/bin/env python3
"""Reproducible Gaussian square-submatrix detection and exact-recovery experiments.

Requires numpy, scipy, matplotlib. All paper statistics are generated from
independent, seeded Monte Carlo samples, NOT from website visitor simulations.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import platform
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
DEFAULT_SEED = 2362026
ALPHA = 0.05

DETECTION = [
    dict(id='D1', n=16, k=4, mus=[0, .4, .8, 1.2, 1.6, 2., 2.4], scan=True),
    dict(id='D2', n=12, k=5, mus=[0, .2, .4, .6, .8, 1., 1.2], scan=True),
    dict(id='D3', n=64, k=16, mus=[0, .15, .3, .45, .6, .75, .9], scan=False),
    dict(id='D4', n=64, k=4, mus=[0, .5, 1., 1.5, 2., 2.5, 3.], scan=False),
]
RECOVERY = [
    dict(id='R1', n=16, k=4, mus=[0, .75, 1., 1.25, 1.5, 1.75, 2., 2.5], exhaustive=True, thresholded=False),
    dict(id='R2', n=12, k=5, mus=[0, .35, .5, .65, .8, 1., 1.2], exhaustive=True, thresholded=False),
    dict(id='R3', n=64, k=16, mus=[0, .4, .6, .8, 1., 1.2, 1.5], exhaustive=False, thresholded=False),
    dict(id='R4', n=64, k=4, mus=[0, 1., 2., 3., 4., 5.], exhaustive=False, thresholded=True),
]


def wilson(successes: int, total: int, z: float = 1.95996398454) -> tuple[float, float]:
    """Two-sided 95% Wilson score interval (including 0 and 1 cases)."""
    p = successes / total
    d = 1 + z*z/total
    center = (p + z*z/(2*total)) / d
    width = z*math.sqrt((p*(1-p)+z*z/(4*total))/total) / d
    return max(0., center-width), min(1., center+width)


def planted_matrix(n: int, k: int, mu: float, rng: np.random.Generator):
    X = rng.standard_normal((n, n))
    S = np.sort(rng.choice(n, k, replace=False))
    T = np.sort(rng.choice(n, k, replace=False))
    X[np.ix_(S, T)] += mu
    return X, S, T


ROW_COMBOS: dict[tuple[int, int], np.ndarray] = {}


def exact_scan(X: np.ndarray, k: int, *, return_indices: bool = False):
    """Exactly maximize sum over ALL k x k row/column index subsets.

    Row subsets are enumerated. For each, maximization over column subsets
    is solved exactly by choosing its k largest aggregated column sums.
    The max is not restricted to the planted location.
    """
    n = X.shape[0]
    key = (n, k)
    if key not in ROW_COMBOS:
        if math.comb(n, k) > 100_000:
            raise ValueError('Exhaustive scan disabled above 100,000 row subsets')
        ROW_COMBOS[key] = np.array(list(itertools.combinations(range(n), k)), dtype=np.int16)
    combos = ROW_COMBOS[key]
    best_value = -np.inf
    best_rows = best_cols = None
    # Batching ensures bounded RAM even for moderately large candidate counts.
    for start in range(0, len(combos), 512):
        rowsets = combos[start:start+512]
        column_sums = X[rowsets].sum(axis=1)
        # Deterministic tie handling; exact continuous Gaussian data ties have probability zero.
        top = np.argsort(column_sums, axis=1)[:, -k:]
        values = np.take_along_axis(column_sums, top, axis=1).sum(axis=1)
        idx = int(np.argmax(values))
        if values[idx] > best_value:
            best_value = float(values[idx])
            best_rows = rowsets[idx].copy()
            best_cols = np.sort(top[idx]).astype(np.int16)
    if return_indices:
        return best_value / k, best_rows, best_cols
    return best_value / k


def test_statistics(X: np.ndarray, k: int, include_scan: bool):
    n = len(X)
    out = {
        'sum': float(X.sum()/n),
        'maximum': float(np.max(X)),
        'column_pool': float(X.sum(axis=0).max()/math.sqrt(n)),
    }
    if include_scan:
        out['scan'] = exact_scan(X, k)
    return out


def empirical_cutoff(values: list[float], alpha: float):
    """Conservative split-calibration order statistic.

    Under the continuous IID null, a fresh null statistic exceeds the
    m-th ascending order statistic with probability (N+1-m)/(N+1).
    """
    N = len(values)
    rank = math.ceil((N+1)*(1-alpha))  # 1-based
    if rank > N:
        raise ValueError('Not enough calibration draws for requested alpha')
    return float(np.sort(values)[rank-1])


def detection_decisions(stats: dict, critical: dict):
    flags = {name: int(val > critical[name]) for name, val in stats.items()}
    if 'scan' in stats:
        flags['sum_or_scan'] = int(stats['sum'] > critical['sum_bonferroni'] or stats['scan'] > critical['scan_bonferroni'])
    return flags


def recovery_spectral(X: np.ndarray, k: int, threshold: float | None = None):
    if threshold is None:
        Y = X
    else:
        Y = np.sign(X) * np.maximum(np.abs(X)-threshold, 0.)
    u, _, vt = np.linalg.svd(Y, full_matrices=False)
    u = u[:, 0].copy()
    v = vt[0, :].copy()
    # The rank-one mean is positive; orient both singular vectors jointly.
    if float(u.sum()) < 0.:
        u *= -1.
        v *= -1.
    row_scores = Y @ v
    col_scores = Y.T @ u
    rows = np.sort(np.argsort(row_scores)[-k:])
    cols = np.sort(np.argsort(col_scores)[-k:])
    return rows, cols


def jaccard(a: np.ndarray, b: np.ndarray):
    match = np.intersect1d(a,b,assume_unique=True).size
    return float(match / (2*len(a)-match))


def export_csv(path: Path, rows: list[dict]):
    if not rows:
        raise ValueError('Refusing to write empty CSV')
    fields = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def run_detection(scenario, rng, cal_reps, reps, summary, trials, cutoffs):
    n, k, scan = scenario['n'], scenario['k'], scenario['scan']
    keys = ['sum','maximum','column_pool'] + (['scan'] if scan else [])
    calibrations = {method: [] for method in keys}
    for c in range(cal_reps):
        X = rng.standard_normal((n,n))
        ss = test_statistics(X,k,scan)
        for name in keys:
            calibrations[name].append(ss[name])
    critical = {name: empirical_cutoff(vals, ALPHA) for name, vals in calibrations.items()}
    if scan:
        critical['sum_bonferroni'] = empirical_cutoff(calibrations['sum'], ALPHA/2)
        critical['scan_bonferroni'] = empirical_cutoff(calibrations['scan'], ALPHA/2)
    for name, val in critical.items():
        cutoffs.append(dict(scenario=scenario['id'], n=n, k=k, statistic=name, calibration_reps=cal_reps,
                            alpha=(ALPHA/2 if 'bonferroni' in name else ALPHA), critical_value=val))
    all_methods = keys + (['sum_or_scan'] if scan else [])
    for mu in scenario['mus']:
        buckets = {name: {'null':0,'signal':0} for name in all_methods}
        runtime = {name: {'null':0.,'signal':0.} for name in keys}
        for condition in ('null','signal'):
            for rep in range(reps):
                if condition == 'null':
                    X = rng.standard_normal((n,n))
                else:
                    X, _, _ = planted_matrix(n,k,mu,rng)
                # All methods see the EXACT same matrix in each replicate.
                stats = {}
                for name in keys:
                    start = time.perf_counter()
                    if name == 'sum': value = float(X.sum()/n)
                    elif name == 'maximum': value = float(X.max())
                    elif name == 'column_pool': value = float(X.sum(axis=0).max()/math.sqrt(n))
                    else: value = exact_scan(X,k)
                    runtime[name][condition] += time.perf_counter() - start
                    stats[name] = value
                flags = detection_decisions(stats, critical)
                for name in all_methods:
                    buckets[name][condition] += flags[name]
                trial = dict(task='detection',scenario=scenario['id'],n=n,k=k,mu=mu,
                             condition=condition,rep=rep)
                for name,value in stats.items(): trial['stat_'+name] = round(value,10)
                for name,flag in flags.items(): trial['reject_'+name] = flag
                trials.append(trial)
        for name in all_methods:
            for condition in ('null','signal'):
                count = buckets[name][condition]
                low, high = wilson(count,reps)
                summary.append(dict(task='detection',scenario=scenario['id'],n=n,k=k,mu=mu,method=name,
                                   metric=('false_positive_rate' if condition=='null' else 'power'),
                                   count=count,replicates=reps,estimate=count/reps,ci_low=low,ci_high=high,
                                   alpha=ALPHA,calibration_reps=cal_reps,
                                   mean_eval_seconds=(runtime[name][condition]/reps if name in runtime else None)))
        print(f"  {scenario['id']} mu={mu:.2f} power: " + ', '.join(f"{x}:{buckets[x]['signal']/reps:.2f}" for x in all_methods), flush=True)


def run_recovery(scenario, rng, reps, summary, trials):
    n, k = scenario['n'], scenario['k']
    methods = (['exhaustive'] if scenario['exhaustive'] else [])+['spectral'] + (['soft_spectral'] if scenario['thresholded'] else [])
    threshold = math.sqrt(2 * math.log(n/(k*k))) if scenario['thresholded'] else None
    for mu in scenario['mus']:
        exact = {name:0 for name in methods}
        overlap = {name:0. for name in methods}
        runtime = {name:0. for name in methods}
        for rep in range(reps):
            X, S, T = planted_matrix(n,k,mu,rng)
            trial = dict(task='recovery',scenario=scenario['id'],n=n,k=k,mu=mu,condition='signal',rep=rep)
            for name in methods:
                start=time.perf_counter()
                if name == 'exhaustive':
                    _, sr, tc = exact_scan(X,k,return_indices=True)
                elif name=='spectral':
                    sr, tc = recovery_spectral(X,k)
                else:
                    sr, tc = recovery_spectral(X,k,threshold=threshold)
                runtime[name] += time.perf_counter()-start
                ok = int(np.array_equal(S,sr) and np.array_equal(T,tc))
                part=(jaccard(S,sr)+jaccard(T,tc))/2
                exact[name]+=ok
                overlap[name]+=part
                trial['exact_'+name]=ok
                trial['jaccard_'+name]=round(part,8)
            trials.append(trial)
        for name in methods:
            lo,hi=wilson(exact[name],reps)
            summary.append(dict(task='recovery',scenario=scenario['id'],n=n,k=k,mu=mu,method=name,
                                metric='exact_recovery',count=exact[name],replicates=reps,estimate=exact[name]/reps,
                                ci_low=lo,ci_high=hi,alpha=None,calibration_reps=None,
                                mean_eval_seconds=runtime[name]/reps))
            summary.append(dict(task='recovery',scenario=scenario['id'],n=n,k=k,mu=mu,method=name,
                                metric='mean_index_jaccard',count=None,replicates=reps,estimate=overlap[name]/reps,
                                ci_low=None,ci_high=None,alpha=None,calibration_reps=None,
                                mean_eval_seconds=runtime[name]/reps))
        print(f"  {scenario['id']} mu={mu:.2f} exact: " + ', '.join(f"{x}:{exact[x]/reps:.2f}" for x in methods), flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['quick','full'],default='full')
    parser.add_argument('--seed',type=int,default=DEFAULT_SEED)
    parser.add_argument('--only',default='',help='Optional scenario ID, e.g. D1 or R2')
    args=parser.parse_args()
    os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
    # Unlike setting these at import time, this may not affect an already-loaded BLAS.
    cal_reps, reps = (250, 100) if args.mode == 'quick' else (1500,600)
    if args.only and args.only not in [s['id'] for s in DETECTION+RECOVERY]:
        parser.error('Unknown scenario ID')
    RESULTS.mkdir(parents=True,exist_ok=True)
    # Separate independent deterministic streams by scenario, so --only D2
    # reproduces the D2 block of a full run byte-for-byte (timings aside).
    def stream(scenario_id):
        return np.random.default_rng(np.random.SeedSequence(
            [args.seed, ord(scenario_id[0]), int(scenario_id[1:])]))
    summary, trials, cutoffs = [],[],[]
    start=time.perf_counter()
    for scenario in DETECTION:
        if args.only and scenario['id']!=args.only: continue
        print('Starting',scenario['id'],'detection',flush=True)
        run_detection(scenario,stream(scenario['id']),cal_reps,reps,summary,trials,cutoffs)
    for scenario in RECOVERY:
        if args.only and scenario['id']!=args.only: continue
        print('Starting',scenario['id'],'recovery',flush=True)
        run_recovery(scenario,stream(scenario['id']),reps,summary,trials)
    name = '' if not args.only else '_'+args.only
    export_csv(RESULTS/f'summary{name}.csv',summary)
    export_csv(RESULTS/f'trials{name}.csv',trials)
    if cutoffs: export_csv(RESULTS/f'calibration{name}.csv',cutoffs)
    meta = dict(seed=args.seed, alpha=ALPHA, mode=args.mode, calibration_reps=cal_reps,
                null_validation_reps_per_signal=reps, alternative_reps_per_signal=reps,
                recovery_reps_per_signal=reps, time_seconds=time.perf_counter()-start,
                scenario_rng='SeedSequence([master_seed, ord(task_letter), scenario_index])',
                numpy_version=np.__version__, platform=platform.platform(), python=platform.python_version(),
                detection_scenarios=DETECTION, recovery_scenarios=RECOVERY,
                note=('The same held-out null matrices are redrawn independently for every mu; '+
                      'the same planted alternative is shared across procedures within a replicate. '+
                      'The exhaustive scan is exact, with arbitrary nonconsecutive supports.'))
    (RESULTS/f'metadata{name}.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Finished in',round(meta['time_seconds'],2),'seconds; wrote',len(summary),'summary rows and',len(trials),'trial rows.',flush=True)

if __name__=='__main__':
    main()