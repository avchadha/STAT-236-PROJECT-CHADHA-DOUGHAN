"""Fixed experiment designs embedded into the standalone notebook by build_notebook.py."""
import hashlib
import io
import json
import math
import platform
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd

from run_experiments import run_detection, run_recovery


def read_results(path):
    # pandas otherwise interprets the experimental condition "null" as missing data.
    return pd.read_csv(path, keep_default_na=False, na_values=[''])


def k_from_gamma(n, gamma):
    value = n ** gamma
    # For example, floating-point 64**(2/3) is just below the integer 16.
    if math.isclose(value, round(value), rel_tol=1e-12):
        value = round(value)
    return max(2, min(n - 1, math.floor(value)))


def reference_scales(n, k):
    scan = 2 * math.sqrt(math.log(n / k) / k)
    scales = {
        'sum': n / k**2,
        'scan': scan,
        'maximum': 2 * math.sqrt(math.log(n)),
        'exact_recovery': max(
            (math.sqrt(2 * math.log(k)) + math.sqrt(2 * math.log(n - k))) / math.sqrt(k),
            scan),
        'spectral': math.sqrt(n) / k + math.sqrt(math.log(n) / k),
    }
    if n > k**2:
        scales['soft_spectral'] = math.sqrt(math.log(n / k**2))
    return scales


def make_config(n, k, signals, requested_gammas=()):
    signals = [dict(s, mu=round(float(s['mu']), 12)) for s in signals]
    return dict(n=n, k=k, gamma_effective=math.log(k) / math.log(n),
                requested_gammas=list(requested_gammas), signals=signals,
                mus=list(dict.fromkeys(s['mu'] for s in signals)))


def broad_plan():
    sizes = [6, 8, 12, 16, 24, 32, 64]
    gammas = [.25, .40, .55, .70, .85]
    mapping = [dict(n=n, gamma_requested=g, k=k_from_gamma(n, g))
               for n in sizes for g in gammas]
    plan = []
    for n, k in sorted({(r['n'], r['k']) for r in mapping}):
        scale = reference_scales(n, k)['scan']
        signals = [dict(reference='scan', parameter=c, mu=c * scale)
                   for c in [0., .65, 1., 1.35, 1.8]]
        targets = [r['gamma_requested'] for r in mapping if (r['n'], r['k']) == (n, k)]
        plan.append(make_config(n, k, signals, targets))
    return plan, pd.DataFrame(mapping)


def threshold_plan():
    # Small cases allow exact scan; larger cases compare efficient methods.
    plan = []
    for n, k in [(12, 2), (12, 5), (64, 4), (64, 16)]:
        signals = [dict(reference='null', parameter=0., mu=0.)]
        for name, scale in reference_scales(n, k).items():
            for factor in [.5, 1., 2.]:
                signals.append(dict(reference=name, parameter=factor, mu=factor * scale))
        plan.append(make_config(n, k, signals))
    return plan


def power_law_plan():
    plan = []
    for n in [12, 24, 64]:
        for gamma in [.4, .5, 2 / 3, .8]:
            k = k_from_gamma(n, gamma)
            effective = math.log(k) / math.log(n)
            stat = max(effective / 2, 2 * effective - 1)
            comp = max(0., 2 * effective - 1)
            betas = {0., .1, .2, .3, .4, .6, .8}
            # Include points on either side and at each attainable exponent boundary.
            for boundary in [stat, comp]:
                betas.update(round(b, 12) for b in [boundary - .04, boundary, boundary + .04]
                             if 0 <= b < 1)
            signals = [dict(reference='null', parameter=0., mu=0.)]
            signals += [dict(reference='beta', parameter=b, mu=n**(-b)) for b in sorted(betas)]
            plan.append(make_config(n, k, signals, [gamma]))
    return plan


def plan_table(plan, exact_cap=2500):
    rows = []
    for cfg in plan:
        n, k = cfg['n'], cfg['k']
        g = cfg['gamma_effective']
        rows.append(dict(n=n, k=k, gamma_effective=g, k_over_sqrt_n=k / math.sqrt(n),
                         detection_regime=('crossover' if math.isclose(g, 2/3) else
                                           'highly sparse' if g < 2/3 else 'moderately sparse'),
                         recovery_regime=('crossover' if k*k == n else
                                          'small block' if k*k < n else 'large block'),
                         exact_scan=math.comb(n, k) <= exact_cap,
                         soft_spectral=n > k*k, signal_values=len(cfg['mus']),
                         mu_min=min(cfg['mus']), mu_max=max(cfg['mus'])))
    return pd.DataFrame(rows)


def analytic_detection_check(summary, calibration):
    """Compare simulated sum/max rates with exact Gaussian rejection probabilities."""
    critical = calibration.set_index(['n', 'k', 'statistic']).critical_value
    rows = []

    def log_cdf(x):
        if x >= 0:
            return math.log1p(-.5 * math.erfc(x / math.sqrt(2)))
        probability = .5 * math.erfc(-x / math.sqrt(2))
        return math.log(probability) if probability else -math.inf

    for r in summary.itertuples():
        if r.task != 'detection' or r.method not in {'sum', 'maximum'}:
            continue
        cutoff = critical.loc[(r.n, r.k, r.method)]
        mu = r.mu if r.metric == 'power' else 0.
        if r.method == 'sum':
            probability = .5 * math.erfc((cutoff - mu*r.k*r.k/r.n) / math.sqrt(2))
        else:
            log_no_rejection = ((r.n*r.n-r.k*r.k) * log_cdf(cutoff) +
                                r.k*r.k * log_cdf(cutoff-mu))
            probability = -math.expm1(log_no_rejection)
        tolerance = 6 * math.sqrt(probability * (1-probability) / r.replicates) + .02
        assert abs(r.estimate - probability) <= tolerance
        rows.append(dict(n=r.n, k=r.k, mu=r.mu, method=r.method, metric=r.metric,
                         observed=r.estimate, exact_probability=probability))
    return pd.DataFrame(rows)


def validate_study(plan, summary, calibration, trials, cal_reps, eval_reps,
                   recovery=True, exact_cap=2500):
    """Reject missing, duplicated, stale, or incompletely written configurations."""
    expected_pairs = {(c['n'], c['k']) for c in plan}
    assert set(zip(summary.n, summary.k)) == expected_pairs
    keys = ['task', 'n', 'k', 'mu', 'method', 'metric']
    assert not summary.duplicated(keys).any()
    assert summary.estimate.between(0, 1).all()
    assert summary.replicates.eq(eval_reps).all()
    expected_summary_rows = 0
    expected_calibration_rows = 0
    expected_trials = 0
    for cfg in plan:
        n, k = cfg['n'], cfg['k']
        exact = math.comb(n, k) <= exact_cap
        methods = {'sum', 'maximum', 'column_pool'} | ({'scan', 'sum_or_scan'} if exact else set())
        rec_methods = {'spectral'} | ({'exhaustive'} if exact else set())
        if n > k*k:
            rec_methods.add('soft_spectral')
        rows = summary[(summary.n == n) & (summary.k == k)]
        for task, names, metrics in [
            ('detection', methods, {'power', 'false_positive_rate'}),
            ('recovery', rec_methods if recovery else set(),
             {'exact_recovery', 'mean_index_jaccard'} if recovery else set())]:
            group = rows[rows.task == task]
            expected = {(mu, m, metric) for mu in cfg['mus'] for m in names for metric in metrics}
            assert set(zip(group.mu.round(12), group.method, group.metric)) == expected
        cut = calibration[(calibration.n == n) & (calibration.k == k)]
        cal_names = methods - {'sum_or_scan'}
        if exact:
            cal_names |= {'sum_bonferroni', 'scan_bonferroni'}
        assert set(cut.statistic) == cal_names and len(cut) == len(cal_names)
        assert cut.calibration_reps.eq(cal_reps).all()
        expected_summary_rows += len(cfg['mus']) * (2*len(methods) + (2*len(rec_methods) if recovery else 0))
        expected_calibration_rows += len(cal_names)
        expected_trials += len(cfg['mus']) * eval_reps * (3 if recovery else 2)
    assert len(summary) == expected_summary_rows and len(calibration) == expected_calibration_rows
    if trials is not None:
        assert len(trials) == expected_trials
        assert not trials.duplicated(['task', 'n', 'k', 'mu', 'condition', 'rep']).any()
        expected = {(task, c['n'], c['k'], mu, condition): eval_reps
                    for c in plan for mu in c['mus']
                    for task, condition in [('detection', 'null'), ('detection', 'signal')] +
                    ([('recovery', 'signal')] if recovery else [])}
        observed = trials.assign(mu=trials.mu.round(12)).groupby(
            ['task', 'n', 'k', 'mu', 'condition']).size().to_dict()
        assert observed == expected
        # Aggregate once so verification stays inexpensive for the full sweep.
        grouped = trials.assign(mu=trials.mu.round(12)).groupby(
            ['task', 'n', 'k', 'mu', 'condition'])
        raw_counts = grouped[[c for c in trials if c.startswith(('reject_', 'exact_'))]].sum()
        raw_overlaps = grouped[[c for c in trials if c.startswith('jaccard_')]].mean()
        for r in summary.itertuples():
            if r.metric == 'mean_index_jaccard':
                actual = raw_overlaps.loc[(r.task, r.n, r.k, round(r.mu, 12), 'signal'),
                                          'jaccard_' + r.method]
                assert abs(actual - r.estimate) < 1e-8
                continue
            condition = 'null' if r.metric == 'false_positive_rate' else 'signal'
            field = ('reject_' if r.task == 'detection' else 'exact_') + r.method
            assert raw_counts.loc[(r.task, r.n, r.k, round(r.mu, 12), condition), field] == r.count
    return expected_trials


def run_study(name, plan, output, seed=2362026, cal_reps=1500, eval_reps=600,
              recovery=True, stream_keys=(101, 202), exact_cap=2500):
    """Execute every listed signal; resume only checkpoints matching this exact design."""
    design = dict(version=2, name=name, seed=seed, alpha=.05, cal_reps=cal_reps,
                  eval_reps=eval_reps, recovery=recovery, stream_keys=stream_keys,
                  exact_cap=exact_cap, plan=plan)
    tag = hashlib.sha256(json.dumps(design, sort_keys=True).encode()).hexdigest()[:12]
    folder = Path(output) / name / tag
    folder.mkdir(parents=True, exist_ok=True)
    summaries, raw_trials, cutoffs = [], [], []
    for index, cfg in enumerate(plan, 1):
        n, k = cfg['n'], cfg['k']
        prefix = folder / f'n{n}_k{k}'
        paths = {key: Path(str(prefix) + '_' + key + '.csv')
                 for key in ['summary', 'trials', 'calibration']}
        marker = Path(str(prefix) + '_complete.json')
        cached = marker.exists() and all(p.exists() for p in paths.values())
        if cached:
            frames = {key: read_results(p) for key, p in paths.items()}
            validate_study([cfg], frames['summary'], frames['calibration'], frames['trials'],
                           cal_reps, eval_reps, recovery, exact_cap)
        else:
            summary, trials, calibration = [], [], []
            exact = math.comb(n, k) <= exact_cap
            detection = dict(id=f'D_{n}_{k}', n=n, k=k, mus=cfg['mus'], scan=exact)
            localization = dict(id=f'R_{n}_{k}', n=n, k=k, mus=cfg['mus'],
                                exhaustive=exact, thresholded=n > k*k)
            rng_d = np.random.default_rng(np.random.SeedSequence([seed, stream_keys[0], n, k]))
            rng_r = np.random.default_rng(np.random.SeedSequence([seed, stream_keys[1], n, k]))
            with redirect_stdout(io.StringIO()):
                run_detection(detection, rng_d, cal_reps, eval_reps, summary, trials, calibration)
                if recovery:
                    run_recovery(localization, rng_r, eval_reps, summary, trials)
            frames = dict(summary=pd.DataFrame(summary), trials=pd.DataFrame(trials),
                          calibration=pd.DataFrame(calibration))
            validate_study([cfg], frames['summary'], frames['calibration'], frames['trials'],
                           cal_reps, eval_reps, recovery, exact_cap)
            for key, frame in frames.items():
                tmp = paths[key].with_suffix('.tmp')
                frame.to_csv(tmp, index=False)
                tmp.replace(paths[key])
            marker.write_text(json.dumps(dict(design=tag)), encoding='utf-8')
        summaries.append(frames['summary'])
        raw_trials.append(frames['trials'])
        cutoffs.append(frames['calibration'])
        print(f'{name} {index}/{len(plan)}: n={n}, k={k}, '
              f'{len(cfg["mus"])} shifts, exact={math.comb(n,k)<=exact_cap}'
              + (' (verified checkpoint)' if cached else ''), flush=True)
    summary = pd.concat(summaries, ignore_index=True)
    trials = pd.concat(raw_trials, ignore_index=True)
    calibration = pd.concat(cutoffs, ignore_index=True)
    total = validate_study(plan, summary, calibration, trials, cal_reps, eval_reps, recovery, exact_cap)
    for key, frame in [('summary', summary), ('trials', trials), ('calibration', calibration),
                       ('configurations', plan_table(plan, exact_cap))]:
        frame.to_csv(folder / f'{key}.csv', index=False)
    (folder / 'design.json').write_text(json.dumps(design, indent=2) + '\n', encoding='utf-8')
    environment = dict(numpy=np.__version__, python=platform.python_version(),
                       platform=platform.platform(), seed=seed,
                       seed_scheme='SeedSequence([seed, task_stream_key, n, k])',
                       note='Timing and SVD/score tie handling can vary across numerical environments.')
    (folder / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n', encoding='utf-8')
    print(f'PASS: all {len(plan)} configurations and {total:,} trial records verified.')
    return summary, trials, calibration, folder


def signal_view(summary, plan, reference):
    # Join by actual mu: a shift shared by reference grids is simulated only once.
    design = pd.DataFrame([dict(n=c['n'], k=c['k'], **s)
                           for c in plan for s in c['signals'] if s['reference'] == reference])
    table = summary.copy()
    table['mu'] = table.mu.round(12)
    return table.merge(design, on=['n', 'k', 'mu'], validate='many_to_many')


def theory_errors(trials, epsilon=.5):
    """Evaluate the writeup's thresholds on the threshold-study matrices."""
    rows = []
    for (n, k, mu), group in trials[trials.task == 'detection'].groupby(['n', 'k', 'mu']):
        log_candidates = 2 * (math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1))
        thresholds = {'sum': mu*k*k/(2*n), 'scan': math.sqrt(2*log_candidates),
                      'maximum': math.sqrt((4+epsilon)*math.log(n))}
        for method, cutoff in thresholds.items():
            field = 'stat_' + method
            if field not in group or group[field].isna().all():
                continue
            null = group.loc[group.condition == 'null', field]
            signal = group.loc[group.condition == 'signal', field]
            fp, miss = float((null > cutoff).mean()), float((signal <= cutoff).mean())
            rows.append(dict(n=n, k=k, mu=mu, method=method, cutoff=cutoff,
                             false_positive_rate=fp, miss_rate=miss, total_error=fp+miss,
                             exact_sum_risk=(math.erfc(mu*k*k/(2*n*math.sqrt(2)))
                                             if method == 'sum' else np.nan)))
    return pd.DataFrame(rows)
