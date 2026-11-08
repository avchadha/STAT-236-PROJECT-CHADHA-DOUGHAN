"""Finite-sample checks of the writeup's signal scales and rejection thresholds."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

from notebook_studies import reference_scales, theory_errors, analytic_detection_check

ROOT = Path(__file__).resolve().parents[1]


def normal_cdf(x):
    return .5 * math.erfc(-x / math.sqrt(2))


def binomial_upper_tail(count, total, probability):
    if count == 0 or probability == 1:
        return 1.
    if probability == 0:
        return 0.
    logs = [math.lgamma(total+1)-math.lgamma(j+1)-math.lgamma(total-j+1)
            + j*math.log(probability) + (total-j)*math.log1p(-probability)
            for j in range(count, total+1)]
    largest = max(logs)
    return min(1., math.exp(largest) * math.fsum(math.exp(v-largest) for v in logs))


def threshold_consistency(trials):
    """One-sided exact binomial checks against the paper's finite-sample bounds."""
    comparisons = []
    errors = theory_errors(trials)
    for (n, k, mu), group in trials[trials.task == 'detection'].groupby(['n','k','mu']):
        log_candidates = 2 * (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))
        scan_cutoff = math.sqrt(2*log_candidates)
        maximum_cutoff = math.sqrt(4.5*math.log(n))
        for condition in ['null','signal']:
            raw = group[group.condition == condition]
            total = len(raw)
            for method in ['sum','scan','maximum']:
                field = 'stat_'+method
                if field not in raw or raw[field].isna().all():
                    continue
                cutoff = {'sum':mu*k*k/(2*n),'scan':scan_cutoff,'maximum':maximum_cutoff}[method]
                count = int((raw[field] > cutoff).sum() if condition == 'null'
                            else (raw[field] <= cutoff).sum())
                if method == 'sum':
                    bound = normal_cdf(-mu*k*k/(2*n))
                    reference = 'exact Gaussian sum error'
                elif method == 'scan':
                    bound = min(1., math.exp(log_candidates)*normal_cdf(-scan_cutoff)) if condition == 'null' \
                            else normal_cdf(scan_cutoff-mu*k)
                    reference = 'scan union bound' if condition == 'null' else 'true-block miss bound'
                else:
                    bound = n**(-.25) if condition == 'null' else math.exp(-.5*max(mu-maximum_cutoff,0.)**2)
                    reference = 'maximum union bound' if condition == 'null' else 'affected-entry miss bound'
                comparisons.append(dict(n=n,k=k,mu=mu,method=method,condition=condition,
                    cutoff=cutoff,count=count,replicates=total,observed=count/total,
                    probability_bound=bound,reference=reference,
                    p_value=binomial_upper_tail(count,total,bound)))
    checks = pd.DataFrame(comparisons)
    checks['familywise_alpha'] = .01
    checks['bonferroni_cutoff'] = .01/len(checks)
    checks['consistent'] = checks.p_value >= checks.bonferroni_cutoff
    assert checks.consistent.all(), checks.loc[~checks.consistent].to_dict('records')
    return errors, checks


def exact_gaussian_consistency(analytic, replicates=600, family_size=None):
    checks=analytic.copy()
    probabilities=[]
    for r in checks.itertuples():
        count=round(r.observed*replicates)
        p=r.exact_probability
        probabilities.append(min(1.,2*min(binomial_upper_tail(count,replicates,p),
                                          binomial_upper_tail(replicates-count,replicates,1-p))))
    checks['p_value']=probabilities
    checks['bonferroni_cutoff']=.01/(family_size or len(checks))
    checks['consistent']=checks.p_value>=checks.bonferroni_cutoff
    assert checks.consistent.all(),checks.loc[~checks.consistent].to_dict('records')
    return checks


def verify_source_thresholds():
    """Check that the experiment references agree with both displayed sources."""
    paper = (ROOT/'STAT236_Project_With_Simulations.tex').read_text(encoding='utf-8')
    website = BeautifulSoup((ROOT/'site/base-v4.html').read_text(encoding='utf-8'),'html.parser')
    web_math = ' '.join(t['data-tex'] for t in website.select('[data-tex]')).replace('\\\\','\\')
    compact_paper = ''.join(paper.split())
    compact_web = ''.join(web_math.split())
    checks = {
        'sum_shift': '\\mu k^2/n',
        'scan_ratio': '\\mu^2k/[4\\log(n/k)]',
        'maximum_scale': '2\\sqrt{\\log n}',
        'recovery_local_rows': '\\sqrt{2\\log k}',
        'recovery_local_noise': '\\sqrt{2\\log(n-k)}',
        'spectral_large_block': '\\sqrt n/k',
        'spectral_small_block': '\\log(n/k^2)',
    }
    rows = []
    for name, expression in checks.items():
        token = ''.join(expression.split())
        in_paper = token in compact_paper
        in_website = token in compact_web
        assert in_paper and in_website, (name,in_paper,in_website)
        rows.append(dict(reference=name,paper=True,website=True))
    # Independent numerical specializations of the displayed square-block formulas.
    a = reference_scales(64,4)
    assert a['sum'] == 4
    assert math.isclose(a['scan'],math.sqrt(math.log(16)))
    assert math.isclose(a['maximum'],2*math.sqrt(math.log(64)))
    assert math.isclose(a['spectral'],2+math.sqrt(math.log(64)/4))
    assert math.isclose(a['soft_spectral'],math.sqrt(math.log(4)))
    assert math.isclose(a['exact_recovery'],(math.sqrt(2*math.log(4))+math.sqrt(2*math.log(60)))/2)
    return rows


def summarize_recovery(summary, plan):
    rows = []
    for cfg in plan:
        n,k = cfg['n'],cfg['k']
        for reference, method in [('exact_recovery','exhaustive'),('spectral','spectral'),
                                  ('soft_spectral','soft_spectral')]:
            for signal in cfg['signals']:
                if signal['reference'] != reference:
                    continue
                ss=summary[(summary.n==n)&(summary.k==k)&(summary.method==method)&
                           (summary.metric=='exact_recovery')&np.isclose(summary.mu,signal['mu'],atol=1e-12,rtol=0)]
                for r in ss.itertuples():
                    rows.append(dict(n=n,k=k,method=method,reference=reference,
                        multiplier=signal['parameter'],mu=r.mu,estimate=r.estimate,
                        ci_low=r.ci_low,ci_high=r.ci_high,
                        regime=('small block' if n>k*k else 'large block'),
                        theorem_applies=('spectral'!=reference or k*k>n)))
    return pd.DataFrame(rows)
