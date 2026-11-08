# STAT 236 - Sparse Submatrix Detection and Recovery

Amrit Chadha and Faris Doughan. The notebook, paper, and website now use one active dataset: [`results/final/`](results/final/). Every final detection configuration uses **1,500 independent null calibration matrices** and **600 fresh null plus 600 planted matrices per signal**. Recovery uses **600 additional planted matrices per signal**.

## Final studies

| Study | Grid | Evaluation trials |
|---|---|---:|
| Matrix size and sparsity | `n = 6,8,12,16,24,32,64`; requested `gamma = .25,.40,.55,.70,.85`; five scan-scaled signals; 31 distinct `(n,k)` pairs | 279,000 |
| Signal reference levels | `(12,2)`, `(12,5)`, `(64,4)`, `(64,16)`; half, once, and twice every applicable signal reference, plus zero | 126,000 |
| Power laws | `n = 12,24,64`; requested `gamma = .4,.5,2/3,.8`; beta values at and around the exponent boundaries, plus a zero-signal check; detection only | 172,800 |
| Total | 47 study/configuration combinations; 70,500 calibration matrices counted separately | **577,800** |

The block size is `floor(n**gamma)`, clipped to `[2,n-1]`. Near-integer floating-point values are handled explicitly, so `64**(2/3)` gives 16. Actual block sizes and effective exponents are recorded. All chosen signal values run automatically. Exact scanning and block maximization run only when `comb(n,k) <= 2500`; omitted methods are uncomputed, and no heuristic is labeled exact.

The detectors are global sum, maximum entry, column pooling, exact scan, and Bonferroni-calibrated sum-or-scan. Recovery uses exact block maximization, spectral projection, and fixed soft-thresholded spectral projection where `k*k < n`. The implementation uses full SVD and known `k`; the faster power-iteration cost discussed in the theory is not its measured runtime. Exact success requires both planted index sets to match. Wilson 95% intervals describe evaluation uncertainty conditional on the calibrated cutoff.

## Read the results

- [Website](https://avchadha.github.io/STAT-236-PROJECT-CHADHA-DOUGHAN/) and [website source](index.html).
- [Executed notebook](notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb), [reading copy](notebooks/STAT236_Sparse_Submatrix_Experiments.html), or [open in Colab](https://colab.research.google.com/github/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/blob/main/notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb).
- [Current report](STAT236_Project_With_Simulations.pdf) and [report source](STAT236_Project_With_Simulations.tex).
- [Final dataset manifest](results/final/manifest.json), including the complete design and verification counts.

Each study folder contains `summary.csv`, `calibration.csv`, `configurations.csv`, `design.json`, numerical environment details, and compressed raw trials in `trials.csv.gz`. The threshold study also includes theoretical-threshold errors, bound checks, and recovery reference comparisons. `paper_selected_results.csv` records every number selected for the paper's results table. Read compressed trials with `pandas.read_csv(path, keep_default_na=False, na_values=[''])` to preserve the condition label `null`.

## What was verified

All 577,800 raw evaluation records are checked against every expected configuration, signal, procedure, replication count, and reported success count. Statistical function definitions embedded in the notebook are checked against the source implementation. Exact scan and maximizing indices are checked against independent brute-force enumeration; both spectral variants are checked on strong planted signals.

The experiment's reference scales are checked against the paper and website. The 1,476 sum and maximum-entry rejection rates are compared with exact Gaussian probabilities using exact binomial diagnostics and a 1% familywise threshold across that family. The threshold study passes 350 one-sided checks against the sum, scan, and maximum-entry error calculations at a separate 1% familywise level. Recovery curves report finite-sample success below, at, and above reference scales. Unknown constants, logarithmic factors, asymptotic impossibility, and clique-based computational lower bounds are not established by these simulations.

## Reproduce

Requires Python 3.10+; compiling the report also requires LaTeX. Colab already supplies the notebook's main scientific libraries.

```bash
python -m pip install -r requirements-notebook.txt
python tests/test_statistics.py
python experiments/build_notebook.py
python tests/test_notebook_studies.py
python experiments/execute_notebook.py
python experiments/run_final_study.py
python experiments/build_final_report.py
pdflatex -interaction=nonstopmode -halt-on-error STAT236_Project_With_Simulations.tex
pdflatex -interaction=nonstopmode -halt-on-error STAT236_Project_With_Simulations.tex
python experiments/build_site.py
python experiments/export_notebook.py
python tests/test_final_study.py
python tests/test_site_layout.py
```

`build_notebook.py` embeds the statistical functions directly from `run_experiments.py`; running the notebook needs no custom project imports or archived downloads. The three designs use independent deterministic `SeedSequence([seed, task_stream_key, n, k])` streams with master seed 2362026. Completed configurations are verified and resumed only for matching settings. `run_final_study.py` reuses those verified checkpoints and exports the active dataset. `build_site.py` writes identical generated pages at the root and `site/index.html`, embedding the same verified summaries and current report. The website's visitor-run simulator is a separate demonstration. Detection and recovery procedures follow separate connected branches; the recovery boundary precedes the comparison table. Layout tests preserve the original text and formulas of all eight procedure descriptions.

The notebook saves checkpoints in `outputs_notebook/` and writes `STAT236_notebook_results.zip`. Colab storage is temporary: download the ZIP to keep a run. Numerical-library tie handling and measured times can vary across machines; recorded environment details accompany the results.

## Earlier results

The original D1-D4/R1-R4 studies and the earlier 200-replication exploratory sweeps are superseded in the current notebook, paper, and website. Their complete files remain available at [the previous repository version](https://github.com/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/tree/295c1702a6f52350dc3a43dafd8fa97b7023f95e). They are historical results, not evidence of contaminated data. `run_experiments.py` retains the original scenario runner for historical reproduction; use `run_final_study.py` for the active study.

## Source literature

Butucea and Ingster (2013); Ma and Wu (2015); Butucea, Ingster and Suslina (2015); Cai, Liang and Rakhlin (2017). See the report for the hypotheses, constants, and scope of each guarantee.
