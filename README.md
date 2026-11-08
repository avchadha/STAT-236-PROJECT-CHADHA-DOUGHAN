# STAT 236 - Sparse Submatrix Detection and Recovery

Amrit Chadha and Faris Doughan.

## Notebook

The [executed notebook](notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb) is a short, self-contained experiment: create labeled matrices with **seed 236**, define the methods, and evaluate their errors. No project imports, downloads, execution modes, or checkpoint setup are needed. The [HTML reading copy](notebooks/STAT236_Sparse_Submatrix_Experiments.html) includes the same code and saved outputs.

- Cross `n = 12,16` with `gamma = .4,.75`, giving `(n,k) = (12,2),(12,6),(16,3),(16,8)`.
- Use 600 null matrices and 600 planted matrices at each `mu = .25,.5,.75,1,1.5,2,3,4`. Each dataset row records the matrix, signal label, signal strength, and planted indices.
- Compare sum, exact scan, maximum entry, and Bonferroni sum+scan detection. Sum and maximum use exact Gaussian cutoffs; the scan uses 1,500 separate calibration matrices per configuration.
- Compare exhaustive, spectral, and small-block soft spectral recovery on the same planted matrices. Exact recovery requires both index sets to match.
- Produce two ordinary error-versus-signal figures, then one growing-size figure with `n = 16,32,64,128,256`, `k = floor(n**.75)`, and `mu = n**(-beta)` for `beta = 0,.25,.75`. This last experiment evaluates sum detection and spectral recovery, with detection level `1/n`.

Detection error is half the sum of the false-positive and false-negative rates; the paper's total testing risk is twice this value. Recovery error is the fraction of planted matrices whose row or column set is wrong. These small experiments illustrate selected finite-sample and asymptotic regimes; they do not verify every theorem in the source papers.

```bash
python -m pip install -r requirements-notebook.txt
python experiments/build_notebook.py
python tests/test_notebook_studies.py
python experiments/execute_notebook.py
python experiments/export_notebook.py
```

The notebook has been simplified independently of the larger published study below. Its seed, design, and results differ from the paper and website's existing results.

## Published paper and website studies

The paper and website use [`results/final/`](results/final/), generated with seed **2362026**. Every detection configuration uses **1,500 independent null calibration matrices** and **600 fresh null plus 600 planted matrices per signal**. Recovery uses **600 additional planted matrices per signal**.

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

All 577,800 published raw evaluation records are checked against every expected configuration, signal, procedure, replication count, and reported success count. The simplified notebook's exact scan and maximizing indices are checked against independent brute-force enumeration; both spectral variants are compared with the published implementation. Matrix generation is checked against the planted truth.

The simplified notebook executes without errors and contains three figures. Its 21,600 labeled matrices have 600 trials in every condition; all 72 sum and maximum-entry rejection rates agree with exact Gaussian probabilities at a 1% familywise testing level.

The experiment's reference scales are checked against the paper and website. The 1,476 sum and maximum-entry rejection rates are compared with exact Gaussian probabilities using exact binomial diagnostics and a 1% familywise threshold across that family. The threshold study passes 350 one-sided checks against the sum, scan, and maximum-entry error calculations at a separate 1% familywise level. Recovery curves report finite-sample success below, at, and above reference scales. Unknown constants, logarithmic factors, asymptotic impossibility, and clique-based computational lower bounds are not established by these simulations.

## Reproduce

Requires Python 3.10+; compiling the report also requires LaTeX. Colab already supplies the notebook's main scientific libraries.

```bash
python -m pip install -r requirements-notebook.txt
python tests/test_statistics.py
python experiments/run_final_study.py
python experiments/build_final_report.py
pdflatex -interaction=nonstopmode -halt-on-error STAT236_Project_With_Simulations.tex
pdflatex -interaction=nonstopmode -halt-on-error STAT236_Project_With_Simulations.tex
python experiments/build_site.py
python tests/test_final_study.py
python tests/test_site_layout.py
```

The three published designs use independent deterministic `SeedSequence([seed, task_stream_key, n, k])` streams with master seed 2362026. Completed configurations are verified and resumed only for matching settings. `run_final_study.py` reuses those verified checkpoints and exports the published dataset. `build_site.py` writes identical generated pages at the root and `site/index.html`, embedding the same verified summaries and current report. The website's visitor-run simulator is a separate demonstration. Detection and recovery procedures follow separate connected branches; the recovery boundary precedes the comparison table. Layout tests preserve the original text and formulas of all eight procedure descriptions.

The simplified notebook keeps its labeled data and per-matrix decisions in memory and saves its figures in the notebook itself. Numerical-library tie handling can vary across machines.

## Earlier results

The original D1-D4/R1-R4 studies and the earlier 200-replication exploratory sweeps are superseded in the current notebook, paper, and website. Their complete files remain available at [the previous repository version](https://github.com/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/tree/295c1702a6f52350dc3a43dafd8fa97b7023f95e). They are historical results, not evidence of contaminated data. `run_experiments.py` retains the original scenario runner for historical reproduction; use `run_final_study.py` for the active study.

## Source literature

Butucea and Ingster (2013); Ma and Wu (2015); Butucea, Ingster and Suslina (2015); Cai, Liang and Rakhlin (2017). See the report for the hypotheses, constants, and scope of each guarantee.
