# STAT 236 — Sparse Submatrix Detection and Recovery

Research project by Amrit Chadha and Faris Doughan. The source of truth for reported empirical results is the offline simulation code and the files in `results/`. The website's separate live simulator is a demonstration and **does not replace** these archived runs.

## What was simulated

An `n × n` observed matrix has independent standard Gaussian entries under the null. Under the alternative, `k` rows and `k` columns are selected uniformly without replacement (the sets need not be consecutive), and the entries of their Cartesian product receive a positive mean shift `mu`. A single planted matrix is reused across the competing procedures for paired evaluation.

### Detection: false-positive rate and power

- Global sum, normalized as `X.sum() / n`.
- Maximum entry (`X.max()`).
- Column pooling, the maximum column sum divided by `sqrt(n)`; included as a secondary comparator rather than a central square-block procedure.
- Exact exhaustive scan for small matrices: largest normalized `k × k` block sum over **every** candidate row and column set. Enumeration of all row sets with top-`k` column selection is mathematically identical to full row-and-column enumeration, and is verified against brute force.
- Sum-or-scan combined test, calibrated by Bonferroni splitting `alpha=0.05` into two level-0.025 component tests. The combined test is only computed in scenarios for which the exact scan is feasible.

All rejection thresholds are calibrated on **1,500 independent null matrices** per configuration and are fixed before evaluating alternatives. We use the conformal-style empirical upper order statistic `ceil((m+1)(1-alpha))` (one-indexed) for each component; this controls the expected marginal null error under independent, exchangeable samples. For each `mu`, **600 new null samples** and **600 independently planted alternatives** are tested. The Wilson score 95% confidence intervals describe the conditional binomial evaluation uncertainty and **do not include variability in the calibration cutoff**. Small fluctuations around 5% false positives are expected. Combined-test calibration uses the 2.5% component quantiles.

Detection cases: `D1 (n=16,k=4)`, `D2 (n=12,k=5)`, `D3 (n=64,k=16)`, and `D4 (n=64,k=4)`. The expensive exact scan appears only in D1 and D2. **No simulated large-matrix scan is approximated and mislabeled “exact.”**

### Recovery: exact index-set identification

We record a correct recovery only if the estimated set of `k` affected rows *and* the estimated set of `k` affected columns both match the planted sets. Reported 95% Wilson intervals come from **600 independent planted matrices** per `mu`. We also record the average of row and column Jaccard overlap. The recovery procedures are:

- Exhaustive block maximization, retaining the maximizing row and column sets. Evaluated only in `R1 (16,4)` and `R2 (12,5)`.
- Spectral projection. The leading singular vectors are jointly sign-oriented using the positive-mean model, and the top-`k` projected row and column scores are retained. Evaluated in all four recovery configurations.
- Fixed soft-thresholded spectral projection for the small-block case `R4 (64,4)`. Entries are transformed by `sign(x) max(abs(x)-t,0)` using the prespecified `t = sqrt(2 log(n/k²))`, then spectral projection is applied. The threshold is not fit on the outcomes or tuned on the alternative data.

`R3` is `(64,16)`. Our finite-sample implementation of the spectral method uses the known `k` to select the top-`k` scores. These experiments are finite-sample illustrations, not a proof of an asymptotic regime or a computational lower bound. The theoretical papers treat broader conditions that should not be claimed to hold unconditionally for this exact simulation.

## Reproduce from code

Requires Python 3.10 or newer, NumPy, and Matplotlib. Use a single BLAS thread for more stable wall-time comparisons. All sampling uses deterministic NumPy `SeedSequence` streams keyed by scenario, so `--only D2` matches D2 in a full run (except for measured runtimes).

```bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python tests/test_statistics.py
OPENBLAS_NUM_THREADS=1 python experiments/run_experiments.py --mode full --seed 2362026
OPENBLAS_NUM_THREADS=1 python experiments/build_figures.py
OPENBLAS_NUM_THREADS=1 python experiments/build_site.py
python tests/test_site.py  # Optional: requires Playwright and a local Chromium binary
```

`--mode quick` uses fewer replications for development only; quick outputs should not be substituted for the paper results. `--only D1` (etc.) limits a run to one scenario, with matching seeded statistical outcomes.

## What the files mean

- `results/trials.csv`: 50,400 raw evaluation trials, with test statistics/decisions or recovery success/partial overlap.
- `results/calibration.csv`: thresholds from independent null calibration (only thresholds, not all 1,500 raw null statistics; the generator reruns those exactly).
- `results/summary.csv`: per-method estimated power, held-out false-positive rate, recovery success, partial overlap, 95% Wilson intervals, and recorded per-evaluation runtime (machine dependent).
- `results/metadata.json`: experimental configurations, random seed scheme, number of trials, and environment versions.
- `figures/detection_power.{pdf,png,svg}` and `figures/exact_recovery.{pdf,png,svg}`: publication-ready figures generated from the summary CSV, not hand drawn.
- `paper_simulation_section.tex`: ready-to-include LaTeX text, a results table and figure references generated from `summary.csv`.
- `STAT236_Project_With_Simulations.tex` and `.pdf`: a full report draft with that section inserted.
- `site/index.html`: complete standalone interactive website with **embedded fixed results** and a separate visitor-run Monte Carlo module.
- `tests/test_statistics.py`: brute-force verification of the optimized exact scan, spectral easy-case verification, and calibration checks.

## Interpretation and limits

At `(n,k)=(16,4)` the exact scan outperforms global aggregation for many moderate shifts. For `(64,16)` global aggregation can be effective on weak signals, whereas in the smaller planted-block `(64,4)` case an extreme-value test may have better power than the global sum for large entrywise mean shifts. At `(64,4)`, a thresholded spectral procedure improves exact recovery over the unthresholded baseline at intermediate shifts. These are **finite-sample observations** subject to Monte Carlo and calibration uncertainty; they do not establish statistical optimality or computational hardness for the tested dimensions.

## Source literature

Butucea & Ingster (2013); Ma & Wu (2015); Butucea, Ingster & Suslina (2015); Cai, Liang & Rakhlin (2017). See the site's references and the project report for the detailed hypotheses and guarantees.

## GitHub raw trial data

GitHub stores the raw trial records in `results/trial_parts/trials_01.csv`, etc., each with the same header. These partitions concatenate row-for-row to the raw `results/trials.csv` in the downloadable archive. To reconstruct the single CSV from a GitHub checkout, run `python experiments/restore_trials.py`. The replayed results can then be regenerated with `experiments/run_experiments.py` and compared with `results/summary.csv`.
## Colab notebook: matrix-size and sparsity-exponent experiments

[Open the configurable STAT 236 experiment notebook in Google Colab](https://colab.research.google.com/github/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/blob/main/notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb).

The notebook contains the original simulation methods and correctness checks, loads all **50,400 archived paper trials** from the ten GitHub CSV partitions, and adds experiments over matrix size `n`, sparsity exponent `gamma` (with `k=floor(n**gamma)`, clipped for small matrices), and signal strength `mu`. It produces detection power, Type I error, exact recovery, runtime figures, confidence intervals, calibration cutoffs, and per-trial CSV exports. It compares requested gamma with the **effective** exponent after integer rounding and deduplicates repeated (n,k) pairs. Exact scan is only enabled for small enough combinatorial searches.

The [single official notebook](notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb) includes executed plots and tables and remains fully rerunnable. The [new sweep summary](results/gamma_sweep/summary.csv), [calibration cutoffs](results/gamma_sweep/calibration.csv), [gamma mapping](results/gamma_sweep/gamma_mapping.csv), and [metadata](results/gamma_sweep/metadata.json) are committed. The new exploratory run comprised **31 distinct configurations** and **93,000 trials**, with 600 independent calibration nulls per configuration and 200 fresh null / 200 planted evaluations per signal level, plus 200 recovery evaluations per signal level. It is a separate exploratory study, **not a replacement for the fixed paper run**. A complete archive containing the newly generated raw trial rows and graphics is also available as a download in the corresponding ChatGPT conversation; the notebook regenerates them with a fixed seed and supports per-configuration checkpoints.

Google Colab's default filesystem is ephemeral. To resume a multi-session study, copy its outputs to persistent Google Drive storage or set the notebook's output directory to a mounted Drive folder before running.
