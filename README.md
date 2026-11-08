# Sparse submatrix detection and recovery

Amrit Chadha and Faris Doughan · STAT 236 · Fall 2026

[Project website](https://avchadha.github.io/STAT-236-PROJECT-CHADHA-DOUGHAN/) · [Report](STAT236_Project_With_Simulations.pdf) · [Open notebook in Colab](https://colab.research.google.com/github/avchadha/STAT-236-PROJECT-CHADHA-DOUGHAN/blob/main/notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb)

We compare detection and exact recovery of a square elevated-mean block in an independent Gaussian noise matrix. The report discusses the sum, scan, and maximum-entry tests, exhaustive recovery, and spectral recovery with and without soft thresholding.

## Website

The website is a static page with no build step or external JavaScript dependencies. From this directory, run

```sh
python -m http.server 8000
```

Open http://localhost:8000. A local server is needed for the scan worker.

The explorer draws one matrix at a time. Changing the mean shift preserves the noise, and the sample button draws fresh noise and a fresh block location. Revealing or grouping the planted rows and columns changes only the display. The statistics never receive the planted indices.

Sum and max use exact level-.05 null cutoffs. The browser scan uses the conservative cutoff `sqrt(2 log(M / alpha))`, where `M = choose(n,k)^2`. This differs from the simulated scan cutoff used in the notebook. The combined decision uses level .025 for each component. Exact scans enumerate row sets and select the best columns for each set. Searches beyond the browser work limit are explicitly skipped, and a combined negative decision is withheld if the scan is unavailable. There is no Monte Carlo runner on the page.

Run the explorer's numerical checks with

```sh
node tests/test_explorer.js
```

## Experiments

The notebook contains two experiments, detection and recovery. Both use 300 trials per condition, with `n=16, k=3,5,8` and `n=512, k=12,42,107`. Detection uses level .05, with .025 per component for the combined test. Scan cutoffs use 1,500 independent noise matrices. Exhaustive scanning and recovery are run only at n=16. All recovery methods receive k.

The seeds are 236 for simulated matrices and 1236 for null calibration, and BLAS uses one thread. The notebook records package versions and implementation settings. Soft thresholding uses `t = sqrt(log(n/k^2))` when `k^2 < n`. The saved outputs and two PNG figures are from the current notebook, with no additional experiments.

To run locally, install `requirements-notebook.txt` and execute the notebook cells in order. Colab supplies these packages. The website links directly to this notebook, so it does not require access to a private Drive file.

## Files

- `index.html` and `assets/` contain the website and matrix explorer
- `notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb` contains the current code and saved results
- `figures/` contains the two figures extracted from those results
- `STAT236_Project_With_Simulations.tex` and `.pdf` contain the report
- `tests/` contains numerical checks for the browser statistics

Compile the report with `pdflatex STAT236_Project_With_Simulations.tex`, twice to resolve references. The indicator is typeset with `\mathbf{1}` so the source does not need the optional dsfont package.

Earlier experiments, duplicate website exports, and study-generation scripts have been removed from the current tree. They remain available in Git history.
