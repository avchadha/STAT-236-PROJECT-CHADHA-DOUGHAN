# Sparse Submatrix Detection

Amrit Chadha and Faris Doughan · STAT 236 · Fall 2026

[Project website](https://avchadha.github.io/STAT-236-PROJECT-CHADHA-DOUGHAN/)

The website introduces detection and exact recovery of a square elevated-mean block in Gaussian noise. It contains a single-matrix explorer, a summary of the methods, and the two current experiment figures.

The report and notebook will be added when their final versions are ready.

## Run the website

There is no build step. From this directory, run

```sh
python -m http.server 8000
```

Open http://localhost:8000. A local server is needed for the scan worker. MathJax loads from jsDelivr to typeset the LaTeX equations, so math rendering requires an internet connection. The matrix explorer runs independently of MathJax.

Changing the mean shift keeps the noise fixed. Sampling gives fresh noise and a fresh block location. Revealing or grouping the planted rows and columns changes only the display. The test statistics never receive the planted indices.

Sum and max use exact level-.05 null cutoffs. The browser scan uses the conservative cutoff `sqrt(2 log(M / alpha))`, where `M = choose(n,k)^2`. The combined decision uses .025 per component. Exact scans enumerate row sets and select the best columns for each. Searches beyond the browser work limit are skipped, and a combined negative decision is withheld if the scan is unavailable.

## Checks

```sh
node tests/test_explorer.js
```

The checks compare the scan with full enumeration on small matrices and exercise sampling, parameter changes, reveal controls, worker results, and skipped scans. They do not require a browser installation.

## Files

- `index.html` contains the text, LaTeX equations, and page structure
- `assets/` contains the stylesheet, matrix explorer, and statistical functions
- `figures/` contains the detection and recovery figures
- `tests/` contains the numerical and controller checks
