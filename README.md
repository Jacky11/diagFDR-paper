# diagFDR-paper

Code and tables behind *diagFDR: Assessing False Discovery Rate Claims in
Proteomics with Scope, Calibration and Stability Diagnostics*.

The `diagFDR` R package itself lives elsewhere: source
<https://github.com/Jacky11/diagFDR> (version 0.1.2, used for the paper) and
CRAN <https://cran.r-project.org/package=diagFDR>. This repository holds the
analyses. Install the package with
`remotes::install_github("Jacky11/diagFDR")` until 0.1.2 is on CRAN.

| Folder | Content | Paper |
|---|---|---|
| `simulations/` | ground-truth simulations (base R only) | Section 4.1, Proposition 1, Notes S1 and S6, Tables S4, S8-S10, Figure S1 |
| `real_data/` | DIA-NN and MS2Rescore analyses | Sections 4.2-4.5, Table 3, Figures 3-5, Notes S4, S5, S8, S9, Table S6 |
| `figures/` | Figure 2 (workflow and report) and the relabelling of Figures 3-5 | Figures 2-5 |
| `audit/` | audit of FDR reporting in 63 DIA articles | Section 2.2, Note S7, Tables S3 and S5, Data S1 |

## Data

Raw files, DIA-NN reports and mokapot outputs:
<https://doi.org/10.5281/zenodo.18621806>. The DIA benchmark is
ProteomeXchange PXD028735. The simulations need no data.

## Quick start

```
cd simulations
SIM_QUICK=1 Rscript run_all.R        # one-minute check
Rscript run_all.R                    # full run, about 2 hours on one core
```

```
cd real_data
export DIAGFDR_DATA=/path/to/zenodo/files
Rscript run_all.R
```

Each folder has its own README.
