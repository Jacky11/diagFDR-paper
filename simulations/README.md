# Ground-truth simulations for diagFDR

R scripts that reproduce the simulation results of
*diagFDR: Assessing False Discovery Rate Claims in Proteomics with Scope,
Calibration and Stability Diagnostics* (main text, "Simulations against known ground truth";
Supplementary Note S6, Tables S4, S8, S9, S10 and Figure S1).

The simulations answer one question: when a target-decoy list is built under
known conditions, do the quantities the list reports, and the diagnostic flags
computed from it, reveal what its true error rate is?

## Requirements

R 4.0 or later. No package outside base R is needed (`parallel` ships with R).

## Running

```
Rscript run_all.R                  # everything, about 2 hours on one core
SIM_CORES=4 Rscript run_all.R      # same results, four cores
SIM_QUICK=1 Rscript run_all.R      # 1-minute check of the pipeline (noisy)
```

Each script can also be run on its own from this folder. Results are written to
`out/tables` and `out/figures`. Every replicate has its own seed, so the output
is identical whether the scripts run serially or in parallel.

## The model

A simulated list has `n = 200,000` units (for example precursors).

* A unit is **null** with probability `pi0 = 0.667`: its target is an incorrect
  match.
* Each unit has a target score and a decoy score. Null scores come from one
  distribution, standardised to mean 0 and variance 1: Gaussian, Gumbel, or
  Student t with 4 degrees of freedom.
* The target of a non-null unit scores `delta = 4.65` higher.
* Unless stated otherwise the two scores **compete**: only the higher one is
  kept, with its label (target or decoy). For a null unit the two labels are
  then equally likely, which is the "equal chance" assumption.
* q-values are computed by target-decoy counting: entries are ranked by score,
  the estimated FDR at each rank is (decoys above) / (targets above), and the
  q-value is the running minimum from the bottom of the list.

Because the null status of every unit is known, the **realised FDP** (incorrect
accepted targets / accepted targets) is known exactly for every list.

`pi0` and `delta` are design values (Supplementary Note S6.2). Both were first
derived from the DIA-NN export, but both derivations are set by bounds:
`pi0 = 0.667` is the export ceiling (`2 q_max / (1 + q_max)` with
`q_max = 0.5`), and `delta = 4.65 = 2 * qnorm(0.99)` comes from a sensitivity
fraction above one (1.656) clipped to 0.99, i.e. 99% of correct targets exceed
the 1% upper quantile of the null. They give `E[D] = 659` at 1% under S0, of
the order of the 832 decoys of the DIA-NN global list. `02_kappa_robustness.R`
varies both.

## Scenarios

| Scenario | What changes | What it tests |
|---|---|---|
| S0 | nothing | negative control |
| S1 | decoy scores shifted by -0.35 | decoys too easy to beat |
| S2 | decoy scores multiplied by 0.85 | decoys with a depleted upper tail |
| S3 | targets and decoys listed separately (no competition); null fraction among targets lowered to 0.45 | a second search pass on a library reduced to first-pass identifications |
| S4a | minimum q-value over 6 independent runs; unit labelled decoy if a decoy won in any run | cross-run aggregation |
| S4b | same runs; unit takes the label of the run that gave the minimum q-value | cross-run aggregation |
| S5 | scores correlated within blocks of 10 units (correlation 0.2, 0.5, 0.8); Gaussian only | dependence between units, equal chance intact |

In S1 and S2 the decoy model is wrong. In S3, S4 and S5 the decoys are correct:
what changes is how the list is assembled (S3, S4) or the dependence between
units (S5).

## Quantities recorded

For each replicate and each threshold `alpha` in {0.05, 0.01, 0.002, 0.0005}:

| Name | Meaning |
|---|---|
| `T`, `D` | accepted targets, accepted decoys |
| `V` | accepted targets that are null |
| `fdp` | realised FDP, `V / T` |
| `decoy_ratio` | `D / T`, what the list reports |
| `pi_D` | decoy fraction among entries with q-value in [0.2, 0.5] (equal-chance statistic) |
| `ipe_tdc`, `ipe_pooled` | internal PEP error against `D_b/T_b` and `D_b/(T_b+D_b)`; computed where the true PEP has a closed form (competed Gaussian lists) |

## Flags

| Flag | Fires when |
|---|---|
| `D` | fewer than 10 accepted decoys |
| `CV` | decoy-support index `D^(-1/2)` above 0.30 |
| `IPE` | internal PEP error above 0.05 |
| `EQ` | `pi_D` more than 0.15 away from 0.5 (at least 100 entries in the band) |

A replicate is **bad** when its realised FDP exceeds 1.5 x `alpha`. The
**sensitivity** of a flag is the share of bad replicates in which it fires.

`kappa` is the coefficient of variation of the realised FDP across replicates,
divided by the mean of `D^(-1/2)`. Proposition 1 of the paper predicts
`kappa = sqrt(2)` under scenario S0.

## Files

| File | Content | Output |
|---|---|---|
| `parameters.R` | parameters shared by all scripts | |
| `R/sim_engine.R` | score generation, q-values, scenarios, flags, summaries | |
| `01_reference_scenarios.R` | S0 to S5 at the reference parameters, 1000 replicates | `reference_summary.csv`, `table_flags.tex` (Table S4) |
| `02_kappa_robustness.R` | `kappa` under S0 for five (`pi0`, `delta`) pairs, 500 replicates | `kappa_robustness.csv`, `table_kappa.tex` (Table S8) |
| `03_misspecification_grid.R` | S1 and S2 at five magnitudes each, 200 replicates | `misspecification_grid.csv`, `table_grid.tex` (Table S9) |
| `04_library_enrichment.R` | S3 for seven library compositions, 200 replicates | `library_enrichment.csv`, `table_enrichment.tex` (Table S10) |
| `05_discrimination_vs_support.R` | accepted decoys against separation under two operating rules, 200 replicates | `discrimination_vs_support.csv`, `figS1_signflip.pdf` (Figure S1) |
| `06_proposition1_check.R` | S0, 1000 replicates per null family (Gaussian, Gumbel, Student-t4): Var(V) against 2E[D] + Var(D), kappa, the conditional law of V given D; share of accepted decoys coming from units carrying a correct target (assumption iii) | `proposition1_check.csv`, `proposition1_decoy_share.csv` (Note S1) |
| `run_all.R` | runs the six scripts in order | |
