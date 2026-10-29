# Real-data analyses

R scripts for the DIA-NN and MS2Rescore case studies of the paper
(Sections 4.2 to 4.5).

## Requirements

R 4.3 or later with `diagFDR` >= 0.2.1, `arrow`, `dplyr`, `ggplot2`, `cp4p`. About 16 GB of memory
for the import step (three parquet reports of about 100 MB each).

## Data

Download the files of <https://doi.org/10.5281/zenodo.18621806> and point
`DIAGFDR_DATA` at the folder (or edit `DATA` in `00_config.R`). The file
names expected inside that folder are listed at the top of `00_config.R`:

```
diann/Basic_search/report_basic.parquet
diann/Basic_search/report_basic-first-pass.parquet
diann/Semi_Tryp_search/report_semi.parquet
diann/Semi_Tryp_search/report_semi-first-pass.parquet
diann/Basic_search_wih_ArabidopsisT/report_entrapment.parquet
ProteoBenchFASTA_MixedSpecies_HYE.fasta
uniprot-proteome-human-contaminants.fasta
uniprotkb_taxonomy_id_3702_AND_reviewed_2026_02_06.fasta
ms2rescore/maxquant/pyro_mq.ms2rescore.mokapot.psms.txt
ms2rescore/maxquant/pyro_mq.ms2rescore.mokapot.decoy.psms.txt
ms2rescore/maxquant/msms.txt
```

## Running

```
export DIAGFDR_DATA=/path/to/zenodo/files
Rscript run_all.R
```

or run the scripts one at a time, in the order below, from this folder.
Results go to `out/tables` and `out/figures`; `out/universes.rds` and `out/msrescore_universes.rds` feed `../figures/make_figures.R`.

## Scripts

| Script | What it does | Paper |
|---|---|---|
| `00_config.R` | paths, thresholds, column names | |
| `01_preflight.R` | checks the package, the files, the columns, the export ceilings | Note S8 |
| `02_import_diann.R` | imports the three searches and builds every universe and label rule | Methods |
| `03_scope_audit.R` | decoy ratios per universe; four label rules; Jaccard overlap | Table 3, Section 4.2, Figure 3A-B |
| `03b_picked_protein.R` | picked-protein competition; label rules on the entrapment search | Table 3 Panel B, Section 4.2, Note S9 |
| `03c_label_rules_package.R` | the four label rules and picked competition recomputed with the package functions (`dfdr_label_rule_audit()`, `dfdr_picked_protein()`), on the filtered and the unfiltered export; stops if any value differs from the paper | Table 3 Panel B, Note S9 |
| `04_calibration_stability.R` | stability curves, IPE under both comparators, equal chance, sum of PEP | Section 4.4, Figure 4 |
| `05_entrapment.R` | size ratio r, raw and size-corrected entrapment FDP, protein level | Sections 4.2-4.3, Note S4, Figure 3C-D |
| `05b_entrapment_discordant.R` | scope-discordant-set error rate, Equation (5) | Section 4.2 |
| `05c_entrapment_r_sensitivity.R` | homology and sensitivity of the size correction | Note S4 |
| `05d_discordant_size_matched.R` | size-matched control: the min_r list is the Global.Q.Value list at a looser threshold | Section 4.2, Note S4 |
| `05e_entrapment_package.R` | entrapment estimates (raw, combined, sample), size ratio r and discordant-set rates recomputed with `dfdr_entrapment()`, `dfdr_entrapment_size_ratio()`, `dfdr_entrapment_discordant()`; stops if any value differs from the paper | Sections 4.2, Note S4 |
| `05f_entrapment_cluster_bootstrap.R` | protein-group cluster bootstrap for the core / added entrapment rates and the alpha = 1e-4 estimate | Note S4 |
| `06_msrescore.R` | three q-value constructions on the same PSMs | Section 4.5, Figure 5 |
| `06b_msrescore_overlap.R` | before/after counts on the spectrum-level winners (5,610/55 before, 7,070/43 and 7,164/70 after, at 1%) and overlap on a common spectrum key; checks every count of Section 4.5 | Section 4.5, Figure 5E, Note S8.3 |
| `06c_mokapot_universe.R` | the universe on which mokapot's q-values are computed | Section 4.5, Note S8 |
| `08_defaults_sensitivity.R` | window width and PEP binning grids | Note S5, Table S6 |
| `09_double_filter.R` | the double filter (Q.Value and Global.Q.Value, as in DIA-NN's matrices) returns exactly the Global.Q.Value list | Section 4.2, Note S4 |
| `R/dfdr_api.R` | thin adapter over `diagFDR` (>= 0.2.1): bridges column names to the paper's notation; label rules and entrapment estimators now live in the package | |
| `R/utils.R` | Wilson intervals, table and figure helpers, number registry | |

## Notes

* `table2.tex` and the `CV` columns give the decoy-support index
  `D^(-1/2)`. The last column of Table 3 in the paper (scope audit) is the
  FDP-CV benchmark, `sqrt(2)` times that value.
* The scripts write one PDF per panel (`fig2*`, `fig5*`, `fig4_msrescore`,
  `fig6E_overlap`). They are not the composite figures of the paper.
* Each script prints the values expected on the Zenodo data in its header,
  so a rerun can be checked against them.
