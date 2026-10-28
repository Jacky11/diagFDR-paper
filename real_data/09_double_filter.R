## =====================================================================
##  09 --- The DIA-NN "double filter" as a fourth scope choice (NEW)
##
##  Reviewers may object that DIA-NN users do not build experiment-wide
##  lists from min_r(Q.Value) alone but apply BOTH Q.Value <= alpha and
##  Global.Q.Value <= alpha per row (the default of the DIA-NN GUI and of
##  --matrix-qvalue). This script adds that construction to the scope
##  audit: rows are kept when both columns pass, then collapsed to one
##  entry per precursor.
##
##  RESULT (verified 2026-10-25 with diagFDR 0.2.0): on all three searches
##  the double-filtered list is IDENTICAL to the Global.Q.Value list
##  (basic T = 83,149 D = 832; entrapment T = 81,643 D = 816 E = 177;
##  Jaccard = 1), because a precursor whose experiment-wide q-value passes
##  alpha always has at least one run passing it. The default is safe; the
##  min_r construction is what filtering on Q.Value ALONE produces.
##
##  Also reports the run-wise-only practice found in the audit (Note S7):
##  of the 36 DIA-NN studies, 26 state no scope (Yes or Partial: 10), and 2 of the 10 that do
##  filter on Q.Value / Lib.Q.Value only (audit ids A29, A50).
##
##  Output: out/tables/double_filter.csv ; numbers registered as dfT, dfD,
##  dfRatio, dfEntrap, dfEntrapCorr, dfJaccardGlobal.
## =====================================================================
source("00_config.R")
load_numbers()

double_filter_universe <- function(rep, alpha) {
  keep <- rep[[COL$q_run]] <= alpha & rep[[COL$q_global]] <= alpha
  sub <- rep[keep, , drop = FALSE]
  ## the retained q is Global.Q.Value (the experiment-wide column); the
  ## row filter on Q.Value is what changes the list
  diagFDR::collapse_universe(sub, COL$precursor, COL$q_global, rule = "winning_row",
                             decoy_col = COL$decoy, pep_col = COL$pep,
                             q_max_export = alpha, unit = "precursor",
                             scope = "double_filter", q_source = "Q.Value & Global.Q.Value",
                             tool = "DIA-NN", two_pass = TRUE)
}

res <- list()
for (nm in c("basic", "semi", "entrap")) {
  path <- switch(nm, basic = DIANN_BASIC, semi = DIANN_SEMI, entrap = DIANN_ENTRAP)
  if (!file.exists(path)) next
  rep <- sanitise_report(diagFDR::read_diann_parquet(path), COL$decoy, nm)
  u   <- double_filter_universe(rep, ALPHA_MAIN)
  g   <- diagFDR::diann_global_precursor(rep, q_col = COL$q_global, q_max_export = QMAX$q_global, two_pass = TRUE)
  h   <- headline_of(u, ALPHA_MAIN)
  row <- data.frame(search = nm, T_alpha = h$T_alpha, D_alpha = h$D_alpha,
                    decoy_ratio = h$decoy_ratio, CV_benchmark = h$CV_benchmark,
                    jaccard_vs_global = diagFDR:::jaccard_vec(diagFDR::dfdr_accepted(u, ALPHA_MAIN),
                                                               diagFDR::dfdr_accepted(g, ALPHA_MAIN)))
  if (nm == "entrap") {
    ACC <- diagFDR::read_fasta_accessions(FASTA_ENTRAP)
    ent <- diagFDR::diann_entrapment_ids(rep, ACC, id_col = COL$precursor, protein_col = COL$protein, rule = ENTRAP_RULE)
    e <- diagFDR::dfdr_entrapment(u, ent, alpha = ALPHA_MAIN, r = 0.432)
    row$E_alpha <- e$E_alpha; row$fdp_raw <- e$fdp_raw; row$lo <- e$lo; row$hi <- e$hi
    row$fdp_corrected <- e$fdp_corrected
    d <- diagFDR::dfdr_entrapment_discordant(diagFDR::dfdr_accepted(g, ALPHA_MAIN),
                                             diagFDR::dfdr_accepted(u, ALPHA_MAIN), ent)
    print(d)
    reg("dfEntrap", e$fdp_raw, "d3", "Sec 4.2 double-filter raw entrapment FDP")
    reg("dfEntrapCorr", e$fdp_corrected, "d3", "Sec 4.2 double-filter corrected")
  }
  res[[nm]] <- row
}
cols <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
print(res); write_tab(res, "double_filter")
stopifnot(all(abs(res$jaccard_vs_global - 1) < 1e-12))
reg("dfT", res$T_alpha[res$search == "basic"], "int", "Sec 4.2 double filter, basic")
reg("dfD", res$D_alpha[res$search == "basic"], "int", "Sec 4.2 double filter, basic")
reg("dfRatio", res$decoy_ratio[res$search == "basic"], "d3", "Sec 4.2 double filter, basic")
reg("dfJaccardGlobal", res$jaccard_vs_global[res$search == "basic"], "d2", "overlap with global list")
