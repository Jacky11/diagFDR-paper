## =====================================================================
##  Thin adapter over diagFDR (version 0.2.1 or later)
##
##  REPLACES real_data/R/dfdr_api.R of the first commit. Since diagFDR
##  0.2.0 the label-propagation rules, picked-protein competition and
##  entrapment estimators are part of the package, so this file no longer
##  implements them: it only bridges the package's column names to the
##  paper's notation and keeps the old helper names as aliases so that
##  scripts 01--08 run unchanged.
## =====================================================================

if (!requireNamespace("diagFDR", quietly = TRUE))
  stop("diagFDR is not installed. install.packages('diagFDR')", call. = FALSE)
library(diagFDR)

DFDR_VERSION <- as.character(utils::packageVersion("diagFDR"))
if (utils::packageVersion("diagFDR") < "0.2.1")
  stop("diagFDR 0.2.1 or later is required (found ", DFDR_VERSION, ").", call. = FALSE)

## Exports of diagFDR 0.2.1 (checked by 01_preflight.R).
DFDR_EXPORTS_EXPECTED <- c("as_dfdr_tbl","collapse_universe","dfdr_accepted","dfdr_alert_basis","dfdr_bh_diagnostics","dfdr_bh_elasticity","dfdr_check_cv_benchmark","dfdr_curve_stability","dfdr_elasticity","dfdr_entrapment","dfdr_entrapment_discordant","dfdr_entrapment_size_ratio","dfdr_equal_chance_qbands","dfdr_headline","dfdr_ipe","dfdr_is_entrapment_group","dfdr_label_rule_audit","dfdr_local_tail","dfdr_p_calibration","dfdr_p_interpretation","dfdr_p_summary","dfdr_pep_decoy_sanity","dfdr_pep_interpretation","dfdr_pep_reliability","dfdr_pep_reliability_tdc","dfdr_pep_summary","dfdr_pi0_storey","dfdr_picked_protein","dfdr_plot_bh_elasticity","dfdr_plot_cv","dfdr_plot_dalpha","dfdr_plot_dwin","dfdr_plot_elasticity","dfdr_plot_equal_chance","dfdr_plot_fdrhat_pi0","dfdr_plot_p_calibration","dfdr_plot_p_calibration2","dfdr_plot_p_density_by_decoy","dfdr_plot_pep_density_by_decoy","dfdr_plot_pep_reliability","dfdr_plot_pep_reliability_tdc","dfdr_plot_pi0","dfdr_plot_scope_disagreement_matrix","dfdr_plot_score_distributions","dfdr_realised_fdp","dfdr_render_report","dfdr_run_all","dfdr_run_jaccard","dfdr_scope_disagreement","dfdr_simulate_tdc","dfdr_summary_headline","dfdr_sumpep","dfdr_write_manifest","dfdr_write_readme","dfdr_write_report","diann_entrapment_ids","diann_global_minrunq","diann_global_precursor","diann_global_proteingroup","diann_runxprecursor","flag_headline","label_is_intrinsic","mokapot_competed_universe","read_dfdr_maxquant_msms","read_dfdr_mzid","read_diann_parquet","read_fasta_accessions","read_mokapot_psms","read_spectronaut","read_spectronaut_efficient","score_to_pvalue","spectronaut_runxprecursor")

## ---------------------------------------------------------------------
##  1. Notation bridge (unchanged)
##      FDR_hat -> decoy ratio D_alpha/T_alpha; support_index -> D^-1/2
## ---------------------------------------------------------------------
headline_of <- function(x, alpha = ALPHA_MAIN) {
  h <- diagFDR::dfdr_headline(x, alpha = alpha)
  data.frame(alpha = h$alpha, T_alpha = h$T_alpha, D_alpha = h$D_alpha,
             decoy_ratio = h$FDR_hat, CV = h$support_index,
             CV_benchmark = h$CV_benchmark,
             FDR_minus1 = h$FDR_minus1, FDR_plus1 = h$FDR_plus1,
             FDR_minusK = h$FDR_minusK, FDR_plusK = h$FDR_plusK,
             FDR_minus2sqrtD = h$FDR_minus2sqrtD, FDR_plus2sqrtD = h$FDR_plus2sqrtD)
}
stability_of <- function(x, alphas = ALPHA_GRID) {
  s <- diagFDR::dfdr_curve_stability(x, alphas = alphas)
  data.frame(alpha = s$alpha, T_alpha = s$T_alpha, D_alpha = s$D_alpha,
             decoy_ratio = s$FDR_hat, CV = s$support_index)
}
local_tail_of <- function(x, alphas = ALPHA_GRID, win_rel = RHO, truncation = "drop") {
  w <- diagFDR::dfdr_local_tail(x, alphas = alphas, win_rel = win_rel, truncation = truncation)
  data.frame(alpha = w$alpha, D_win = w$D_alpha_win, n_win = w$n_win,
             decoy_frac_win = w$decoy_frac_win, truncated = w$truncated)
}
elasticity_of <- function(x, alphas = ALPHA_GRID, eps = EPSILON) {
  e <- diagFDR::dfdr_elasticity(x, alphas = alphas, eps = eps)
  data.frame(alpha = e$alpha, eps = e$eps, S = e$jaccard, n1 = e$n1, n2 = e$n2)
}
scope_disagreement_of <- function(x1, x2, alphas = ALPHA_GRID, label1 = "A", label2 = "B") {
  d <- diagFDR::dfdr_scope_disagreement(x1, x2, alphas = alphas, label1 = label1,
                                        label2 = label2, normalize_ids = TRUE)
  data.frame(alpha = d$alpha, J = d$jaccard, n1 = d$n1, n2 = d$n2)
}
sumpep_of <- function(x, alphas = ALPHA_GRID) {
  s <- diagFDR::dfdr_sumpep(x, alphas = alphas)
  data.frame(alpha = s$alpha, n_accepted = s$n_targets, sum_pep = s$sum_PEP, mean_pep = s$mean_PEP)
}
equal_chance_of <- function(x, breaks = EQ_BREAKS, low_conf = EQ_BAND, min_N = 2000) {
  e <- diagFDR::dfdr_equal_chance_qbands(x, breaks = breaks, low_conf = low_conf, min_N = min_N)
  list(bands = as.data.frame(e$bands), pooled = as.data.frame(e$pooled))
}
pep_reliability_of <- function(x, denom = c("tdc", "pooled"), n_min = NMIN_BIN,
                               pep_max = PEP_MAX, binwidth = PEP_BINWIDTH, add_decoy = 0L) {
  denom <- match.arg(denom)
  if (denom == "pooled") {
    r <- diagFDR::dfdr_pep_reliability(x, binwidth = binwidth, n_min = n_min, pep_max = pep_max)
    b <- as.data.frame(r$bins); b <- b[is.finite(b$pep_mean), , drop = FALSE]
    return(list(bins = data.frame(bin = as.character(b$bin), n = b$n, pep_mean = b$pep_mean,
                                  e_b = b$decoy_rate), IPE = r$IPE, denom = "pooled"))
  }
  r <- diagFDR::dfdr_ipe(x, binwidth = binwidth, pep_max = pep_max, n_min = n_min, add_decoy = add_decoy)
  b <- as.data.frame(r$bins)
  list(bins = data.frame(bin = b$bin, n = b$n, n_target = b$n_target, n_decoy = b$n_decoy,
                         pep_mean = b$pep_mean_targets, e_b = b$e_b, used = b$used),
       IPE = r$IPE, denom = "tdc")
}

## ---------------------------------------------------------------------
##  2. Label propagation --- now in the package. Aliases for scripts 02-03.
## ---------------------------------------------------------------------
as_decoy_logical <- diagFDR:::as_decoy_logical

sanitise_report <- function(rep, decoy_col = "Decoy", label = "") {
  d <- rep[[decoy_col]]; bad <- is.na(d)
  if (any(bad)) {
    message(sprintf("  [clean] %s: dropping %s row(s) with NA %s", label,
                    format(sum(bad), big.mark = ","), decoy_col))
    rep <- rep[!bad, , drop = FALSE]
  }
  rep[[decoy_col]] <- as.integer(as_decoy_logical(rep[[decoy_col]], decoy_col))
  rep
}

## Same signature as before; the work is done by diagFDR::collapse_universe().
collapse_universe <- function(rep, key_col, q_col,
                              rule = c("any_decoy", "winning_row", "majority", "split"),
                              pep_col = "PEP", q_max_export = NA_real_,
                              unit = "precursor", scope = "global", drop_empty_key = TRUE) {
  rule <- match.arg(rule)
  diagFDR::collapse_universe(rep, key_col, q_col, rule = rule, decoy_col = "Decoy",
                             pep_col = pep_col, q_max_export = q_max_export,
                             unit = unit, scope = scope, q_source = q_col,
                             tool = "DIA-NN", drop_empty_key = drop_empty_key)
}
label_is_intrinsic <- function(rep, key_col) diagFDR::label_is_intrinsic(rep, key_col)

## Picked-protein competition (03b) --- now diagFDR::dfdr_picked_protein().
picked_protein <- function(df, alpha = ALPHA_MAIN, offset = 0) {
  p <- diagFDR::dfdr_picked_protein(df, alpha = alpha, by = "q", add_decoy = offset)
  p$ratio <- p$decoy_ratio; p$cv <- p$CV_benchmark
  p
}
