## =====================================================================
##  Thin adapter over diagFDR (version 0.1.2 or later)
##
##  Every signature below was checked against the package source:
##     class_dfdr_tbl.R, io_diann.R, io_mokapot.R,
##     diag_stability.R, diag_calibration.R, diag_scope.R, run_all.R
##
##  It exists for two reasons only:
##    1. the package returns column names (FDR_hat, support_index,
##       D_alpha_win, jaccard) that differ from the paper's notation;
##    2. the package implements ONE label-propagation rule
##       (decoy-if-any). The paper needs three, so the alternatives are
##       built here.
##
##  There are no fallbacks. If diagFDR is missing, the scripts stop.
## =====================================================================

if (!requireNamespace("diagFDR", quietly = TRUE))
  stop("diagFDR is not installed. install.packages('diagFDR')", call. = FALSE)

library(diagFDR)

DFDR_VERSION <- as.character(utils::packageVersion("diagFDR"))
if (utils::packageVersion("diagFDR") < "0.1.2")
  stop("diagFDR 0.1.2 or later is required (found ", DFDR_VERSION, ").", call. = FALSE)


## Exports of diagFDR used or checked by these scripts.
DFDR_EXPORTS_EXPECTED <- c(
  "as_dfdr_tbl","dfdr_bh_diagnostics","dfdr_bh_elasticity","dfdr_curve_stability",
  "dfdr_elasticity","dfdr_equal_chance_qbands","dfdr_headline","dfdr_ipe","dfdr_local_tail",
  "dfdr_p_calibration","dfdr_pep_decoy_sanity","dfdr_pep_reliability",
  "dfdr_pep_reliability_tdc","dfdr_pi0_storey","dfdr_plot_bh_elasticity",
  "dfdr_plot_cv","dfdr_plot_dalpha","dfdr_plot_dwin","dfdr_plot_elasticity",
  "dfdr_plot_equal_chance","dfdr_plot_fdrhat_pi0","dfdr_plot_p_calibration",
  "dfdr_plot_p_calibration2","dfdr_plot_p_density_by_decoy",
  "dfdr_plot_pep_density_by_decoy","dfdr_plot_pep_reliability",
  "dfdr_plot_pep_reliability_tdc","dfdr_plot_pi0",
  "dfdr_plot_scope_disagreement_matrix","dfdr_plot_score_distributions",
  "dfdr_render_report","dfdr_run_all","dfdr_run_jaccard","dfdr_scope_disagreement",
  "dfdr_summary_headline","dfdr_sumpep","dfdr_write_manifest","dfdr_write_readme",
  "dfdr_write_report","diann_global_minrunq","diann_global_precursor",
  "diann_runxprecursor","flag_headline","mokapot_competed_universe",
  "read_dfdr_maxquant_msms","read_dfdr_mzid","read_diann_parquet",
  "read_mokapot_psms","read_spectronaut","read_spectronaut_efficient",
  "score_to_pvalue","spectronaut_runxprecursor")

## ---------------------------------------------------------------------
##  1. Notation bridge
##
##  package        paper
##  FDR_hat    ->  decoy ratio on the reported list, D_alpha / T_alpha
##  support_index -> decoy-support index D_alpha^{-1/2} (column CV here)
##  D_alpha_win->  D_{alpha,win}
##  jaccard    ->  S_alpha(eps)  /  J_run/global
## ---------------------------------------------------------------------

headline_of <- function(x, alpha = ALPHA_MAIN) {
  h <- diagFDR::dfdr_headline(x, alpha = alpha)
  data.frame(alpha = h$alpha, T_alpha = h$T_alpha, D_alpha = h$D_alpha,
             decoy_ratio = h$FDR_hat, CV = h$support_index,
             FDR_minus1 = h$FDR_minus1, FDR_plus1 = h$FDR_plus1,
             FDR_minusK = h$FDR_minusK, FDR_plusK = h$FDR_plusK,
             FDR_minus2sqrtD = h$FDR_minus2sqrtD,
             FDR_plus2sqrtD = h$FDR_plus2sqrtD)
}

stability_of <- function(x, alphas = ALPHA_GRID) {
  s <- diagFDR::dfdr_curve_stability(x, alphas = alphas)
  data.frame(alpha = s$alpha, T_alpha = s$T_alpha, D_alpha = s$D_alpha,
             decoy_ratio = s$FDR_hat, CV = s$support_index)
}

## NOTE: dfdr_local_tail() returns NA and warns when
## alpha * (1 + win_rel) exceeds the universe's q_max_export. Build the
## universes with the correct q_max_export or the wide-alpha rows vanish.
local_tail_of <- function(x, alphas = ALPHA_GRID, win_rel = RHO,
                          truncation = "drop") {
  w <- diagFDR::dfdr_local_tail(x, alphas = alphas, win_rel = win_rel,
                                truncation = truncation)
  data.frame(alpha = w$alpha, D_win = w$D_alpha_win, n_win = w$n_win,
             decoy_frac_win = w$decoy_frac_win, truncated = w$truncated)
}

elasticity_of <- function(x, alphas = ALPHA_GRID, eps = EPSILON) {
  e <- diagFDR::dfdr_elasticity(x, alphas = alphas, eps = eps)
  data.frame(alpha = e$alpha, eps = e$eps, S = e$jaccard, n1 = e$n1, n2 = e$n2)
}

scope_disagreement_of <- function(x1, x2, alphas = ALPHA_GRID,
                                  label1 = "A", label2 = "B") {
  ## normalize_ids = TRUE strips a "run||" prefix, so a run x precursor
  ## universe and a global one are compared on the same identifier space.
  d <- diagFDR::dfdr_scope_disagreement(x1, x2, alphas = alphas,
                                        label1 = label1, label2 = label2,
                                        normalize_ids = TRUE)
  data.frame(alpha = d$alpha, J = d$jaccard, n1 = d$n1, n2 = d$n2)
}

sumpep_of <- function(x, alphas = ALPHA_GRID) {
  s <- diagFDR::dfdr_sumpep(x, alphas = alphas)
  data.frame(alpha = s$alpha, n_accepted = s$n_targets,
             sum_pep = s$sum_PEP, mean_pep = s$mean_PEP)
}

equal_chance_of <- function(x, breaks = EQ_BREAKS, low_conf = EQ_BAND,
                            min_N = 2000) {
  e <- diagFDR::dfdr_equal_chance_qbands(x, breaks = breaks,
                                         low_conf = low_conf, min_N = min_N)
  list(bands = as.data.frame(e$bands), pooled = as.data.frame(e$pooled))
}

## ---------------------------------------------------------------------
##  2. The IPE
##
##  dfdr_pep_reliability()      -> list(bins, IPE, params);
##      bins$decoy_rate = D_b / (T_b + D_b)   <- the POOLED comparator
##  dfdr_pep_reliability_tdc()  -> bin table only, no IPE;
##      err_hat_target = (D_b + add_decoy) / T_b   <- the CORRECT one
##
##  The paper defines IPE as a weighted mean absolute deviation over
##  bins with mean PEP <= pep_max and n_b >= n_min. The package supplies
##  that weighting only for the pooled version, so we apply the same rule
##  to the TDC bins here. Both are returned so the difference can be
##  reported.
## ---------------------------------------------------------------------

pep_reliability_of <- function(x, denom = c("tdc", "pooled"),
                               n_min = NMIN_BIN, pep_max = PEP_MAX,
                               binwidth = PEP_BINWIDTH, add_decoy = 0L) {
  denom <- match.arg(denom)

  if (denom == "pooled") {
    r <- diagFDR::dfdr_pep_reliability(x, binwidth = binwidth,
                                       n_min = n_min, pep_max = pep_max)
    b <- as.data.frame(r$bins)
    b <- b[is.finite(b$pep_mean), , drop = FALSE]
    return(list(bins = data.frame(bin = as.character(b$bin), n = b$n,
                                  pep_mean = b$pep_mean, e_b = b$decoy_rate),
                IPE = r$IPE, denom = "pooled"))
  }

  ## IPE of the paper (Eq. 1), computed by diagFDR::dfdr_ipe()
  r <- diagFDR::dfdr_ipe(x, binwidth = binwidth, pep_max = pep_max,
                         n_min = n_min, add_decoy = add_decoy)
  b <- as.data.frame(r$bins)
  list(bins = data.frame(bin = b$bin, n = b$n, n_target = b$n_target,
                         n_decoy = b$n_decoy, pep_mean = b$pep_mean_targets,
                         e_b = b$e_b, used = b$used),
       IPE = r$IPE, denom = "tdc")
}

## ---------------------------------------------------------------------
##  3. Label propagation
##
##  diann_global_precursor(), diann_runxprecursor() and
##  diann_global_minrunq() all collapse with  is_decoy = any(is_decoy).
##  That is the paper's PRIMARY rule. The two alternatives are built here.
##
##  WARNING, and it matters. At PRECURSOR level the decoy label is
##  effectively intrinsic to Precursor.Id, so the three rules coincide.
##  At PROTEIN-GROUP level it is not: a group routinely contains both
##  target and decoy precursors, and decoy-if-any then labels almost
##  every group a decoy. Use rule = "split" at protein level.
## ---------------------------------------------------------------------

## Robust coercion of a DIA-NN Decoy column to logical.
## arrow may hand back integer, integer64 or double; a re-exported TSV may
## give character or factor. as.integer() on a factor returns LEVEL CODES
## (1/2), which would silently invert the labels, so factors are handled
## explicitly.
as_decoy_logical <- function(v, where = "Decoy") {
  if (is.logical(v)) return(v)
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) {
    u <- unique(v[!is.na(v)])
    if (!all(u %in% c("0","1","TRUE","FALSE","True","False","true","false")))
      stop(where, ": unexpected character values: ", paste(head(u, 5), collapse = ", "))
    return(v %in% c("1","TRUE","True","true"))
  }
  vn <- suppressWarnings(as.numeric(v))
  u  <- unique(vn[!is.na(vn)])
  if (!all(u %in% c(0, 1)))
    stop(where, ": expected 0/1, found ", paste(head(sort(u), 5), collapse = ", "))
  vn == 1
}

## Drop rows the package would choke on, and say how many. Passing a report
## with NA in Decoy straight into diann_*() produces `any(is_decoy) = NA`
## and a validate_dfdr_tbl() abort several steps later, far from the cause.
sanitise_report <- function(rep, decoy_col = "Decoy", label = "") {
  d <- rep[[decoy_col]]
  bad <- is.na(d)
  if (any(bad)) {
    message(sprintf("  [clean] %s: dropping %s row(s) with NA %s",
                    label, format(sum(bad), big.mark = ","), decoy_col))
    rep <- rep[!bad, , drop = FALSE]
  }
  rep[[decoy_col]] <- as.integer(as_decoy_logical(rep[[decoy_col]], decoy_col))
  rep
}

collapse_universe <- function(rep, key_col, q_col,
                              rule = c("any_decoy", "winning_row", "majority", "split"),
                              pep_col = "PEP", q_max_export = NA_real_,
                              unit = "precursor", scope = "global",
                              drop_empty_key = TRUE) {
  rule <- match.arg(rule)
  need <- c(key_col, "Decoy", q_col)
  miss <- setdiff(need, names(rep))
  if (length(miss)) stop("collapse_universe: missing columns ",
                         paste(miss, collapse = ", "))

  key <- as.character(rep[[key_col]])
  q   <- suppressWarnings(as.numeric(rep[[q_col]]))
  dc  <- as_decoy_logical(rep[["Decoy"]])
  pep <- if (!is.null(pep_col) && pep_col %in% names(rep))
    suppressWarnings(as.numeric(rep[[pep_col]])) else NULL

  keep <- !is.na(key) & !is.na(dc) & is.finite(q) & q >= 0 & q <= 1
  if (drop_empty_key) keep <- keep & nzchar(key)
  n_drop <- sum(!keep)
  if (n_drop > 0 && n_drop / length(keep) > 0.5)
    warning(sprintf("collapse_universe(%s, %s): dropping %.0f%% of rows",
                    key_col, q_col, 100 * n_drop / length(keep)), call. = FALSE)
  key <- key[keep]; q <- q[keep]; dc <- dc[keep]
  if (!is.null(pep)) pep <- pep[keep]
  if (!length(key)) stop("collapse_universe(", key_col, ", ", q_col,
                         "): no usable rows after filtering.")

  ## `method = "radix"` sorts in C byte order. It is locale-independent, so
  ## identical keys are guaranteed adjacent (a locale collation that treats
  ## two distinct strings as equal can otherwise interleave their rows and
  ## break the "first row per key is the minimum q" assumption). It is also
  ## several times faster on a million rows.
  if (rule == "split") {
    ## Target and decoy entries of the same key are DIFFERENT hypotheses.
    ## DIA-NN's own convention at protein-group level.
    grp_key <- paste0(ifelse(dc, "DECOY:", ""), key)
  } else {
    grp_key <- key
  }
  o <- order(grp_key, q, method = "radix")
  grp_key <- grp_key[o]; q <- q[o]; dc <- dc[o]
  if (!is.null(pep)) pep <- pep[o]

  first <- !duplicated(grp_key)
  idk   <- grp_key[first]

  ## Positional group aggregation. `grp` is 1..G with no gaps and in the
  ## same order as `idk`, so rowsum(reorder = FALSE) lines up element by
  ## element. No name lookup, so no possibility of an unmatched key
  ## silently becoming NA.
  lab <- if (rule == "split" || rule == "winning_row") {
    dc[first]
  } else {
    grp   <- cumsum(first)
    n_dec <- as.vector(rowsum(as.integer(dc), grp, reorder = FALSE))
    n_tot <- tabulate(grp, nbins = length(idk))
    if (rule == "any_decoy") n_dec > 0L else (n_dec / n_tot) > 0.5
  }

  if (length(lab) != length(idk))
    stop("collapse_universe(", key_col, ", ", q_col, ", rule=", rule,
         "): label length ", length(lab), " != group count ", length(idk))
  if (anyNA(lab))
    stop("collapse_universe(", key_col, ", ", q_col, ", rule=", rule, "): ",
         sum(is.na(lab)), " of ", length(lab),
         " groups got an NA label. This should be impossible after ",
         "filtering -- please report the offending file.")

  out <- data.frame(id = idk, is_decoy = as.logical(lab), q = q[first],
                    stringsAsFactors = FALSE)
  out$pep   <- if (!is.null(pep)) pep[first] else NA_real_
  out$run   <- NA_character_
  out$score <- NA_real_
  out <- out[, c("id", "is_decoy", "q", "pep", "run", "score")]

  diagFDR::as_dfdr_tbl(out, unit = unit, scope = scope, q_source = q_col,
                       q_max_export = q_max_export,
                       provenance = list(tool = "DIA-NN", label_rule = rule,
                                         key = key_col, n_dropped = n_drop))
}

## Is the decoy label intrinsic to the grouping key? If yes, the three
## propagation rules coincide and the corresponding audit is unnecessary.
label_is_intrinsic <- function(rep, key_col) {
  key <- as.character(rep[[key_col]])
  dc  <- as.integer(as_decoy_logical(rep[["Decoy"]]))
  ok  <- !is.na(key) & !is.na(dc) & nzchar(key)
  n_lab <- tapply(dc[ok], key[ok], function(z) length(unique(z)))
  list(intrinsic = sum(n_lab > 1) == 0,
       n_mixed = sum(n_lab > 1), n_keys = length(n_lab))
}
