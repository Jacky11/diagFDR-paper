## =====================================================================
##  06 --- MS2Rescore case study: Figure 5, Section 4.5
##
##  Three q-value constructions on the same PSMs:
##    (a) BEFORE  target-decoy counting on the Andromeda score, via
##                read_dfdr_maxquant_msms() on msms.txt
##    (b) MOKAPOT the exported "mokapot q-value" -- a (D+1)/T count over all
##        31,314 candidate PSMs (up to ten per spectrum), applied here to the
##        spectrum-level winners (see 06c_mokapot_universe.R)
##    (c) AFTER   counting recomputed from the "mokapot score"
##
##  (b) versus (c) is the point of the section. In (c) the q-value IS a
##  decoy ratio on the reported universe (the winners), so D_alpha is
##  pinned to alpha*T_alpha. In (b) it is a decoy ratio on another
##  universe -- all candidate PSMs -- so on the winners it is not: at 1%
##  it accepts 43 decoys for 7,070 targets (0.61%). The loss of support is
##  a scope effect, not a consequence of the improved discrimination.
##
##  The mokapot export carries NO pre-rescoring score column (verified:
##  index, spectrum_id, peptide, is_target, run, calcmass, expmass,
##  retention_time, charge, mokapot score, mokapot q-value, mokapot PEP,
##  protein_list), which is why (a) comes from msms.txt.
##
##  UNIVERSE CAVEAT, state it in the paper. (b) and (c) are competed to
##  one winner per run x spectrum_id by mokapot_competed_universe().
##  read_dfdr_maxquant_msms() returns PSM-level rows keyed on the MaxQuant
##  `id`, which is not the same universe. Counts are therefore comparable
##  in trend but not row-for-row.
## =====================================================================

source("00_config.R")
load_numbers()

need_file(MOKAPOT_TARGETS, "mokapot target PSMs")
need_file(MOKAPOT_DECOYS,  "mokapot decoy PSMs")

## --- (b) and (c): the competed mokapot universe -----------------------
raw <- diagFDR::read_mokapot_psms(target_path = MOKAPOT_TARGETS,
                                  decoy_path  = MOKAPOT_DECOYS)
say("%s mokapot PSM rows", format(nrow(raw), big.mark = ","))

x_mok <- diagFDR::mokapot_competed_universe(
  raw, run_col = "run", spectrum_id_col = "spectrum_id",
  score_col = "mokapot score", q_col = "mokapot q-value", pep_col = "mokapot PEP",
  id_mode = "runxid", unit = "psm", scope = "global",
  q_source = "mokapot q-value", q_max_export = 1)
say("competed to %s winners (%s decoys)", format(nrow(x_mok), big.mark = ","),
    format(sum(x_mok$is_decoy), big.mark = ","))

## (c) recount from the rescored score on the SAME competed universe
x_after <- x_mok
## offset = 1 matches read_dfdr_maxquant_msms(add_decoy = 1L), so arms (a)
## and (c) use the same estimator and differ only in the score.
x_after$q <- tdc_qvalues(x_after$score, x_after$is_decoy, offset = 1L)
x_after <- diagFDR::as_dfdr_tbl(
  as.data.frame(x_after[, c("id","is_decoy","q","pep","run","score")]),
  unit = "psm", scope = "global", q_source = "TDC on mokapot score",
  q_max_export = 1, provenance = list(tool = "mokapot", note = "recounted"))

## --- (a) before rescoring, from msms.txt ------------------------------
x_before <- NULL
if (file.exists(MAXQUANT_MSMS)) {
  x_before <- diagFDR::read_dfdr_maxquant_msms(
    MAXQUANT_MSMS, pep_mode = "drop", exclude_contaminants = TRUE, add_decoy = 1L,
    unit = "psm", scope = "global", provenance = list(tool = "MaxQuant"))
  say("msms.txt: %s PSMs (%s decoys)", format(nrow(x_before), big.mark = ","),
      format(sum(x_before$is_decoy), big.mark = ","))
} else {
  message("[note] msms.txt not found -- the before-rescoring arm is omitted.")
}

TB <- list()
if (!is.null(x_before)) TB[["before rescoring (Andromeda, TDC)"]] <- x_before
TB[["mokapot q-value (exported, all candidates)"]] <- x_mok
TB[["after rescoring (mokapot score, TDC)"]] <- x_after
N <- names(TB)

## --- headline, stability, local support -------------------------------
res <- do.call(rbind, lapply(N, function(n) {
  s <- stability_of(TB[[n]], ALPHA_GRID)
  w <- local_tail_of(TB[[n]], ALPHA_GRID, RHO, truncation = "drop")
  cbind(construction = n, merge(s, w[, c("alpha","D_win","truncated")], by = "alpha"))
}))
print(res[res$alpha %in% c(5e-4, 0.01, 0.1), ])
write_tab(res, "msrescore_stability")

gm <- function(n, a, col) { v <- res[res$construction == n & res$alpha == a, col]
                            if (!length(v)) NA else v[1] }
nb <- if (length(N) == 3) N[1] else NA
nm <- N[length(N) - 1]; na_ <- N[length(N)]
if (!is.na(nb)) {
  reg("psmBeforeOne",    gm(nb, 0.01, "T_alpha"), "int", "Sec 4.5")
  reg("psmBeforeStrict", gm(nb, 5e-4, "T_alpha"), "int", "Sec 4.5")
  reg("DbeforeOne",      gm(nb, 0.01, "D_alpha"), "int", "Sec 4.5")
  reg("CVbeforeOne",     gm(nb, 0.01, "CV"), "d3", "Sec 4.5")
}
reg("psmMokapotOne",    gm(nm, 0.01, "T_alpha"), "int", "Sec 4.5")
reg("psmScoreOne",      gm(na_, 0.01, "T_alpha"), "int", "Sec 4.5")
reg("psmMokapotStrict", gm(nm, 5e-4, "T_alpha"), "int", "Sec 4.5")
reg("psmScoreStrict",   gm(na_, 5e-4, "T_alpha"), "int", "Sec 4.5")
reg("DmokapotOne",      gm(nm, 0.01, "D_alpha"), "int", "Sec 4.5: exported q counted on all candidates, NOT pinned on winners")
reg("DafterOne",        gm(na_, 0.01, "D_alpha"), "int", "Sec 4.5: counting q, pinned to alpha*T")
reg("CVmokapotOne",     gm(nm, 0.01, "CV"), "d3", "Sec 4.5")
reg("CVafterOne",       gm(na_, 0.01, "CV"), "d3", "Sec 4.5")
reg("ratioMokapotOne",  gm(nm, 0.01, "decoy_ratio"), "pct2",
    "Sec 4.5: mokapot q (computed on all candidates) gives 0.61% on the winners")
reg("DwinMokapot",      gm(nm, 0.01, "D_win"), "int", "Sec 4.5")
reg("DwinScore",        gm(na_, 0.01, "D_win"), "int", "Sec 4.5")

## --- overlap before vs after ------------------------------------------
if (!is.na(nb)) {
  J1 <- scope_disagreement_of(TB[[nb]], TB[[nm]], ALPHA_GRID, "before", "mokapot q")
  J2 <- scope_disagreement_of(TB[[nb]], TB[[na_]], ALPHA_GRID, "before", "rescored TDC")
  ov <- data.frame(alpha = J1$alpha, J_mokapot = J1$J, J_score = J2$J)
  print(ov); write_tab(ov, "msrescore_overlap")
  reg("JrescoreOne", mean(c(ov$J_mokapot[ov$alpha==0.01], ov$J_score[ov$alpha==0.01])), "d2",
      "Sec 4.5 list overlap at 1%")
  reg("JrescoreStrict", mean(c(ov$J_mokapot[ov$alpha==5e-4], ov$J_score[ov$alpha==5e-4])), "d2",
      "Sec 4.5 list overlap at 5e-4")
} else ov <- NULL

## --- calibration -------------------------------------------------------
rel <- pep_reliability_of(x_mok, "tdc")
if (nrow(rel$bins)) {
  reg("mokapotPEPbin", rel$bins$pep_mean[1], "d5", "Sec 4.5 lowest PEP bin")
  reg("mokapotEb",     rel$bins$e_b[1], "d5", "Sec 4.5, D_b/T_b comparator")
  reg("mokapotIPE",    rel$IPE, "d4", "Sec 4.5 IPE, correct comparator")
  write_tab(rel$bins, "msrescore_pep_reliability")
}
spm <- sumpep_of(x_mok, ALPHA_MAIN)
reg("mokapotSumPEP", round(spm$sum_pep), "int", "Sec 4.5")
eqm <- equal_chance_of(x_mok)
write_tab(eqm$pooled, "msrescore_equal_chance")
write_tab(eqm$bands,  "msrescore_equal_chance_bands")
if ("pi_D_hat" %in% names(eqm$pooled)) {
  reg("mokapotPiD",  eqm$pooled$pi_D_hat[1], "d3", "Sec 4.5 decoy fraction in q in [0.2,0.5]")
  reg("mokapotPiDp", eqm$pooled$p_value_binom[1], NULL, "Sec 4.5 binomial p against 0.5")
}

## --- Figure 5 ----------------------------------------------------------
p <- open_pdf("fig4_msrescore", 11, 6.6)
par(mfrow = c(2, 3), mar = c(4.2, 4.6, 2.2, 1))
dt <- density(x_mok$score[!x_mok$is_decoy]); dd <- density(x_mok$score[x_mok$is_decoy])
plot(dt, main = "A: mokapot score", xlab = "score", ylim = range(c(dt$y, dd$y)))
lines(dd, col = 2, lty = 2)
legend("topright", c("target","decoy"), col = 1:2, lty = 1:2, bty = "n", cex = .7)
plot(range(res$alpha), range(res$CV[is.finite(res$CV)]), type = "n", log = "x",
     xlab = expression(alpha), ylab = expression(hat(CV)[alpha]), main = "B")
for (k in seq_along(N)) { d <- res[res$construction == N[k], ]; lines(d$alpha, d$CV, col = k, lwd = 1.6) }
legend("topright", N, col = seq_along(N), lwd = 1.6, bty = "n", cex = .55)
dw <- res[is.finite(res$D_win), ]
plot(range(dw$alpha), range(pmax(dw$D_win,1)), type = "n", log = "xy",
     xlab = expression(alpha), ylab = expression(D[list(alpha,win)]), main = "C")
for (k in seq_along(N)) { d <- dw[dw$construction == N[k], ]; if (nrow(d)) lines(d$alpha, pmax(d$D_win,1), col = k, lwd = 1.6) }
plot(range(res$alpha), range(res$T_alpha), type = "n", log = "x",
     xlab = expression(alpha), ylab = "accepted target PSMs", main = "D")
for (k in seq_along(N)) { d <- res[res$construction == N[k], ]; lines(d$alpha, d$T_alpha, col = k, lwd = 1.6) }
if (!is.null(ov)) {
  plot(ov$alpha, ov$J_mokapot, type = "b", pch = 19, log = "x", ylim = c(0,1),
       xlab = expression(alpha), ylab = "Jaccard vs before", main = "E")
  lines(ov$alpha, ov$J_score, type = "b", pch = 17, col = 3)
} else plot.new()
b <- eqm$bands
plot(b$q_mean, b$decoy_frac, type = "b", pch = 19, ylim = c(0, .6),
     xlab = "mean q in band", ylab = "decoy fraction", main = "F")
abline(h = 0.5, lty = 2)
dev.off(); message("  -> ", p)

save_numbers()
message("\nNext: 06b_msrescore_overlap.R")
