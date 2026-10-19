## =====================================================================
##  04 --- Calibration and stability: the panels of Figure 3 (A-F) and
##         every IPE of Section 4.3
##
##  Both comparators are computed so the difference can be reported:
##    dfdr_pep_reliability()      decoy_rate     = D_b/(T_b+D_b)   POOLED
##    dfdr_pep_reliability_tdc()  err_hat_target = D_b/T_b         CORRECT
##  The package returns the TDC bin table WITHOUT an IPE; the weighting
##  (bins with mean PEP <= PEP_MAX and n_b >= NMIN_BIN) is applied in
##  R/dfdr_api.R so that both IPEs use an identical rule.
##
##  NOTE on the equal-chance flag. These are --reanalyse (two-pass)
##  exports, so the low-confidence region is conditioned on passing an
##  intermediate criterion. The statistic is reported; the flag is marked
##  NOT EVALUABLE (Section 3.2 of the paper).
## =====================================================================

source("00_config.R")
load_numbers()
U <- readRDS(file.path(OUT, "universes.rds"))

UNIVERSES <- list()
for (s in c("basic","semi")) for (u in c("runxprec","global","minrunq")) {
  if (!is.null(U[[s]][[u]]))
    UNIVERSES[[paste(s, switch(u, runxprec="run x precursor", global="global",
                               minrunq="global min_r"))]] <- U[[s]][[u]]
}
for (s in c("basic_fp","semi_fp")) if (!is.null(U[[s]]$global))
  UNIVERSES[[paste(sub("_fp","",s), "global first pass")]] <- U[[s]]$global
say("%d universes", length(UNIVERSES))

## ---------------------------------------------------------------------
##  1. Stability curves --- Figure 3A/B/C
## ---------------------------------------------------------------------
message("\n--- stability ---")
stab <- do.call(rbind, lapply(names(UNIVERSES), function(n) {
  x <- UNIVERSES[[n]]
  s <- stability_of(x, ALPHA_GRID)
  ## truncation = "drop": rows whose window exceeds q_max_export return NA
  ## rather than a silently wrong count.
  w <- local_tail_of(x, ALPHA_GRID, RHO, truncation = "drop")
  e <- elasticity_of(x, ALPHA_GRID, EPSILON)
  m <- merge(merge(s, w[, c("alpha","D_win","truncated")], by = "alpha"),
             e[, c("alpha","S")], by = "alpha")
  cbind(universe = n, m)
}))
write_tab(stab, "stability_curves")
if (any(stab$truncated, na.rm = TRUE))
  say("%d (universe, alpha) rows dropped: window above q_max_export",
      sum(stab$truncated, na.rm = TRUE))

sb <- function(u, a, col) {
  v <- stab[stab$universe == u & stab$alpha == a, col]; if (!length(v)) NA else v[1]
}
reg("DalphaBasicGlobal", sb("basic global", 0.01, "D_alpha"), "int", "Sec 4.3")
reg("DalphaSemiGlobal",  sb("semi global",  0.01, "D_alpha"), "int", "Sec 4.3")
reg("CVbasicGlobal",     sb("basic global", 0.01, "CV"), "d3", "Sec 4.3")
reg("DalphaBasicStringent", sb("basic global", 1e-4, "D_alpha"), "int", "Sec 4.3 alpha=1e-4")
reg("CVbasicStringent",     sb("basic global", 1e-4, "CV"), "d2", "Sec 4.3")
reg("DwinBasic", sb("basic global", 0.01, "D_win"), "int", "Sec 4.3 local support")
reg("DwinSemi",  sb("semi global",  0.01, "D_win"), "int", "Sec 4.3")
reg("ElasBasic", sb("basic global", 0.01, "S"), "d3", "Sec 4.3 elasticity at 1%")

p <- open_pdf("fig2ABC_stability", 11, 3.8)
par(mfrow = c(1, 3), mar = c(4.2, 4.6, 2.2, 1))
cols <- seq_along(UNIVERSES)
plot(range(stab$alpha), range(pmax(stab$D_alpha, 1)), type = "n", log = "xy",
     xlab = expression(alpha), ylab = expression(D[alpha]), main = "A")
for (k in seq_along(UNIVERSES)) { d <- stab[stab$universe == names(UNIVERSES)[k], ]
  lines(d$alpha, pmax(d$D_alpha, 1), col = cols[k], lwd = 1.6) }
legend("topleft", names(UNIVERSES), col = cols, lwd = 1.6, bty = "n", cex = .5)
dw <- stab[is.finite(stab$D_win), ]
plot(range(dw$alpha), range(pmax(dw$D_win, 1)), type = "n", log = "xy",
     xlab = expression(alpha), ylab = expression(D[list(alpha, win)]), main = "B")
for (k in seq_along(UNIVERSES)) { d <- dw[dw$universe == names(UNIVERSES)[k], ]
  if (nrow(d)) lines(d$alpha, pmax(d$D_win, 1), col = cols[k], lwd = 1.6) }
plot(range(stab$alpha), c(0, 1), type = "n", log = "x",
     xlab = expression(alpha), ylab = expression(S[alpha](0.2)), main = "C")
for (k in seq_along(UNIVERSES)) { d <- stab[stab$universe == names(UNIVERSES)[k], ]
  lines(d$alpha, d$S, col = cols[k], lwd = 1.6) }
dev.off(); message("  -> ", p)

## ---------------------------------------------------------------------
##  2. PEP reliability and the IPE
## ---------------------------------------------------------------------
message("\n--- PEP reliability: both comparators ---")
ipe_rows <- NULL; rel_store <- list()
for (n in names(UNIVERSES)) {
  x <- UNIVERSES[[n]]
  if (!"pep" %in% names(x) || all(is.na(x$pep))) { say("%-28s no PEP, skipped", n); next }
  rt <- pep_reliability_of(x, "tdc")
  rp <- tryCatch(pep_reliability_of(x, "pooled"),
                 error = function(e) list(IPE = NA_real_, bins = NULL))
  rel_store[[n]] <- rt$bins
  ipe_rows <- rbind(ipe_rows, data.frame(universe = n, IPE_tdc = rt$IPE,
                                         IPE_pooled = rp$IPE,
                                         ratio = rt$IPE / rp$IPE,
                                         bins_used = sum(rt$bins$used, na.rm = TRUE)))
  say("%-28s IPE(D/T)=%.4f  IPE(pooled)=%.4f  x%.2f", n, rt$IPE, rp$IPE, rt$IPE/rp$IPE)
}
if (!is.null(ipe_rows)) {
  write_tab(ipe_rows, "ipe_both_comparators")
  gi <- function(u) { v <- ipe_rows$IPE_tdc[ipe_rows$universe == u]; if (!length(v)) NA else v[1] }
  reg("IPEbasicRunPrec", gi("basic run x precursor"), "d4",
      "Sec 4.3 and the verifiable-claim example. D_b/T_b comparator.")
  reg("IPEbasicGlobal",  gi("basic global"), "d4", "Sec 4.3")
  reg("IPEbasicMinrun",  gi("basic global min_r"), "d4", "Sec 4.3")
  reg("IPEsemiGlobal",   gi("semi global"), "d4", "Sec 4.3")

  lb <- function(n) { b <- rel_store[[n]]
    if (is.null(b) || !nrow(b)) c(NA, NA) else c(b$pep_mean[1], b$e_b[1]) }
  v <- lb("basic run x precursor"); reg("pepBinRunPrec", v[1], "d5", "Sec 4.3 lowest PEP bin")
  reg("ebBinRunPrec", v[2], "d5", "Sec 4.3")
  v <- lb("basic global min_r");    reg("pepBinMinrun", v[1], "d5", "Sec 4.3")
  reg("ebBinMinrun", v[2], "d5", "Sec 4.3")

  p <- open_pdf("fig2D_pep_reliability", 5.8, 4.8)
  ymax <- max(0.6, max(unlist(lapply(rel_store, function(b) b$e_b)), na.rm = TRUE))
  plot(c(0, PEP_MAX), c(0, ymax), type = "n",
       xlab = "mean predicted PEP among targets",
       ylab = expression(hat(e)[b] == D[b] / T[b]), main = "D")
  abline(0, 1, lty = 2); k <- 0
  for (n in names(rel_store)) { k <- k + 1; b <- rel_store[[n]]
    lines(b$pep_mean, b$e_b, col = k, lwd = 1.6); points(b$pep_mean, b$e_b, col = k, pch = 19, cex = .6) }
  legend("topleft", names(rel_store), col = seq_along(rel_store), lwd = 1.6, bty = "n", cex = .5)
  dev.off(); message("  -> ", p)
}

## ---------------------------------------------------------------------
##  3. Equal-chance q-bands --- Figure 3E, gated NOT EVALUABLE
## ---------------------------------------------------------------------
message("\n--- equal-chance ---")
eq <- do.call(rbind, lapply(names(UNIVERSES), function(n) {
  e <- equal_chance_of(UNIVERSES[[n]])
  cbind(universe = n, e$pooled)
}))
## Two-pass (--reanalyse) exports: conditional universe, so the flag is
## reported but NOT evaluated.
eq$multistage  <- TRUE
eq$flag_status <- "not evaluable (two-pass export)"
print(eq[, intersect(names(eq), c("universe","N_test","pi_D_hat","ci95_lo","ci95_hi","p_value_binom","flag_status"))])
write_tab(eq, "equal_chance")
gp <- function(u) { v <- eq$pi_D_hat[eq$universe == u]; if (!length(v)) NA else v[1] }
reg("piDbasicGlobal",  gp("basic global"), "d3", "Sec 4.3")
reg("piDbasicRunPrec", gp("basic run x precursor"), "d3", "Sec 4.3")
reg("piDbasicMinrun",  gp("basic global min_r"), "d3", "Sec 4.3")

bands <- do.call(rbind, lapply(names(UNIVERSES), function(n)
  cbind(universe = n, equal_chance_of(UNIVERSES[[n]])$bands)))
write_tab(bands, "equal_chance_bands")
p <- open_pdf("fig2E_equal_chance", 5.8, 4.4)
plot(c(0, 0.5), c(0, 1), type = "n", xlab = "mean q in band",
     ylab = "decoy fraction", main = "E (not evaluable: two-pass exports)")
abline(h = 0.5, lty = 2); k <- 0
for (n in unique(bands$universe)) { k <- k + 1; d <- bands[bands$universe == n, ]
  lines(d$q_mean, d$decoy_frac, col = k, lwd = 1.6, type = "b", pch = 19, cex = .6) }
dev.off(); message("  -> ", p)

## ---------------------------------------------------------------------
##  4. Absolute error burden --- Figure 3F
## ---------------------------------------------------------------------
message("\n--- sum(PEP) ---")
sp <- do.call(rbind, lapply(names(UNIVERSES), function(n) {
  x <- UNIVERSES[[n]]
  if (!"pep" %in% names(x) || all(is.na(x$pep))) return(NULL)
  cbind(universe = n, sumpep_of(x, ALPHA_GRID))
}))
if (!is.null(sp)) {
  write_tab(sp, "sumpep")
  gs <- function(u) { v <- sp$sum_pep[sp$universe == u & sp$alpha == ALPHA_MAIN]
                      if (!length(v)) NA else v[1] }
  reg("sumPEPbasicGlobal",  round(gs("basic global")), "int", "Sec 4.3")
  reg("sumPEPsemiGlobal",   round(gs("semi global")),  "int", "Sec 4.3")
  reg("sumPEPbasicRunPrec", round(gs("basic run x precursor")), "int", "Sec 4.3")
  reg("sumPEPsemiRunPrec",  round(gs("semi run x precursor")),  "int", "Sec 4.3")
  p <- open_pdf("fig2F_sumpep", 5.8, 4.4)
  plot(range(sp$alpha), range(pmax(sp$sum_pep, 1)), type = "n", log = "xy",
       xlab = expression(alpha), ylab = expression(sum(PEP)), main = "F")
  k <- 0; for (n in unique(sp$universe)) { k <- k + 1; d <- sp[sp$universe == n, ]
    lines(d$alpha, pmax(d$sum_pep, 1), col = k, lwd = 1.6) }
  dev.off(); message("  -> ", p)
}

## ---------------------------------------------------------------------
##  5. IPE under the two comparators
## ---------------------------------------------------------------------
if (!is.null(ipe_rows)) {
  message("\n=== IPE: pooled comparator against D_b/T_b ===")
  print(ipe_rows[, c("universe","IPE_pooled","IPE_tdc","ratio")])
}

save_numbers()
message("\nNext: 05_entrapment.R")
