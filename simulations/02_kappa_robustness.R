## =====================================================================
##  02 --- How stable is the variability factor kappa?
##
##  Proposition 1 predicts CV(FDP) = sqrt(2) * D^(-1/2) under correct
##  specification, i.e. kappa = sqrt(2) = 1.414. This script repeats
##  scenario S0 over five (pi0, delta) pairs and the three null families.
##  Cells with fewer than 10 accepted decoys on average are reported as
##  not evaluable.
##
##  It then checks that the inflation of the FDP under S1 and S2 does not
##  depend on the choice of (pi0, delta).
##
##  Writes out/tables/kappa_robustness.csv, table_kappa.tex and
##         out/tables/misspecification_by_parameters.csv
## =====================================================================

source("parameters.R")
cl <- make_cluster(N_CORES)

settings <- data.frame(pi0   = c(0.550, 0.600, 0.667, 0.667, 0.750),
                       delta = c(5.50,  4.00,  3.50,  4.65,  5.00))
MIN_MEAN_D <- 10

rows <- list(); k <- 0
for (fam in FAMILIES) for (j in seq_len(nrow(settings))) {
  k <- k + 1
  message(sprintf("[%2d/15] %-8s pi0=%.3f delta=%.2f", k, fam,
                  settings$pi0[j], settings$delta[j]))
  d <- run_cell(B_KAPPA, SEED + 2e6 + 1e4 * k, "S0", fam, N_UNITS,
                settings$pi0[j], settings$delta[j], ALPHAS, cl = cl)
  s <- summarise_cells(d)
  rows[[k]] <- cbind(pi0 = settings$pi0[j], delta = settings$delta[j],
                     s[c("family", "alpha", "B", "mean_D", "mean_fdp",
                         "cv_fdp", "mean_support_index", "kappa")])
}

## Do the misspecification results depend on (pi0, delta)? Scenarios S1
## and S2 with Gaussian nulls at alpha = 1%, for the same five pairs.
mis <- list(); m <- 0
for (sc in c("S1", "S2")) for (j in seq_len(nrow(settings))) {
  m <- m + 1
  message(sprintf("[%2d/10] %s pi0=%.3f delta=%.2f", m, sc,
                  settings$pi0[j], settings$delta[j]))
  d <- run_cell(B_GRID, SEED + 3e6 + 1e4 * m, sc, "gaussian", N_UNITS,
                settings$pi0[j], settings$delta[j], ALPHA_MAIN, cl = cl)
  s <- summarise_cells(d)
  mis[[m]] <- cbind(pi0 = settings$pi0[j], delta = settings$delta[j],
                    s[c("scenario", "mean_fdp", "mean_decoy_ratio",
                        "frac_bad", "fire_any")])
}
write_table(do.call(rbind, mis), "misspecification_by_parameters")
if (!is.null(cl)) parallel::stopCluster(cl)

kappa <- do.call(rbind, rows)
kappa$evaluable <- kappa$mean_D >= MIN_MEAN_D
write_table(kappa, "kappa_robustness")

ev <- kappa$kappa[kappa$evaluable]
message(sprintf("  evaluable cells: %d; kappa range %.2f-%.2f, median %.2f",
                length(ev), min(ev), max(ev), median(ev)))

family_label <- c(gaussian = "Gaussian", gumbel = "Gumbel", t4 = "Student-$t_4$")
alphas_sorted <- sort(ALPHAS)
lines <- c("\\begin{tabular}{lrrrrrr}", "\\toprule",
           "null family & $\\pi_0$ & $\\delta$ & $\\alpha=5\\times10^{-4}$ & $0.002$ & $0.01$ & $0.05$ \\\\",
           "\\midrule")
for (fam in FAMILIES) {
  for (j in seq_len(nrow(settings))) {
    x <- kappa[kappa$family == fam & kappa$pi0 == settings$pi0[j] &
               kappa$delta == settings$delta[j], ]
    x <- x[match(alphas_sorted, x$alpha), ]
    v <- ifelse(x$evaluable, formatC(x$kappa, format = "f", digits = 2), "---")
    lines <- c(lines, sprintf("%s & %.3f & %.2f & %s \\\\", family_label[fam],
                              settings$pi0[j], settings$delta[j],
                              paste(v, collapse = " & ")))
  }
  if (fam != FAMILIES[length(FAMILIES)]) lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(OUT_TABLES, "table_kappa.tex"))
message("  wrote ", file.path(OUT_TABLES, "table_kappa.tex"))
