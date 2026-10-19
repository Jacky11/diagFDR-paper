## =====================================================================
##  01 --- Reference scenarios S0 to S5
##
##  Runs every scenario at the reference parameters and writes
##    out/tables/reference_replicates.rds   one row per replicate x alpha
##    out/tables/reference_summary.csv      one row per cell x alpha
##    out/tables/table_flags.tex            flag rates at alpha = 1%
##
##  Usage:  Rscript 01_reference_scenarios.R
## =====================================================================

source("parameters.R")
cl <- make_cluster(N_CORES)

cells <- list()
add <- function(scenario, family, label = scenario, ...)
  cells[[length(cells) + 1]] <<- list(scenario = scenario, family = family,
                                      label = label, args = list(...))
for (fam in FAMILIES) {
  add("S0", fam); add("S1", fam); add("S2", fam); add("S3", fam)
  add("S4a", fam); add("S4b", fam)
}
for (rho in c(0.2, 0.5, 0.8))
  add("S5", "gaussian", label = sprintf("S5 (rho=%.1f)", rho), rho = rho)

results <- list()
for (i in seq_along(cells)) {
  cell <- cells[[i]]
  ## S4a and S4b share a seed: they are the same simulated runs, labelled
  ## by two different rules.
  seed_id <- if (cell$scenario == "S4b") i - 1 else i
  message(sprintf("[%2d/%d] %-8s %s", i, length(cells), cell$family, cell$label))
  d <- do.call(run_cell, c(list(B = B_REFERENCE, seed = SEED + 1e4 * seed_id,
                                scenario = cell$scenario, family = cell$family,
                                n = N_UNITS, pi0 = PI0, delta = DELTA,
                                alphas = ALPHAS, cl = cl), cell$args))
  d$scenario <- cell$label
  results[[i]] <- d
}
if (!is.null(cl)) parallel::stopCluster(cl)

replicates <- do.call(rbind, results)
saveRDS(replicates, file.path(OUT_TABLES, "reference_replicates.rds"))
summary_table <- summarise_cells(replicates)
write_table(summary_table, "reference_summary")

## ---------------------------------------------------------------------
##  LaTeX table: flag rates at the main threshold
## ---------------------------------------------------------------------
family_label <- c(gaussian = "Gaussian", gumbel = "Gumbel", t4 = "Student-$t_4$")
num <- function(x, digits = 2)
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))

s <- summary_table[summary_table$alpha == ALPHA_MAIN, ]
s <- s[order(match(s$scenario, unique(replicates$scenario)),
             match(s$family, FAMILIES)), ]
lines <- c(
  "\\begin{tabular}{llrrrrrrrrrr}", "\\toprule",
  paste("scenario & null family & FDP (\\%) & $D/T$ (\\%) & $\\overline{D_\\alpha}$ &",
        "$\\widehat\\pi_D$ & \\% bad &",
        "\\multicolumn{4}{c}{share of replicates flagged} & sens. \\\\"),
  " & & & & & & & $D$ & $D^{-1/2}$ & IPE & eq. & any \\\\", "\\midrule",
  sprintf("%s & %s & %s & %s & %.0f & %s & %.0f & %s & %s & %s & %s & %s \\\\",
          s$scenario, family_label[s$family],
          num(100 * s$mean_fdp), num(100 * s$mean_decoy_ratio), s$mean_D,
          num(s$mean_pi_D, 3), 100 * s$frac_bad,
          num(s$fire_D), num(s$fire_CV),
          ifelse(is.nan(s$mean_ipe_tdc), "n/e", num(s$fire_IPE)),
          num(s$fire_EQ), num(s$sens_any)),
  "\\bottomrule", "\\end{tabular}")
lines <- gsub("rho=", "$\\\\rho=$", lines)
writeLines(lines, file.path(OUT_TABLES, "table_flags.tex"))
message("  wrote ", file.path(OUT_TABLES, "table_flags.tex"))
