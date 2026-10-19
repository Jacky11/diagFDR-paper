## =====================================================================
##  03 --- Decoy misspecification over a range of magnitudes
##
##  Scenarios S1 (shifted decoys) and S2 (narrowed decoys) at five
##  magnitudes each, for the three null families, at alpha = 1%.
##
##  Writes out/tables/misspecification_grid.csv and table_grid.tex
## =====================================================================

source("parameters.R")
cl <- make_cluster(N_CORES)

SHIFTS <- c(-0.10, -0.20, -0.35, -0.50, -0.70)
SCALES <- c(0.95, 0.90, 0.85, 0.80, 0.70)

rows <- list(); k <- 0
run_one <- function(scenario, fam, magnitude, ...) {
  k <<- k + 1
  message(sprintf("[%2d/30] %-8s %s %+.2f", k, fam, scenario, magnitude))
  d <- run_cell(B_GRID, SEED + 4e6 + 1e4 * k, scenario, fam, N_UNITS, PI0,
                DELTA, ALPHA_MAIN, ..., cl = cl)
  s <- summarise_cells(d)
  rows[[k]] <<- cbind(magnitude = magnitude, s)
}
for (fam in FAMILIES) for (s in SHIFTS) run_one("S1", fam, s, shift = s)
for (fam in FAMILIES) for (s in SCALES) run_one("S2", fam, s, scale = s)
if (!is.null(cl)) parallel::stopCluster(cl)

grid <- do.call(rbind, rows)
grid$fold_nominal <- grid$mean_fdp / ALPHA_MAIN
write_table(grid, "misspecification_grid")

family_label <- c(gaussian = "Gaussian", gumbel = "Gumbel", t4 = "Student-$t_4$")
num <- function(x) ifelse(is.na(x), "---", formatC(x, format = "f", digits = 2))
block <- function(scenario, symbol, title) {
  g <- grid[grid$scenario == scenario, ]
  c(sprintf("\\multicolumn{8}{l}{\\textbf{%s}} \\\\", title),
    sprintf("Null family & %s & realised FDP & $\\times$ nominal & $\\overline{D_\\alpha}$ & $\\widehat\\pi_D$ & \\%% bad & sens.\\ any flag \\\\", symbol),
    "\\midrule",
    sprintf("%s & $%.2f$ & $%.2f\\%%$ & $%.2f$ & $%.0f$ & $%.3f$ & $%.0f$ & %s \\\\",
            family_label[g$family], g$magnitude, 100 * g$mean_fdp,
            g$fold_nominal, g$mean_D, g$mean_pi_D, 100 * g$frac_bad,
            num(g$sens_any)))
}
lines <- c("\\begin{tabular}{llrrrrrr}", "\\toprule",
           block("S1", "$s$", "S1: decoy scores shifted by $s$"), "\\midrule",
           block("S2", "$c$", "S2: decoy scores multiplied by $c$"),
           "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(OUT_TABLES, "table_grid.tex"))
message("  wrote ", file.path(OUT_TABLES, "table_grid.tex"))
