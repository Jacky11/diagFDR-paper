## =====================================================================
##  04 --- Scenario S3 over a range of library compositions
##
##  Targets and decoys are listed separately and the null fraction among
##  targets, pi0_library, is varied. With no failure of the decoy model,
##  the decoy fraction among low-confidence entries moves away from 0.5
##  and the decoy ratio over-states the realised FDP by about
##  1 / pi0_library.
##
##  Writes out/tables/library_enrichment.csv and table_enrichment.tex
## =====================================================================

source("parameters.R")
cl <- make_cluster(N_CORES)

PI0_LIBRARY <- c(0.90, 0.80, 0.667, 0.50, 0.45, 0.40, 0.30)

rows <- list()
for (k in seq_along(PI0_LIBRARY)) {
  message(sprintf("[%d/%d] pi0_library = %.3f", k, length(PI0_LIBRARY), PI0_LIBRARY[k]))
  d <- run_cell(B_GRID, SEED + 6e6 + 1e4 * k, "S3", "gaussian", N_UNITS, PI0,
                DELTA, ALPHA_MAIN, pi0_library = PI0_LIBRARY[k], cl = cl)
  rows[[k]] <- cbind(pi0_library = PI0_LIBRARY[k], summarise_cells(d))
}
if (!is.null(cl)) parallel::stopCluster(cl)

enrich <- do.call(rbind, rows)
write_table(enrich, "library_enrichment")

lines <- c("\\begin{tabular}{rrrrrr}", "\\toprule",
           "$\\pi_0^{\\mathrm{lib}}$ & realised FDP & decoy ratio & ratio / FDP & $\\widehat\\pi_D$ & eq.-chance flag \\\\",
           "\\midrule",
           sprintf("$%.3f$ & $%.2f\\%%$ & $%.2f\\%%$ & $%.2f$ & $%.3f$ & $%.2f$ \\\\",
                   enrich$pi0_library, 100 * enrich$mean_fdp,
                   100 * enrich$mean_decoy_ratio,
                   enrich$mean_decoy_ratio / enrich$mean_fdp,
                   enrich$mean_pi_D, enrich$fire_EQ),
           "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(OUT_TABLES, "table_enrichment.tex"))
message("  wrote ", file.path(OUT_TABLES, "table_enrichment.tex"))
