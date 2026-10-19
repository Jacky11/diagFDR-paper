## Runs every simulation script in order.
##   Rscript run_all.R                 full run
##   SIM_QUICK=1 Rscript run_all.R     fast check of the pipeline
##   SIM_CORES=4 Rscript run_all.R     use four cores
for (script in c("01_reference_scenarios.R", "02_kappa_robustness.R",
                 "03_misspecification_grid.R", "04_library_enrichment.R",
                 "05_discrimination_vs_support.R", "06_proposition1_check.R")) {
  message("\n=== ", script, " ===")
  t0 <- Sys.time()
  source(script, local = new.env())
  message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
