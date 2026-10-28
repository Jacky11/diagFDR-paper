## =====================================================================
##  07 --- Reference scenarios reproduced with diagFDR::dfdr_simulate_tdc()
##
##  The simulation engine of sim_engine.R is now shipped in the package
##  (dfdr_simulate_tdc, dfdr_realised_fdp, dfdr_check_cv_benchmark). This
##  script reproduces the S0, S1, S2, S3 and S5 rows of Table S4 at
##  alpha = 1% (Gaussian nulls) and the kappa check of Note S1 from the
##  package alone. SIM_QUICK=1 runs 50 replicates instead of 1000.
## =====================================================================
suppressPackageStartupMessages(library(diagFDR))
B <- if (nzchar(Sys.getenv("SIM_QUICK"))) 50 else 1000
n <- 2e5; alpha <- 0.01
set.seed(20260811)

scen <- list(
  S0 = list(),
  S1 = list(decoy_shift = -0.35),
  S2 = list(decoy_scale = 0.85),
  S3 = list(compete = FALSE, pi0_lib = 0.45),
  `S5 (rho=0.5)` = list(rho = 0.5))

rows <- lapply(names(scen), function(s) {
  args <- c(list(B = B, alpha = alpha, n = n, pi0 = 0.667, delta = 4.65), scen[[s]])
  k <- do.call(dfdr_check_cv_benchmark, args)
  data.frame(scenario = s, FDP = 100 * k$mean_FDP, D_alpha = k$mean_D, kappa = k$kappa,
             share_FDP_above_alpha = k$share_FDP_above_alpha,
             Var_V = k$Var_V, two_ED_plus_VarD = k$two_ED_plus_VarD)
})
out <- do.call(rbind, rows)
print(out, digits = 3)
dir.create("out/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "out/tables/package_reproduction.csv", row.names = FALSE)
## Expected (B = 1000): S0 FDP 0.98, D 659, kappa 1.41; S1 2.59; S2 2.83; S3 0.45; S5 kappa 1.80
