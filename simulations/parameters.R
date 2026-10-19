## =====================================================================
##  Parameters shared by all simulation scripts
## =====================================================================

N_UNITS <- 2e5                              # units per simulated list
## pi0 and delta are DESIGN VALUES (Note S6.2). They were first derived from
## the DIA-NN run x precursor export (T_tot = 680,203, D_tot = 340,262), but
## both derivations are set by bounds rather than by the data:
##   pi0   = 2 D_tot / (T_tot + D_tot) = 2 q_max / (1 + q_max) = 0.667 for an
##           export capped at q_max = 0.5 (--qvalue 0.5), whatever the data;
##   delta = qnorm(frac) + qnorm(0.99) with frac = T_0.01 / (T_tot (1 - pi0))
##           = 1.656 > 1, clipped to 0.99, hence delta = 2 * qnorm(0.99) = 4.65:
##           99% of correct targets exceed the 1% upper quantile of the null.
## With these values S0 gives E[D_0.01] = 659, of the order of the 832 decoys
## of the DIA-NN global list. 02_kappa_robustness.R varies (pi0, delta).
PI0     <- 0.667                            # null fraction (design value)
DELTA   <- 4.65                             # separation of correct targets (design value)
ALPHAS  <- c(0.05, 0.01, 0.002, 5e-4)       # nominal FDR levels
ALPHA_MAIN <- 0.01
FAMILIES   <- c("gaussian", "gumbel", "t4")

## Replicates. Set the environment variable SIM_QUICK=1 for a fast check
## of the whole pipeline (results are then too noisy to report).
QUICK <- nzchar(Sys.getenv("SIM_QUICK"))
B_REFERENCE <- if (QUICK) 20 else 1000      # reference scenarios S0-S5
B_GRID      <- if (QUICK) 10 else 200       # misspecification grid
B_KAPPA     <- if (QUICK) 10 else 500       # robustness of kappa
B_SIGNFLIP  <- if (QUICK) 5  else 200       # discrimination vs support

N_CORES <- as.integer(Sys.getenv("SIM_CORES", "1"))
SEED    <- 20260101L

OUT_TABLES  <- file.path("out", "tables")
OUT_FIGURES <- file.path("out", "figures")
dir.create(OUT_TABLES,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGURES, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "sim_engine.R"))

write_table <- function(x, name) {
  path <- file.path(OUT_TABLES, paste0(name, ".csv"))
  write.csv(x, path, row.names = FALSE)
  message("  wrote ", path)
}
