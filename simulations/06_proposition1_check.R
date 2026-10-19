## =====================================================================
##  06 --- What Proposition 1 approximates (Note S1)
##
##  Proposition 1 gives sqrt(2) * D^(-1/2) as a benchmark for the
##  coefficient of variation of the realised FDP. Its idealised stopping
##  representation gives E[V] ~ E[D] and Var(V) ~ 2 E[D] + Var(D)
##  (law of total variance), so CV(FDP) ~ sqrt(2) D^(-1/2) only when D_alpha
##  and T_alpha vary little between experiments. This script checks, under
##  the correct specification S0 for the three null families:
##    (1) the marginal agreement: Var(V) against 2 E[D] + Var(D) and 2 E[D],
##        and kappa = CV(FDP) / mean(D^(-1/2)) against sqrt(2);
##    (2) the conditional law: among replicates with nearly the same D_alpha,
##        is V ~ NegBin(D, 1/2) (mean D, variance 2D)?
##    (3) assumption (iii): the share of accepted decoys that come from
##        units carrying a correct target (a decoy that out-scores it).
##  Expected (Gaussian, alpha = 1%): kappa ~ 1.4, Var(V) ~ 2E[D]; the
##  conditional variance is 63-86% of 2D and E[V | D] rises 4-12 per decoy;
##  about 2% of accepted decoys come from correct units (about 26% with
##  Gumbel and 18% with Student-t4 nulls, where kappa is 1.7-1.8).
##
##  Usage:  Rscript 06_proposition1_check.R      (SIM_QUICK=1 for a dry run)
##  Output: out/tables/proposition1_check.csv
##          out/tables/proposition1_decoy_share.csv
## =====================================================================

source("parameters.R")
B <- if (QUICK) 50 else 1000
cl <- make_cluster(N_CORES)

summ <- function(d, family) do.call(rbind, lapply(ALPHAS, function(a) {
  x <- d[d$alpha == a & d$D > 0 & d$T > 0, ]
  if (nrow(x) < 20 || mean(x$D) < 10) return(NULL)
  fdp <- x$V / x$T
  w <- max(3, 0.025 * mean(x$D))
  fit <- stats::lm(V ~ D, data = x)
  near <- abs(x$D - mean(x$D)) <= w
  data.frame(
    family = family, alpha = a, B = nrow(x), mean_D = mean(x$D),
    mean_V = mean(x$V), var_V = var(x$V), two_mean_D = 2 * mean(x$D),
    two_mean_D_plus_var_D = 2 * mean(x$D) + var(x$D),
    kappa = (sd(fdp) / mean(fdp)) / mean(x$D^-0.5),
    cv_fdp = sd(fdp) / mean(fdp), bench = sqrt(2) * mean(x$D^-0.5),
    slope_V_on_D = unname(coef(fit)[2]),
    var_V_given_D = mean(residuals(fit)[near]^2),
    ratio_var_given_D_to_2D = mean(residuals(fit)[near]^2) / (2 * mean(x$D[near])))
}))

res <- do.call(rbind, lapply(seq_along(FAMILIES), function(k) {
  d <- run_cell(B = B, seed = SEED + 6e5 + 1e4 * (k - 1), scenario = "S0",
                family = FAMILIES[k], n = N_UNITS, pi0 = PI0, delta = DELTA,
                alphas = ALPHAS, cl = cl)
  summ(d, FAMILIES[k])
}))
if (!is.null(cl)) parallel::stopCluster(cl)
print(res, digits = 3)
write_table(res, "proposition1_check")

## (3) share of accepted decoys at alpha = 1% coming from correct units
B3 <- if (QUICK) 10 else 200
share <- do.call(rbind, lapply(FAMILIES, function(fm) {
  r <- vapply(seq_len(B3), function(b) {
    set.seed(SEED + 7e5 + b)
    is_null <- runif(N_UNITS) < PI0
    X <- draw_null(N_UNITS, fm) + ifelse(is_null, 0, DELTA)
    Y <- draw_null(N_UNITS, fm)
    tg <- X > Y; q <- tdc_qvalues(pmax(X, Y), tg)
    dec <- q <= ALPHA_MAIN & !tg
    c(sum(dec & !is_null) / max(sum(dec), 1), sum(dec))
  }, numeric(2))
  data.frame(family = fm, alpha = ALPHA_MAIN, B = B3,
             share_decoys_from_correct_units = mean(r[1, ]), mean_D = mean(r[2, ]))
}))
print(share, digits = 3)
write_table(share, "proposition1_decoy_share")
