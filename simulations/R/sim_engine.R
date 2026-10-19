## =====================================================================
##  Ground-truth simulation engine for target-decoy FDR diagnostics
##
##  Every simulated list records which targets are truly incorrect, so the
##  realised false discovery proportion (FDP) of any accepted list is known
##  exactly and can be compared with what the list itself reports.
##
##  Vocabulary used throughout
##    unit       one tested hypothesis (e.g. a precursor)
##    null unit  a unit whose target is an incorrect match
##    T, D       accepted targets and accepted decoys at a threshold alpha
##    V          accepted targets that are null (unobservable in real data)
##    FDP        V / T, the realised false discovery proportion
##    decoy ratio  D / T, the quantity a target-decoy list reports
## =====================================================================

## ---------------------------------------------------------------------
##  Null score distributions, all standardised to mean 0 and variance 1
##  so that the separation `delta` means the same thing in each family.
## ---------------------------------------------------------------------
EULER_GAMMA <- 0.5772156649015329

draw_null <- function(n, family) {
  switch(family,
    gaussian = rnorm(n),
    gumbel   = (-log(-log(runif(n))) - EULER_GAMMA) / (pi / sqrt(6)),
    t4       = rt(n, df = 4) / sqrt(2),
    stop("unknown null family: ", family))
}

## ---------------------------------------------------------------------
##  Target-decoy q-values.
##  Entries are ranked by decreasing score; the estimated FDR at each rank
##  is (#decoys) / (#targets) above it, and the q-value is the running
##  minimum of that estimate taken from the bottom of the list.
## ---------------------------------------------------------------------
tdc_qvalues <- function(score, is_target) {
  o   <- order(score, decreasing = TRUE)
  tgt <- is_target[o]
  fdr <- cumsum(!tgt) / pmax(cumsum(tgt), 1)
  q   <- numeric(length(score))
  q[o] <- rev(cummin(rev(fdr)))
  q
}

## ---------------------------------------------------------------------
##  Internal PEP error (IPE): size-weighted mean absolute difference
##  between the mean PEP of targets in a PEP bin and a decoy-derived
##  reference for that bin. Two references are computed:
##    "tdc"    D_b / T_b            (the one used in the paper)
##    "pooled" D_b / (T_b + D_b)
##  Only bins with mean PEP <= 0.5 and at least `min_bin` entries count.
## ---------------------------------------------------------------------
ipe <- function(pep, is_target, reference = c("tdc", "pooled"),
                n_bins = 10, min_bin = 200) {
  reference <- match.arg(reference)
  bin <- cut(pep, seq(0, 0.5, length.out = n_bins + 1), include.lowest = TRUE)
  num <- 0; den <- 0
  for (ix in split(seq_along(pep), bin)) {
    if (length(ix) < min_bin) next
    Tb <- sum(is_target[ix]); Db <- length(ix) - Tb
    if (Tb == 0) next
    ref <- if (reference == "tdc") min(1, Db / Tb) else Db / (Tb + Db)
    num <- num + length(ix) * abs(mean(pep[ix][is_target[ix]]) - ref)
    den <- den + length(ix)
  }
  if (den > 0) num / den else NA_real_
}

## ---------------------------------------------------------------------
##  Summarise one simulated list at each threshold in `alphas`.
##    q          q-value of each entry
##    is_target  TRUE for targets, FALSE for decoys
##    is_null    TRUE where a target is an incorrect match
##    pep        optional true posterior error probabilities
##  The equal-chance statistic pi_D is the decoy fraction among entries
##  with q in `eq_band`; it does not depend on alpha.
## ---------------------------------------------------------------------
summarise_list <- function(q, is_target, is_null, alphas, pep = NULL,
                           eq_band = c(0.2, 0.5)) {
  in_band <- q >= eq_band[1] & q <= eq_band[2]
  n_band  <- sum(in_band)
  pi_D    <- if (n_band > 0) mean(!is_target[in_band]) else NA_real_
  ipe_tdc <- ipe_pool <- NA_real_
  if (!is.null(pep)) {
    ipe_tdc  <- ipe(pep, is_target, "tdc")
    ipe_pool <- ipe(pep, is_target, "pooled")
  }
  out <- t(vapply(alphas, function(a) {
    acc <- q <= a
    Ta  <- sum(acc & is_target)
    Da  <- sum(acc) - Ta
    Va  <- sum(acc & is_target & is_null)
    c(alpha = a, T = Ta, D = Da, V = Va,
      fdp = if (Ta > 0) Va / Ta else NA_real_,
      decoy_ratio = if (Ta > 0) Da / Ta else NA_real_)
  }, numeric(6)))
  cbind(out, pi_D = pi_D, n_band = n_band,
        ipe_tdc = ipe_tdc, ipe_pooled = ipe_pool)
}

## ---------------------------------------------------------------------
##  One replicate of one scenario. Returns one row per alpha.
##
##  Common to all scenarios: n units; a unit is null with probability pi0;
##  the target of a non-null unit scores `delta` higher than a null score.
##
##  S0   correct specification. Each unit has one target and one decoy
##       score from the same null distribution; the higher one is kept
##       (paired competition), so a null unit yields a target or a decoy
##       with probability 1/2.
##  S1   as S0, with decoy scores shifted by `shift` (< 0: decoys are
##       easier to beat than incorrect targets).
##  S2   as S0, with decoy scores multiplied by `scale` (< 1: the upper
##       tail of the decoys is depleted).
##  S3   enriched library without competition. Targets and decoys are
##       listed separately (no pairing), and the null fraction among
##       targets is `pi0_library` instead of pi0. This is what a second
##       search pass produces when its library has been reduced to the
##       identifications of a first pass. Nothing is wrong with the
##       decoys: null targets and decoys have the same distribution.
##  S4a, S4b  cross-run aggregation. `n_runs` independent runs of S0 are
##       searched separately; each unit keeps the smallest of its run-level
##       q-values. Its label is "decoy if a decoy won in any run" (S4a) or
##       "the label of the run that supplied the smallest q-value" (S4b).
##  S5   as S0 with dependence between units: scores share a block effect
##       within blocks of `block_size` units (within-block correlation
##       `rho`). Target and decoy block effects are independent, so equal
##       chance still holds for every null unit. Gaussian family only.
## ---------------------------------------------------------------------
sim_one <- function(scenario, family, n, pi0, delta, alphas,
                    shift = -0.35, scale = 0.85, pi0_library = 0.45,
                    n_runs = 6, rho = 0.5, block_size = 10) {

  competed <- function(X, Y) list(score = pmax(X, Y), is_target = X > Y)

  if (scenario %in% c("S0", "S1", "S2")) {
    is_null <- runif(n) < pi0
    X <- draw_null(n, family) + ifelse(is_null, 0, delta)
    Y <- draw_null(n, family)
    if (scenario == "S1") Y <- Y + shift
    if (scenario == "S2") Y <- Y * scale
    L <- competed(X, Y)
    q <- tdc_qvalues(L$score, L$is_target)
    pep <- NULL
    if (family == "gaussian") {
      ## True PEP of a target that won with score w under S0. The factor
      ## pnorm(w), the chance of beating the decoy, cancels out.
      f0 <- dnorm(L$score); f1 <- dnorm(L$score - delta)
      pep <- pi0 * f0 / (pi0 * f0 + (1 - pi0) * f1)
    }
    return(summarise_list(q, L$is_target, is_null, alphas, pep))
  }

  if (scenario == "S3") {
    is_null <- runif(n) < pi0_library
    X <- draw_null(n, family) + ifelse(is_null, 0, delta)
    Y <- draw_null(n, family)
    score     <- c(X, Y)
    is_target <- rep(c(TRUE, FALSE), each = n)
    q <- tdc_qvalues(score, is_target)
    return(summarise_list(q, is_target, c(is_null, rep(FALSE, n)), alphas))
  }

  if (scenario %in% c("S4a", "S4b")) {
    is_null <- runif(n) < pi0
    q_run   <- matrix(NA_real_, n, n_runs)
    tgt_run <- matrix(NA, n, n_runs)
    for (r in seq_len(n_runs)) {
      X <- draw_null(n, family) + ifelse(is_null, 0, delta)
      Y <- draw_null(n, family)
      L <- competed(X, Y)
      q_run[, r]   <- tdc_qvalues(L$score, L$is_target)
      tgt_run[, r] <- L$is_target
    }
    best <- max.col(-q_run, ties.method = "first")
    q    <- q_run[cbind(seq_len(n), best)]
    is_target <- if (scenario == "S4a") rowSums(tgt_run) == n_runs
                 else tgt_run[cbind(seq_len(n), best)]
    return(summarise_list(q, is_target, is_null, alphas))
  }

  if (scenario == "S5") {
    if (family != "gaussian") stop("S5 is defined for the gaussian family")
    n_blocks <- ceiling(n / block_size)
    block    <- rep(seq_len(n_blocks), each = block_size)[seq_len(n)]
    is_null  <- runif(n) < pi0
    X <- sqrt(rho) * rnorm(n_blocks)[block] + sqrt(1 - rho) * rnorm(n) +
         ifelse(is_null, 0, delta)
    Y <- sqrt(rho) * rnorm(n_blocks)[block] + sqrt(1 - rho) * rnorm(n)
    L <- competed(X, Y)
    q <- tdc_qvalues(L$score, L$is_target)
    return(summarise_list(q, L$is_target, is_null, alphas))
  }

  stop("unknown scenario: ", scenario)
}

## ---------------------------------------------------------------------
##  Run B replicates of one cell. Replicate b uses seed `seed + b`, so a
##  cell gives the same result whether it is run serially or in parallel.
##  `cl` is an optional cluster created with parallel::makeCluster().
## ---------------------------------------------------------------------
run_cell <- function(B, seed, scenario, family, n, pi0, delta, alphas,
                     ..., cl = NULL) {
  args <- c(list(scenario, family, n, pi0, delta, alphas), list(...))
  seed <- as.integer(seed)
  one <- function(b) {
    set.seed(seed + b)
    cbind(rep = b, do.call(sim_one, args))
  }
  reps <- if (is.null(cl)) lapply(seq_len(B), one)
          else parallel::parLapply(cl, seq_len(B), one)
  d <- as.data.frame(do.call(rbind, reps))
  d$scenario <- scenario; d$family <- family
  d
}

make_cluster <- function(n_cores) {
  if (n_cores <= 1) return(NULL)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("draw_null", "tdc_qvalues", "ipe",
                                "summarise_list", "sim_one", "EULER_GAMMA"),
                          envir = globalenv())
  cl
}

## ---------------------------------------------------------------------
##  Flags evaluated on each replicate (thresholds as in the paper).
##    D     fewer than 10 accepted decoys
##    CV    decoy-support index D^(-1/2) above 0.30
##    IPE   internal PEP error above 0.05 (where a true PEP is available)
##    EQ    decoy fraction in the equal-chance band more than 0.15 away
##          from 0.5 (needs at least `min_band` entries in the band)
## ---------------------------------------------------------------------
FLAG_THRESHOLDS <- list(D = 10, CV = 0.30, IPE = 0.05, EQ = 0.15, min_band = 100)
BAD_FACTOR <- 1.5      # a replicate is "bad" when FDP > BAD_FACTOR * alpha

add_flags <- function(d, thr = FLAG_THRESHOLDS) {
  d$flag_D   <- d$D < thr$D
  d$flag_CV  <- d$D == 0 | 1 / sqrt(pmax(d$D, 1)) > thr$CV
  d$flag_IPE <- !is.na(d$ipe_tdc) & d$ipe_tdc > thr$IPE
  d$eq_evaluable <- d$n_band >= thr$min_band
  d$flag_EQ  <- d$eq_evaluable & abs(d$pi_D - 0.5) > thr$EQ
  d$flag_any <- d$flag_D | d$flag_CV | d$flag_IPE | d$flag_EQ
  d$bad      <- !is.na(d$fdp) & d$fdp > BAD_FACTOR * d$alpha
  d
}

## ---------------------------------------------------------------------
##  One summary row per (group columns, alpha).
##  kappa is the ratio of the observed coefficient of variation of the FDP
##  to the mean decoy-support index D^(-1/2); Proposition 1 predicts
##  sqrt(2) under independence.
## ---------------------------------------------------------------------
summarise_cells <- function(d, by = c("family", "scenario")) {
  d <- add_flags(d)
  key <- interaction(d[c(by, "alpha")], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(d, key), function(x) {
    ok   <- !is.na(x$fdp) & x$D > 0
    sens <- function(f) if (any(x$bad)) mean(f[x$bad]) else NA_real_
    cv   <- sd(x$fdp[ok]) / mean(x$fdp[ok])
    dsi  <- mean(x$D[ok]^-0.5)
    data.frame(
      x[1, by, drop = FALSE], alpha = x$alpha[1], B = nrow(x),
      mean_T = mean(x$T), mean_D = mean(x$D),
      mean_V = mean(x$V), var_V = var(x$V),
      mean_fdp = mean(x$fdp, na.rm = TRUE),
      mean_decoy_ratio = mean(x$decoy_ratio, na.rm = TRUE),
      p_fdp_above_alpha = mean(x$fdp > x$alpha, na.rm = TRUE),
      frac_bad = mean(x$bad), n_bad = sum(x$bad),
      cv_fdp = cv, mean_support_index = dsi, kappa = cv / dsi,
      frac_with_decoys = mean(ok),
      mean_pi_D = mean(x$pi_D, na.rm = TRUE),
      mean_n_band = mean(x$n_band),
      frac_eq_evaluable = mean(x$eq_evaluable),
      mean_ipe_tdc = mean(x$ipe_tdc, na.rm = TRUE),
      mean_ipe_pooled = mean(x$ipe_pooled, na.rm = TRUE),
      fire_D = mean(x$flag_D), fire_CV = mean(x$flag_CV),
      fire_IPE = mean(x$flag_IPE), fire_EQ = mean(x$flag_EQ),
      fire_any = mean(x$flag_any),
      sens_D = sens(x$flag_D), sens_CV = sens(x$flag_CV),
      sens_IPE = sens(x$flag_IPE), sens_EQ = sens(x$flag_EQ),
      sens_any = sens(x$flag_any),
      row.names = NULL)
  })
  do.call(rbind, rows)
}
