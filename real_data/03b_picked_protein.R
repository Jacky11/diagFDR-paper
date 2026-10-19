## =====================================================================
##  03b --- Picked-protein competition at protein-group level
##
##  A procedural comparator for the four label-propagation rules.
##
##  WHY. The four rules of Note S9 all assume a protein group must
##  INHERIT a label from its constituent rows. Picked-protein FDR
##  (Savitski et al. 2015; The et al. 2021) declines the premise: the
##  target and decoy versions of the SAME protein compete directly, the
##  better of the pair is kept and the other discarded, and q-values are
##  then recomputed by target-decoy counting over the retained entries.
##  Because every retained entry comes from a matched pair, the null is
##  explicit and no label has to be propagated.
##
##  WHY IT WORKS ON THIS EXPORT. DIA-NN writes the SAME Protein.Group
##  accession string for the target and the decoy precursors of a group,
##  which is exactly what makes the pairing exact. That same fact is what
##  makes the four label rules disagree in the first place.
##
##  ONE LIMITATION. The competition uses
##  Global.PG.Q.Value as the evidence measure rather than a raw score,
##  because this export carries no protein-level score column. That is a
##  monotone transform of the underlying evidence, so the ranking is
##  right, but the pairing inherits whatever protein inference DIA-NN
##  performed upstream.
##
##  RECOMPUTATION IS OVER THE FULL UNIVERSE, not the 1%-filtered export:
##  the decoy tail needed for the q-value recomputation is only present
##  in the unfiltered table.
##
##  Expected on the Zenodo data (basic search, alpha = 1%):
##    identifiers 13,751; both labels present 10,262
##    picked list 13,751 (12,022 target, 1,729 decoy)
##    T = 9,146   D = 91   decoy ratio = 0.995%
##  and on the entrapment search:
##    T = 8,825   D = 88   ratio 0.997%   entrap FDP 0.771% [0.608, 0.976]
##
##  Outputs, in tables/:
##    picked_protein_basic.csv     the Panel B row
##    picked_protein_entrapment.csv
##    label_rules_vs_entrapment.csv   label rules on the entrapment search
## =====================================================================

source("00_config.R")
load_numbers()

suppressPackageStartupMessages({
  library(dplyr)
})

## ---------------------------------------------------------------------
##  Core: picked competition on a (id, is_decoy, q) table
## ---------------------------------------------------------------------

##' Picked-protein competition.
##'
##' @param df data.frame with columns id (character), is_decoy (logical),
##'   q (numeric). One row per (id, is_decoy) after taking the group
##'   minimum; the function will take the minimum itself if given more.
##' @param offset added to the decoy count in the FDR estimate; 0
##'   reproduces the paper, 1 gives the (D+1)/T convention.
##' @return list(tbl, T_alpha, D_alpha, ratio, n_pairs, n_ids)
picked_protein <- function(df, alpha = ALPHA_MAIN, offset = 0) {

  stopifnot(all(c("id", "is_decoy", "q") %in% names(df)))
  df <- df[!is.na(df$q) & !is.na(df$is_decoy), , drop = FALSE]

  ## group minimum per (id, label)
  best <- df %>%
    group_by(id, is_decoy) %>%
    summarise(q = min(q), .groups = "drop")

  wide <- best %>%
    tidyr::pivot_wider(names_from = is_decoy, values_from = q,
                       names_prefix = "q_") %>%
    rename(qT = q_FALSE, qD = q_TRUE)
  if (!"qT" %in% names(wide)) wide$qT <- NA_real_
  if (!"qD" %in% names(wide)) wide$qD <- NA_real_

  n_pairs <- sum(!is.na(wide$qT) & !is.na(wide$qD))

  ## the pick: whichever of the pair has the SMALLER q-value survives
  pick_is_decoy <- ifelse(is.na(wide$qT), TRUE,
                   ifelse(is.na(wide$qD), FALSE, wide$qD < wide$qT))
  pick_q <- ifelse(pick_is_decoy,
                   ifelse(is.na(wide$qD), Inf, wide$qD),
                   ifelse(is.na(wide$qT), Inf, wide$qT))
  keep <- is.finite(pick_q)
  P <- data.frame(id = wide$id[keep], is_decoy = pick_is_decoy[keep],
                  q_in = pick_q[keep], stringsAsFactors = FALSE)

  ## recompute q-values by target-decoy counting on the picked list.
  ## Stable sort so that ties never place a target above an equal decoy.
  o  <- order(P$q_in, method = "radix")
  P  <- P[o, , drop = FALSE]
  Tc <- cumsum(!P$is_decoy)
  Dc <- cumsum(P$is_decoy)
  fdr <- (Dc + offset) / pmax(Tc, 1)
  P$q_picked <- rev(cummin(rev(fdr)))

  acc <- P$q_picked <= alpha
  Ta <- sum(acc & !P$is_decoy)
  Da <- sum(acc &  P$is_decoy)

  list(tbl = P, T_alpha = Ta, D_alpha = Da,
       ratio = if (Ta > 0) Da / Ta else NA_real_,
       cv = if (Da > 0) sqrt(2) / sqrt(Da) else NA_real_,
       n_pairs = n_pairs, n_ids = nrow(wide),
       n_pick_T = sum(!P$is_decoy), n_pick_D = sum(P$is_decoy),
       accepted_targets = P$id[acc & !P$is_decoy])
}

## ---------------------------------------------------------------------
##  A --- basic search: the Panel B row of main-text Table 2
## ---------------------------------------------------------------------
message("\n=== A. picked-protein competition, basic search ===")

rep_basic <- sanitise_report(diagFDR::read_diann_parquet(DIANN_BASIC),
                             COL$decoy, "basic")

pg_basic <- data.frame(
  id       = as.character(rep_basic[[COL$protein]]),
  is_decoy = as.logical(rep_basic[[COL$decoy]] == 1L),
  q        = as.numeric(rep_basic[[COL$q_pg_glob]]),
  stringsAsFactors = FALSE)

pk <- picked_protein(pg_basic)

say("identifiers %s; carrying both labels %s",
    format(pk$n_ids, big.mark = ","), format(pk$n_pairs, big.mark = ","))
say("picked list %s (target %s, decoy %s)",
    format(pk$n_pick_T + pk$n_pick_D, big.mark = ","),
    format(pk$n_pick_T, big.mark = ","), format(pk$n_pick_D, big.mark = ","))
say("at alpha=%.3g: T=%s D=%s ratio=%.3f%%  sqrt(2)/sqrt(D)=%.4f",
    ALPHA_MAIN, format(pk$T_alpha, big.mark = ","), pk$D_alpha,
    100 * pk$ratio, pk$cv)

reg("pickedT",     pk$T_alpha, "int", "Table 2 Panel B: picked-protein")
reg("pickedD",     pk$D_alpha, "int", "Table 2 Panel B: picked-protein")
reg("pickedRatio", pk$ratio,   "pct2","Table 2 Panel B: picked-protein")
reg("pickedCV",    pk$cv,      "d4",  "Table 2 Panel B: sqrt(2) D^-1/2")

write_tab(data.frame(rule = "picked_protein",
                     T_alpha = pk$T_alpha, D_alpha = pk$D_alpha,
                     ratio = pk$ratio, cv = pk$cv,
                     n_ids = pk$n_ids, n_pairs = pk$n_pairs),
          "picked_protein_basic")

## --- the four label rules on the 1%-filtered protein list -------------
## Table 2 Panel B and Note S9.3: rows with Global.PG.Q.Value <= alpha are
## collapsed to one unit per Protein.Group and the unit is labelled by
## each rule. `split` keeps the target and decoy entries apart.
label_rules_filtered <- function(pg, alpha = ALPHA_MAIN) {
  f <- pg[!is.na(pg$q) & !is.na(pg$is_decoy) & pg$q <= alpha, , drop = FALSE]
  agg <- f %>%
    group_by(id) %>%
    summarise(qmin = min(q), n = n(), nD = sum(is_decoy),
              win_is_D = is_decoy[which.min(q)], .groups = "drop")
  lab <- list(any_decoy   = agg$nD > 0,
              winning_row = agg$win_is_D,
              majority    = agg$nD / agg$n > 0.5)
  out <- do.call(rbind, lapply(names(lab), function(nm) {
    D <- sum(lab[[nm]]); Tn <- sum(!lab[[nm]])
    data.frame(rule = nm, units = nrow(agg), T_alpha = Tn, D_alpha = D,
               ratio = D / Tn, support_index = 1 / sqrt(D),
               fdp_cv_benchmark = sqrt(2 / D), stringsAsFactors = FALSE)
  }))
  sp <- f %>% group_by(id, is_decoy) %>% summarise(q = min(q), .groups = "drop")
  D <- sum(sp$is_decoy); Tn <- sum(!sp$is_decoy)
  out <- rbind(out, data.frame(rule = "split", units = nrow(sp),
               T_alpha = Tn, D_alpha = D, ratio = D / Tn,
               support_index = 1 / sqrt(D), fdp_cv_benchmark = sqrt(2 / D),
               stringsAsFactors = FALSE))
  attr(out, "n_groups") <- nrow(agg)
  attr(out, "n_mixed")  <- sum(agg$nD > 0 & agg$nD < agg$n)
  out
}

message("\n=== A2. label rules on the 1%-filtered protein list, basic search ===")
lr <- label_rules_filtered(pg_basic)
print(lr)
say("groups accepted: %s; retaining both a target and a decoy row: %s",
    format(attr(lr, "n_groups"), big.mark = ","), attr(lr, "n_mixed"))
write_tab(lr, "table2_panelB_filtered")

## Sanity: picked must land near nominal on a correctly scoped export.
## If it does not, the pairing is wrong -- check that Protein.Group really
## carries the same string for target and decoy rows in your DIA-NN
## version, because that is the assumption the whole procedure rests on.
if (is.finite(pk$ratio) && abs(pk$ratio - ALPHA_MAIN) > 0.5 * ALPHA_MAIN)
  warning("Picked-protein ratio is ", round(100 * pk$ratio, 3),
          "% against a nominal ", 100 * ALPHA_MAIN,
          "%. Check the target/decoy pairing on Protein.Group.",
          call. = FALSE)

## ---------------------------------------------------------------------
##  B --- entrapment search: label rules and picked, side by side,
##        with the external error rate of each accepted TARGET list.
##        This table is the
##        result that separates "changes what is reported" from
##        "changes what is accepted".
## ---------------------------------------------------------------------
message("\n=== B. label rules vs picked, with entrapment (Table S) ===")

if (!file.exists(DIANN_ENTRAP)) {
  message("  no entrapment report at ", DIANN_ENTRAP, " -- part B skipped.")
} else {
  need_file(FASTA_ENTRAP, "entrapment FASTA")
  ACC <- entrap_accessions(FASTA_ENTRAP)

  is_pure <- function(pg) {
    vapply(strsplit(ifelse(is.na(pg), "", pg), ";"), function(z) {
      z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) }, logical(1))
  }

  rep_ent <- sanitise_report(diagFDR::read_diann_parquet(DIANN_ENTRAP),
                             COL$decoy, "entrapment")
  pg_ent <- data.frame(
    id       = as.character(rep_ent[[COL$protein]]),
    is_decoy = as.logical(rep_ent[[COL$decoy]] == 1L),
    q        = as.numeric(rep_ent[[COL$q_pg_glob]]),
    stringsAsFactors = FALSE)

  ## --- the four label rules, on the 1%-filtered export --------------
  f <- pg_ent[pg_ent$q <= ALPHA_MAIN, , drop = FALSE]
  agg <- f %>%
    group_by(id) %>%
    summarise(q = min(q), n = n(), nD = sum(is_decoy),
              win_is_D = is_decoy[which.min(q)], .groups = "drop")

  rule_label <- list(
    any_decoy   = agg$nD > 0,
    winning_row = agg$win_is_D,
    majority    = agg$nD / agg$n > 0.5)

  rows <- lapply(names(rule_label), function(nm) {
    isD <- rule_label[[nm]]
    tgt <- agg$id[!isD]
    E <- sum(is_pure(tgt)); Tn <- length(tgt); Dn <- sum(isD)
    w <- wilson(E, Tn)
    data.frame(rule = nm, T_alpha = Tn, D_alpha = Dn,
               ratio = Dn / Tn, E_alpha = E,
               entrap = unname(w["est"]), lo = unname(w["lo"]),
               hi = unname(w["hi"]), stringsAsFactors = FALSE)
  })

  ## split: target and decoy entries are two hypotheses
  sp <- pg_ent %>% group_by(id, is_decoy) %>%
    summarise(q = min(q), .groups = "drop")
  acc_sp <- sp[sp$q <= ALPHA_MAIN, ]
  tgt <- acc_sp$id[!acc_sp$is_decoy]
  E <- sum(is_pure(tgt)); Tn <- length(tgt); Dn <- sum(acc_sp$is_decoy)
  w <- wilson(E, Tn)
  rows[[length(rows) + 1]] <- data.frame(
    rule = "split", T_alpha = Tn, D_alpha = Dn, ratio = Dn / Tn,
    E_alpha = E, entrap = unname(w["est"]), lo = unname(w["lo"]),
    hi = unname(w["hi"]), stringsAsFactors = FALSE)

  ## picked
  pke <- picked_protein(pg_ent)
  E <- sum(is_pure(pke$accepted_targets))
  w <- wilson(E, pke$T_alpha)
  rows[[length(rows) + 1]] <- data.frame(
    rule = "picked_protein", T_alpha = pke$T_alpha, D_alpha = pke$D_alpha,
    ratio = pke$ratio, E_alpha = E, entrap = unname(w["est"]),
    lo = unname(w["lo"]), hi = unname(w["hi"]), stringsAsFactors = FALSE)

  out <- do.call(rbind, rows)
  print(out)
  write_tab(out, "label_rules_vs_entrapment")
  write_tab(data.frame(rule = "picked_protein", T_alpha = pke$T_alpha,
                       D_alpha = pke$D_alpha, ratio = pke$ratio,
                       entrap = out$entrap[out$rule == "picked_protein"]),
            "picked_protein_entrapment")

  reg("pickedEntrapT",   pke$T_alpha, "int", "Table S: picked, entrapment search")
  reg("pickedEntrapFDP", out$entrap[out$rule == "picked_protein"], "d3",
      "Table S: picked-protein entrapment FDP")
  reg("splitEntrapFDP",  out$entrap[out$rule == "split"], "d3",
      "Table S: split-rule entrapment FDP")

  ## The headline reading, printed so it cannot be missed.
  lab <- out[out$rule != "picked_protein", ]
  cat("\n--- reading for Section 4.2 ---\n")
  cat(sprintf(
"The four label rules span %.3f%%-%.3f%% in the reported decoy ratio but only\n%.3f%%-%.3f%% in the entrapment FDP of the accepted target list: the rule changes\nwhat is REPORTED. Picked competition accepts a different list (T=%d against\n%d-%d) at an entrapment FDP of %.3f%%: the competition scheme changes what is\nACCEPTED.\n\n",
    100 * min(lab$ratio), 100 * max(lab$ratio),
    100 * min(lab$entrap), 100 * max(lab$entrap),
    pke$T_alpha, min(lab$T_alpha), max(lab$T_alpha),
    100 * out$entrap[out$rule == "picked_protein"]))
}

save_numbers()
message("\nNext: 04_calibration_stability.R")
