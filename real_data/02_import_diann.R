## =====================================================================
##  02 --- Import the DIA-NN searches and build every universe
##
##  Universes built, per search:
##    runxprec     diann_runxprecursor(q_col = "Q.Value", id_mode="runxid")
##    global       diann_global_precursor(q_col = "Global.Q.Value")
##    minrunq      diann_global_minrunq(q_col = "Q.Value")     <- scope misuse
##    pg_global    protein groups on Global.PG.Q.Value
##    pg_run       protein groups on PG.Q.Value
##  plus the precursor and protein universes under all three
##  label-propagation rules.
##
##  Output: out/universes.rds
##
##  TWO THINGS TO KNOW, both verified on report_basic.parquet:
##
##  1. diann_global_precursor() DEFAULTS TO q_col = "Q.Value". The global
##     universe therefore requires q_col = "Global.Q.Value" explicitly.
##     Getting this wrong silently produces a run-wise list labelled
##     "global" -- the exact error the paper is about.
##
##  2. At PROTEIN level the decoy-if-any rule used by the package's
##     collapse functions is wrong: 10,262 of 13,751 groups contain both
##     target and decoy precursors, so almost every group is labelled a
##     decoy and the ratio reads ~840% instead of ~1%. Protein universes
##     are therefore built with PROTEIN_LABEL_RULE = "split".
## =====================================================================

source("00_config.R")

need_file(DIANN_BASIC, "DIA-NN basic report")
need_file(DIANN_SEMI,  "DIA-NN semi-tryptic report")

LABEL_RULES <- c("any_decoy", "winning_row", "majority")

## q_max_export is a property of the export, so read it from the data
## rather than trusting 00_config.R. Getting it wrong silently drops
## local-tail windows and equal-chance bands.
qmax_of <- function(rep) {
  lapply(c(q_run = COL$q_run, q_global = COL$q_global,
           q_pg = COL$q_pg, q_pg_glob = COL$q_pg_glob),
         function(cn) if (cn %in% names(rep))
           suppressWarnings(max(as.numeric(rep[[cn]]), na.rm = TRUE)) else NA_real_)
}

import_one <- function(path, tag) {
  message("\n== ", tag, " : ", basename(path))
  rep <- diagFDR::read_diann_parquet(path)
  say("%s rows, %d columns", format(nrow(rep), big.mark = ","), ncol(rep))

  miss <- setdiff(unlist(COL[c("run","precursor","protein","q_run","q_global","decoy")]),
                  names(rep))
  if (length(miss))
    stop("Missing columns in ", tag, ": ", paste(miss, collapse = ", "),
         "\nRun 01_preflight.R to see what the export actually contains.")

  ## Drop rows the package cannot handle BEFORE it sees them. diann_*()
  ## computes is_decoy = any(is_decoy) without filtering NA, so a single
  ## NA in Decoy becomes an NA label and aborts inside validate_dfdr_tbl()
  ## several steps later, far from the cause.
  rep <- sanitise_report(rep, COL$decoy, tag)

  ## q_max_export read from the data, not from config: it is a property of
  ## the export and getting it wrong silently drops local-tail windows and
  ## equal-chance bands.
  qmax <- qmax_of(rep)
  say("q ceilings: %s", paste(sprintf("%s=%.4g", names(qmax), unlist(qmax)), collapse = "  "))
  for (k in names(qmax)) if (is.finite(qmax[[k]]) && is.finite(QMAX[[k]]) &&
                             abs(qmax[[k]] - QMAX[[k]]) > 1e-3)
    warning(sprintf("%s: QMAX$%s is %.4g in 00_config.R but %.4g in the data; using the data.",
                    tag, k, QMAX[[k]], qmax[[k]]), call. = FALSE)

  ## --- label-intrinsic check, at BOTH levels --------------------------
  li_prec <- label_is_intrinsic(rep, COL$precursor)
  li_prot <- label_is_intrinsic(rep, COL$protein)
  say("precursors with both labels: %s of %s%s",
      format(li_prec$n_mixed, big.mark = ","), format(li_prec$n_keys, big.mark = ","),
      if (li_prec$intrinsic || li_prec$n_mixed / li_prec$n_keys < 1e-4)
        "  -> effectively intrinsic; the three rules coincide" else "")
  say("protein groups with both labels: %s of %s%s",
      format(li_prot$n_mixed, big.mark = ","), format(li_prot$n_keys, big.mark = ","),
      if (!li_prot$intrinsic) "  -> decoy-if-any is INVALID here; using 'split'" else "")

  ## --- precursor universes, via the package ---------------------------
  u <- list(
    runxprec = diagFDR::diann_runxprecursor(
      rep, q_col = COL$q_run, run_col = COL$run, pep_col = COL$pep,
      id_mode = "runxid", q_max_export = qmax$q_run,
      unit = "runxprecursor", scope = "runwise"),
    global = diagFDR::diann_global_precursor(
      rep, q_col = COL$q_global, pep_col = COL$pep,
      q_max_export = qmax$q_global, unit = "precursor", scope = "global"),
    minrunq = diagFDR::diann_global_minrunq(
      rep, run_col = COL$run, q_col = COL$q_run, pep_col = COL$pep,
      q_max_export = qmax$q_run, unit = "precursor", scope = "aggregated"))

  ## --- protein-group universes, split rule ----------------------------
  u$pg_global <- collapse_universe(rep, COL$protein, COL$q_pg_glob,
                                   rule = PROTEIN_LABEL_RULE, pep_col = NULL,
                                   q_max_export = qmax$q_pg_glob,
                                   unit = "protein group", scope = "global")
  u$pg_run <- collapse_universe(rep, COL$protein, COL$q_pg,
                                rule = PROTEIN_LABEL_RULE, pep_col = NULL,
                                q_max_export = qmax$q_pg,
                                unit = "protein group", scope = "runwise")

  ## --- all label rules, at precursor AND protein level --------------
  u$minrunq_by_rule <- setNames(lapply(LABEL_RULES, function(r)
    collapse_universe(rep, COL$precursor, COL$q_run, rule = r,
                      pep_col = COL$pep, q_max_export = qmax$q_run,
                      unit = "precursor", scope = "aggregated")), LABEL_RULES)
  u$pg_by_rule <- setNames(lapply(c(LABEL_RULES, "split"), function(r)
    collapse_universe(rep, COL$protein, COL$q_pg_glob, rule = r, pep_col = NULL,
                      q_max_export = qmax$q_pg_glob,
                      unit = "protein group", scope = "global")),
    c(LABEL_RULES, "split"))

  for (n in c("runxprec", "global", "minrunq", "pg_global", "pg_run")) {
    h <- headline_of(u[[n]])
    say("%-10s T=%9s D=%7s ratio=%7.3f%% CV=%.4f", n,
        format(h$T_alpha, big.mark = ","), format(h$D_alpha, big.mark = ","),
        100 * h$decoy_ratio, h$CV)
  }

  attr(u, "label_prec") <- li_prec
  attr(u, "label_prot") <- li_prot
  attr(u, "n_rows") <- nrow(rep)
  attr(u, "n_runs") <- length(unique(rep[[COL$run]]))
  u
}

U <- list(basic = import_one(DIANN_BASIC, "basic"),
          semi  = import_one(DIANN_SEMI,  "semi-tryptic"))

## --- first-pass reports, for the Figure 3 first-pass vs reanalysis panel
for (nm in c("basic", "semi")) {
  fp <- if (nm == "basic") DIANN_BASIC_FP else DIANN_SEMI_FP
  if (!is.na(fp) && file.exists(fp)) {
    message("\n== ", nm, " first pass")
    rep <- sanitise_report(diagFDR::read_diann_parquet(fp), COL$decoy, paste(nm, "first pass"))
    qmax <- qmax_of(rep)
    U[[paste0(nm, "_fp")]] <- list(
      global = diagFDR::diann_global_precursor(
        rep, q_col = COL$q_global, pep_col = COL$pep,
        q_max_export = qmax$q_global, unit = "precursor", scope = "global"),
      runxprec = diagFDR::diann_runxprecursor(
        rep, q_col = COL$q_run, run_col = COL$run, pep_col = COL$pep,
        id_mode = "runxid", q_max_export = qmax$q_run,
        unit = "runxprecursor", scope = "runwise"))
    h <- headline_of(U[[paste0(nm, "_fp")]]$global)
    say("first-pass global T=%s D=%s ratio=%.3f%%",
        format(h$T_alpha, big.mark = ","), format(h$D_alpha, big.mark = ","),
        100 * h$decoy_ratio)
  } else {
    message("\n[note] first-pass report not found for ", nm,
            " -- Figure 3's first-pass/reanalysis comparison will be omitted.")
  }
}

## --- entrapment search -------------------------------------------------
if (file.exists(DIANN_ENTRAP)) {
  message("\n== entrapment search")
  de <- sanitise_report(diagFDR::read_diann_parquet(DIANN_ENTRAP), COL$decoy, "entrapment")
  qmax <- qmax_of(de)
  say("%s rows", format(nrow(de), big.mark = ","))
  U$entrap <- list(
    runxprec = diagFDR::diann_runxprecursor(
      de, q_col = COL$q_run, run_col = COL$run, pep_col = COL$pep,
      id_mode = "runxid", q_max_export = qmax$q_run),
    global = diagFDR::diann_global_precursor(
      de, q_col = COL$q_global, pep_col = COL$pep, q_max_export = qmax$q_global),
    minrunq = diagFDR::diann_global_minrunq(
      de, run_col = COL$run, q_col = COL$q_run, pep_col = COL$pep,
      q_max_export = qmax$q_run),
    pg_global = collapse_universe(de, COL$protein, COL$q_pg_glob,
                                  rule = PROTEIN_LABEL_RULE, pep_col = NULL,
                                  q_max_export = qmax$q_pg_glob,
                                  unit = "protein group", scope = "global"))
  ## precursor -> protein map, for entrapment attribution in 05
  keep <- !duplicated(de[[COL$precursor]])
  U$entrap$prot_map <- data.frame(
    id      = as.character(de[[COL$precursor]][keep]),
    protein = as.character(de[[COL$protein]][keep]),
    protein_ids = as.character(de[[COL$protein_ids]][keep]),
    stringsAsFactors = FALSE)
  say("entrapment universes built")
} else {
  message("\n[note] no entrapment report at ", DIANN_ENTRAP)
}

saveRDS(U, file.path(OUT, "universes.rds"))
message("\nSaved -> ", file.path(OUT, "universes.rds"))
message("Next: 03_scope_audit.R")
