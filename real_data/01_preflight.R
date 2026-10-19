## =====================================================================
##  01 --- PREFLIGHT.  Run this first. It changes nothing.
##
##  Checks the installed diagFDR against the exports these scripts use, every
##  data file, every column name the scripts use, the q-value export
##  ceilings, and the entrapment accession attribution.
## =====================================================================

source("00_config.R")

log <- character(0)
note <- function(...) { s <- sprintf(...); log <<- c(log, s); cat(s, "\n") }
ok   <- function(...) note("  [ok]   %s", sprintf(...))
bad  <- function(...) note("  [FAIL] %s", sprintf(...))
warn_<- function(...) note("  [warn] %s", sprintf(...))

note("=== diagFDR paper preflight ===")
note("R %s | diagFDR %s | %s", getRversion(), DFDR_VERSION, format(Sys.time()))
note("")

## --- 1. package --------------------------------------------------------
note("--- 1. diagFDR ---")
exp_ <- getNamespaceExports("diagFDR")
miss <- setdiff(DFDR_EXPORTS_EXPECTED, exp_)
if (!length(miss)) {
  ok("all %d expected exports present", length(DFDR_EXPORTS_EXPECTED))
} else {
  bad("%d expected export(s) missing:", length(miss))
  for (m in miss) note("           %s", m)
  warn_("Edit the matching wrapper in R/dfdr_api.R.")
}
extra <- setdiff(exp_, DFDR_EXPORTS_EXPECTED)
if (length(extra)) note("  [info] %d export(s) not used here: %s", length(extra),
                        paste(head(extra, 20), collapse = ", "))
reg("diagFDRversion", DFDR_VERSION, note = "package version used for the analysis")

## Signature check on the functions the scripts depend on most.
for (fn in intersect(c("dfdr_headline","dfdr_curve_stability","dfdr_local_tail",
                       "dfdr_elasticity","dfdr_pep_reliability_tdc",
                       "dfdr_equal_chance_qbands","dfdr_sumpep",
                       "diann_global_precursor","diann_runxprecursor",
                       "diann_global_minrunq","mokapot_competed_universe"), exp_)) {
  a <- names(formals(getExportedValue("diagFDR", fn)))
  note("  [info] %s(%s)", fn, paste(a, collapse = ", "))
}
note("")

## --- 2. other packages -------------------------------------------------
note("--- 2. dependencies ---")
for (p in c("arrow","dplyr","tibble","purrr","readr","ggplot2","rlang","cp4p","irr")) {
  if (requireNamespace(p, quietly = TRUE)) {
    ok("%-8s %s", p, utils::packageVersion(p))
  } else if (p %in% c("irr")) {
    warn_("%s not installed (only for a kappa cross-check)", p)
  } else {
    bad("%s NOT installed -- required", p)
  }
}
note("")

## --- 3. data files -----------------------------------------------------
note("--- 3. data files ---")
files <- list(
  "DIA-NN basic" = DIANN_BASIC, "DIA-NN semi" = DIANN_SEMI,
  "DIA-NN entrapment" = DIANN_ENTRAP,
  "DIA-NN basic 1st pass" = DIANN_BASIC_FP, "DIA-NN semi 1st pass" = DIANN_SEMI_FP,
  "FASTA target" = FASTA_TARGET, "FASTA contaminants" = FASTA_CONTAM,
  "FASTA entrapment" = FASTA_ENTRAP,
  "mokapot targets" = MOKAPOT_TARGETS, "mokapot decoys" = MOKAPOT_DECOYS,
  "MaxQuant msms.txt" = MAXQUANT_MSMS)
optional <- c("DIA-NN basic 1st pass","DIA-NN semi 1st pass")
for (nm in names(files)) {
  f <- files[[nm]]
  if (!is.na(f) && file.exists(f)) {
    ok("%-22s %-42s %8.1f MB", nm, basename(f), file.info(f)$size / 1e6)
  } else if (nm %in% optional) {
    warn_("%-22s not found (optional): %s", nm, f)
  } else {
    bad("%-22s NOT FOUND: %s", nm, f)
  }
}
note("")

## --- 4. DIA-NN columns and export ceilings -----------------------------
note("--- 4. DIA-NN export ---")
if (file.exists(DIANN_BASIC)) {
  rep <- tryCatch(diagFDR::read_diann_parquet(DIANN_BASIC),
                  error = function(e) { bad("cannot read: %s", conditionMessage(e)); NULL })
  if (!is.null(rep)) {
    note("  %s rows, %d columns, %d runs", format(nrow(rep), big.mark = ","),
         ncol(rep), length(unique(rep[[COL$run]])))
    reg("diannRows", nrow(rep), "int", "rows in the basic export, for the Implementation sentence")
    reg("diannRuns", length(unique(rep[[COL$run]])), "int", "number of runs")
    for (k in names(COL)) {
      v <- COL[[k]]
      if (v %in% names(rep)) ok("%-12s -> %s", k, v) else bad("%-12s -> '%s' NOT PRESENT", k, v)
    }
    dcv <- rep[[COL$decoy]]
    note("  Decoy column: class %s, %s NA, values {%s}", paste(class(dcv), collapse="/"),
         format(sum(is.na(dcv)), big.mark = ","),
         paste(head(sort(unique(dcv)), 4), collapse = ","))
    if (anyNA(dcv))
      bad("Decoy contains %s NA. 02 drops those rows; check that is intended.",
          format(sum(is.na(dcv)), big.mark = ","))
    if (!"CScore" %in% names(rep))
      warn_("no CScore column: `score` stays NA")
    for (k in names(QMAX)) {
      cn <- COL[[k]]
      if (!is.null(cn) && cn %in% names(rep)) {
        m <- suppressWarnings(max(as.numeric(rep[[cn]]), na.rm = TRUE))
        ## Tolerance is relative: parquet stores these as float32, so the
        ## max of a column capped at 0.5 need not equal the double 0.5
        ## exactly. 02 reads the ceiling from the data anyway; this is only
        ## a sanity cross-check against 00_config.R.
        if (abs(m - QMAX[[k]]) <= 1e-4 * max(1, QMAX[[k]])) {
          ok("q_max_export %-18s = %.10g", cn, m)
        } else {
          warn_("q_max_export %s: config says %.10g, data says %.10g. 02 uses the data value; update QMAX only if this surprises you.",
                cn, QMAX[[k]], m)
        }
      }
    }
    li_p <- label_is_intrinsic(rep, COL$precursor)
    li_g <- label_is_intrinsic(rep, COL$protein)
    note("  precursors with both labels: %s / %s", format(li_p$n_mixed, big.mark = ","),
         format(li_p$n_keys, big.mark = ","))
    note("  protein groups with both labels: %s / %s", format(li_g$n_mixed, big.mark = ","),
         format(li_g$n_keys, big.mark = ","))
    if (!li_g$intrinsic && PROTEIN_LABEL_RULE != "split") {
      bad("protein groups are label-mixed but PROTEIN_LABEL_RULE is '%s'. Set it to 'split'.",
          PROTEIN_LABEL_RULE)
    } else if (!li_g$intrinsic) {
      ok("PROTEIN_LABEL_RULE = 'split' (correct for this export)")
    }
    rm(rep); invisible(gc())
  }
}
note("")

## --- 5. entrapment attribution ----------------------------------------
note("--- 5. entrapment attribution ---")
if (file.exists(FASTA_ENTRAP)) {
  h <- readLines(FASTA_ENTRAP, n = 2, warn = FALSE)
  note("  first header: %s", substr(h[1], 1, 90))
  accs <- entrap_accessions(FASTA_ENTRAP)
  ok("%s entrapment accessions parsed", format(length(accs), big.mark = ","))
  reg("entrapAccessions", length(accs), "int", "entrapment database size")
  if (file.exists(DIANN_ENTRAP)) {
    de <- diagFDR::read_diann_parquet(DIANN_ENTRAP)
    pgv <- as.character(de[[COL$protein]])
    n_hit <- sum(vapply(strsplit(pgv, ";"), function(z) any(z %in% accs), logical(1)))
    if (n_hit > 0) {
      ok("%s rows carry >=1 entrapment accession in %s", format(n_hit, big.mark = ","), COL$protein)
    } else {
      bad("NO row carries an entrapment accession in %s -- check ENTRAP_COL", COL$protein)
    }
    if (any(grepl("ARATH", pgv, fixed = TRUE))) {
      note("  [info] species tags present in Protein.Group")
    } else warn_("Protein.Group holds bare accessions (no species tag): a string prefix such as 'ARATH' would match NOTHING. Attribution is by accession membership, which is what 05_entrapment.R does.")
    rm(de); invisible(gc())
  }
} else bad("entrapment FASTA not found")
note("")

## --- 6. mokapot columns ------------------------------------------------
note("--- 6. mokapot export ---")
if (file.exists(MOKAPOT_TARGETS)) {
  nms <- names(read.delim(MOKAPOT_TARGETS, nrows = 1, check.names = FALSE))
  req <- c("run","spectrum_id","mokapot score","mokapot q-value","mokapot PEP","is_target")
  for (r in req) { if (r %in% nms) ok("%s", r) else bad("missing mokapot column: %s", r) }
  note("  columns: %s", paste(nms, collapse = ", "))
  warn_("no pre-rescoring score in the mokapot export: the 'before' universe comes from msms.txt via read_dfdr_maxquant_msms()")
}
note("")

## --- 7. environment, for Note S8 --------------------------------------
note("--- 7. environment (Supplementary Note S8) ---")
si <- sessionInfo()
note("  %s", si$R.version$version.string)
note("  platform: %s", si$platform)
note("  running:  %s", si$running)
note("  cores: %s", tryCatch(parallel::detectCores(), error = function(e) NA))
note("")

save_numbers()
writeLines(log, file.path(OUT, "preflight.txt"))
nfail <- sum(grepl("\\[FAIL\\]", log))
if (nfail == 0) {
  message("\nPREFLIGHT PASSED. Next: 02_import_diann.R")
} else {
  message("\nPREFLIGHT: ", nfail, " problem(s). See ", file.path(OUT, "preflight.txt"))
}
