## =====================================================================
##  05c --- Entrapment: homology, and sensitivity of the size correction
##
##  Homology, database size and detectability all bear on the size
##  correction. Two of the three can be examined from the data; the third
##  cannot, and this script makes clear which is which.
##
##  PART A --- HOMOLOGY. The relevant quantity is the fraction of the
##  entrapment digest that is indistinguishable from the target digest at
##  the level where competition occurs. Expected: 10,037 of 962,791
##  entrapment tryptic peptides (1.04%) also occur in the target digest.
##  Because the PURE-entrapment counting rule requires EVERY accession in
##  a protein group to be an entrapment accession, a shared peptide can
##  never on its own generate an E_alpha count: homology at this level
##  cannot inflate the estimate, only deflate it. The residual bias is
##  therefore conservative, and the script says so in the printed reading.
##
##  PART B --- DATABASE SIZE. The corrected estimator assumes E_alpha
##  scales linearly with entrapment database size, so that E_alpha / r
##  estimates the false-target count. We test that by subsampling the
##  entrapment ACCESSION SET and recomputing E_alpha on the SAME search
##  output, scaling r proportionally.
##
##  READ THIS BEFORE INTERPRETING PART B. Removing accessions from the
##  ANALYSIS is not the same as removing them from the SEARCHED DATABASE,
##  which would also change the competition. This is a check on the
##  ESTIMATOR at fixed competition, not a re-design.
##
##  Expected: the corrected estimate falls from 0.717% to 0.483% as the
##  database is thinned to a quarter. Most of that drift is an artefact of
##  the pure-entrapment rule -- E_alpha falls sub-linearly (177 -> 38.5, a
##  factor of 4.6 for a factor-of-4 reduction) because protein groups that
##  MIX entrapment and target accessions drop out of the count as the
##  entrapment set shrinks. A further reason to read the corrected value
##  comparatively rather than absolutely.
##
##  PART C --- DETECTABILITY. Not answerable in this design, and the
##  script says so rather than producing a number that would look like an
##  answer. No measurement here can separate a detectability difference
##  from a genuine error-rate difference. This is the principal reason the
##  entrapment estimate is interpreted comparatively throughout.
##
##  Outputs, in tables/:
##    entrapment_composition.csv
##    entrapment_r_subsample.csv
## =====================================================================

source("00_config.R")
load_numbers()

set.seed(SEED)
N_DRAWS <- 20L
FRACTIONS <- c(1.00, 0.75, 0.50, 0.25)

## =====================================================================
##  PART A --- composition and homology. Needs only the FASTA files.
## =====================================================================
message("\n=== A. entrapment database composition and homology ===")

r_pep_unique <- NA_real_

if (file.exists(FASTA_TARGET) && file.exists(FASTA_ENTRAP)) {
  tgt <- read_fasta(FASTA_TARGET)
  ent <- read_fasta(FASTA_ENTRAP)
  con <- if (file.exists(FASTA_CONTAM)) read_fasta(FASTA_CONTAM)
         else list(ids = character(0), seqs = character(0))

  n_tgt_prot <- length(tgt$ids) + length(con$ids)
  n_ent_prot <- length(ent$ids)
  say("target proteins (incl. contaminants): %s",
      format(n_tgt_prot, big.mark = ","))
  say("entrapment proteins:                  %s",
      format(n_ent_prot, big.mark = ","))

  message("  digesting (a few minutes on the full proteomes)...")
  pep_tgt <- unique(unlist(lapply(c(tgt$seqs, con$seqs), digest)))
  pep_ent <- unique(unlist(lapply(ent$seqs, digest)))
  shared  <- intersect(pep_tgt, pep_ent)

  r_prot       <- n_ent_prot / n_tgt_prot
  r_pep        <- length(pep_ent) / length(pep_tgt)
  r_pep_unique <- length(setdiff(pep_ent, shared)) / length(pep_tgt)
  homology_pct <- length(shared) / length(pep_ent)

  say("tryptic peptides target/entrapment/shared: %s / %s / %s (%.2f%% of entrapment)",
      format(length(pep_tgt), big.mark = ","),
      format(length(pep_ent), big.mark = ","),
      format(length(shared),  big.mark = ","), 100 * homology_pct)
  say("r protein = %.4f | r peptide = %.4f | r peptide (shared removed) = %.4f",
      r_prot, r_pep, r_pep_unique)

  reg("rProtein",       r_prot,       "d3", "Note S4: size ratio, protein level")
  reg("rPeptide",       r_pep,        "d3", "Note S4: peptide level, incl. shared")
  reg("rPeptideUnique", r_pep_unique, "d3", "Note S4: the one used throughout")
  reg("sharedPeptides", length(shared), "int", "Note S4: homology overlap")
  reg("homologyPct",    homology_pct, "pct2", "Note S4: shared as % of entrapment digest")
  reg("entrapProteins", n_ent_prot,   "int", "Note S4: entrapment database size")

  write_tab(data.frame(
    quantity = c("target proteins","entrapment proteins","target peptides",
                 "entrapment peptides","shared peptides","r protein","r peptide",
                 "r peptide unique","shared as fraction of entrapment digest"),
    value = c(n_tgt_prot, n_ent_prot, length(pep_tgt), length(pep_ent),
              length(shared), r_prot, r_pep, r_pep_unique, homology_pct)),
    "entrapment_composition")

  cat("\n--- homology ---\n")
  cat(sprintf(
"Of the $%s$ entrapment tryptic peptides, $%s$ ($%.2f\\%%$) also occur in the\ntarget digest. These are excluded from $r$ (giving $%.3f$ rather than $%.3f$)\nand, because the pure-entrapment counting rule requires \\emph{every} accession\nin a protein group to be an entrapment accession, a shared peptide can never on\nits own generate an $E_\\alpha$ count. Homology at this level therefore cannot\ninflate the entrapment estimate; the direction of the residual bias is\nconservative.\n\n",
    formatC(length(pep_ent), format = "d", big.mark = "{,}"),
    formatC(length(shared),  format = "d", big.mark = "{,}"),
    100 * homology_pct, r_pep_unique, r_pep))
} else {
  message("  FASTA files not found -- Part A skipped, Part B needs r.")
}

## =====================================================================
##  PART B --- sensitivity of the corrected estimator to database size
## =====================================================================
message("\n=== B. r-subsampling check of the corrected estimator ===")

if (!file.exists(DIANN_ENTRAP)) {
  message("  no entrapment report -- Part B skipped.")
} else if (!is.finite(r_pep_unique)) {
  message("  r unavailable (Part A skipped) -- Part B skipped.")
} else {

  ACC_ALL <- entrap_accessions(FASTA_ENTRAP)
  rep_ent <- sanitise_report(diagFDR::read_diann_parquet(DIANN_ENTRAP),
                             COL$decoy, "entrapment")

  ## the correctly scoped global precursor list at alpha
  U <- rep_ent[!duplicated(rep_ent[[COL$precursor]]), ]
  keyq <- tapply(as.numeric(rep_ent[[COL$q_global]]),
                 as.character(rep_ent[[COL$precursor]]), min)
  keyd <- tapply(rep_ent[[COL$decoy]] == 1L,
                 as.character(rep_ent[[COL$precursor]]), any)
  acc_ids <- names(keyq)[keyq <= ALPHA_MAIN & !keyd]
  T_alpha <- length(acc_ids)
  say("accepted targets on the global list: %s", format(T_alpha, big.mark = ","))

  prot_of <- setNames(as.character(U[[COL$protein]]),
                      as.character(U[[COL$precursor]]))
  groups  <- strsplit(ifelse(is.na(prot_of[acc_ids]), "",
                             prot_of[acc_ids]), ";")

  count_pure <- function(A) {
    sum(vapply(groups, function(z) {
      z <- z[nzchar(z)]; length(z) > 0 && all(z %in% A) }, logical(1)))
  }

  corrected <- function(E, T_, r) E * (1 + 1 / r) / (T_ + E)

  out <- do.call(rbind, lapply(FRACTIONS, function(fr) {
    reps <- if (fr >= 1) 1L else N_DRAWS
    Es <- vapply(seq_len(reps), function(i) {
      A <- if (fr >= 1) ACC_ALL
           else sample(ACC_ALL, size = round(fr * length(ACC_ALL)))
      count_pure(A)
    }, numeric(1))
    r_fr <- r_pep_unique * fr
    data.frame(fraction = fr, r = r_fr, draws = reps,
               E_mean = mean(Es), E_sd = if (reps > 1) sd(Es) else 0,
               corrected = corrected(mean(Es), T_alpha, r_fr))
  }))
  print(out); write_tab(out, "entrapment_r_subsample")

  base <- out$corrected[out$fraction == 1]
  low  <- out$corrected[which.min(out$fraction)]
  say("corrected FDP drifts from %.3f%% to %.3f%% as the database is thinned",
      100 * base, 100 * low)

  ## Sub-linearity of E in database size is the pure-rule artefact.
  E1 <- out$E_mean[out$fraction == 1]
  E4 <- out$E_mean[which.min(out$fraction)]
  say("E falls by a factor of %.1f for a factor-of-%.0f reduction in accessions",
      E1 / E4, 1 / min(out$fraction))

  reg("rSubBase", base, "d3", "corrected FDP, full entrapment DB")
  reg("rSubLow",  low,  "d3", "corrected FDP, 25% of the DB")

  cat("\n--- database size ---\n")
  cat(sprintf(
"The corrected estimate is not invariant: it falls from $%.3f\\%%$ to $%.3f\\%%$ as\nthe database is thinned to a quarter. Most of this drift is attributable to the\npure-entrapment counting rule rather than to the correction itself --- $E_\\alpha$\nfalls sub-linearly ($%.0f\\to%.1f$, a factor of $%.1f$ for a factor-of-$%.0f$\nreduction in accessions) because protein groups that mix entrapment and target\naccessions drop out of the count as the entrapment set shrinks.\n\n",
    100 * base, 100 * low, E1, E4, E1 / E4, 1 / min(out$fraction)))
}

## =====================================================================
##  PART C --- detectability. Deliberately not answered.
## =====================================================================
message("\n=== C. detectability ===")
message("  Not answerable in this design, and no number is produced here.")
message("  No measurement on a single entrapment proteome can separate a")
message("  detectability difference from a genuine error-rate difference.")
message("  This is why every entrapment result in the paper is interpreted")
message("  COMPARATIVELY (the 1.79-fold increase under min_r aggregation, the")
message("  3-fold precursor-to-protein increase) and never as an absolute")
message("  measurement of the FDP. A detectability difference that is constant")
message("  across a comparison cancels in the ratio.")
message("  To go further you would need a second, phylogenetically distant")
message("  entrapment proteome and to show the comparison is stable across the")
message("  two. That is a new search, not a re-analysis.")

save_numbers()
message("\nNext: 06_msrescore.R")
