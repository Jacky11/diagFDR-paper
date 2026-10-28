## =====================================================================
##  05e --- Entrapment (Section 4.2, Note S4) recomputed with diagFDR >= 0.2.0
##
##  Reproduces the raw and corrected estimates, the Wilson intervals, the
##  size ratio r and the scope-discordant-set rate from the package
##  functions, and ASSERTS the paper's numbers.
##
##  Expected at alpha = 1% (entrapment search):
##    global   T = 81,643  E = 177  raw 0.217% [0.187, 0.251]  corr 0.719%
##    min_r    T = 88,064  E = 342  raw 0.388% [0.349, 0.432]  corr 1.287%
##    added set 6,421 precursors, E = 165, 2.57%; ratio to core 11.9
##    protein groups (split): T = 8,751  E = 58  raw 0.663%
##    r: protein 0.155, peptide 0.436, peptide without shared 0.432
##  Verified 2026-10-27 with diagFDR 0.2.1 (corrected estimate divides by T_alpha,
##  which already counts the entrapment matches; 0.2.0 divided by T + E) (counts, CIs, corrected values,
##  discordant set and alpha grid); see manuscript/verification_report.md.
## =====================================================================
source("00_config.R")
load_numbers()

ACC <- diagFDR::read_fasta_accessions(FASTA_ENTRAP)
say("%s entrapment accessions", format(length(ACC), big.mark = ","))

rep <- sanitise_report(diagFDR::read_diann_parquet(DIANN_ENTRAP), COL$decoy, "entrapment")
ent_prec <- diagFDR::diann_entrapment_ids(rep, ACC, id_col = COL$precursor,
                                          protein_col = COL$protein, rule = ENTRAP_RULE)
ent_pg   <- diagFDR::diann_entrapment_ids(rep, ACC, id_col = COL$protein,
                                          protein_col = COL$protein, rule = ENTRAP_RULE)
say("%s entrapment precursors, %s entrapment protein groups", length(ent_prec), length(ent_pg))

## --- size ratio r (slow: digests both proteomes; cache the result) ----
r_file <- file.path(TABDIR, "entrapment_size_ratio_package.csv")
if (!file.exists(r_file)) {
  rr <- diagFDR::dfdr_entrapment_size_ratio(c(FASTA_TARGET, FASTA_CONTAM), FASTA_ENTRAP,
                                            sites = c("K", "R"), no_cut_before = character(0),
                                            missed = 1, min_len = 6, max_len = 30)   # same rule as utils.R digest()
  write.csv(rr, r_file, row.names = FALSE)
} else rr <- read.csv(r_file)
print(rr)
r_pep <- rr$value[rr$quantity == "r_peptide_unique"]
r_pro <- rr$value[rr$quantity == "r_protein"]
stopifnot(abs(r_pep - 0.432) < 0.002, abs(r_pro - 0.155) < 0.002)

## --- universes -----------------------------------------------------
U <- list(
  global  = diagFDR::diann_global_precursor(rep, q_col = COL$q_global, q_max_export = QMAX$q_global, two_pass = TRUE),
  minrunq = diagFDR::diann_global_minrunq(rep, q_col = COL$q_run, q_max_export = QMAX$q_run, two_pass = TRUE),
  pg      = diagFDR::diann_global_proteingroup(rep, q_col = COL$q_pg_glob, rule = PROTEIN_LABEL_RULE,
                                               q_max_export = QMAX$q_pg_glob, two_pass = TRUE))

e_glob <- diagFDR::dfdr_entrapment(U$global,  ent_prec, alpha = ALPHA_MAIN, r = r_pep)
e_minr <- diagFDR::dfdr_entrapment(U$minrunq, ent_prec, alpha = ALPHA_MAIN, r = r_pep)
e_pg   <- diagFDR::dfdr_entrapment(U$pg, ent_pg, alpha = ALPHA_MAIN, r = r_pep)
e_pg_rprot <- diagFDR::dfdr_entrapment(U$pg, ent_pg, alpha = ALPHA_MAIN, r = r_pro)
out <- rbind(cbind(universe = "global", e_glob), cbind(universe = "minrunq", e_minr),
             cbind(universe = "protein_group(r_pep)", e_pg), cbind(universe = "protein_group(r_prot)", e_pg_rprot))
print(out); write_tab(out, "entrapment_package")

stopifnot(e_glob$T_alpha == 81643, e_glob$E_alpha == 177,
          e_minr$T_alpha == 88064, e_minr$E_alpha == 342,
          e_pg$T_alpha == 8751, e_pg$E_alpha == 58,
          abs(100 * e_glob$fdp_corrected - 0.719) < 0.001,   # paper: 0.719%, 1.29%, 2.20%, 4.94%
          abs(100 * e_minr$fdp_corrected - 1.288) < 0.002,
          abs(100 * e_pg$fdp_corrected - 2.197) < 0.002,
          abs(100 * e_pg_rprot$fdp_corrected - 4.94) < 0.01)   # r_protein = 0.1548 exactly
say("Section 4.2 / 4.3 entrapment counts reproduced by the package.")

## --- discordant set, across the grid ----------------------------------
disc <- do.call(rbind, lapply(ALPHA_GRID[ALPHA_GRID <= 0.05], function(a) {
  d <- diagFDR::dfdr_entrapment_discordant(diagFDR::dfdr_accepted(U$global, a),
                                           diagFDR::dfdr_accepted(U$minrunq, a), ent_prec)
  data.frame(alpha = a, n_core = d$n[d$set == "core"], E_core = d$E[d$set == "core"],
             fdp_core = d$fdp[d$set == "core"], n_added = d$n[d$set == "added"],
             E_added = d$E[d$set == "added"], fdp_added = d$fdp[d$set == "added"],
             lo_added = d$lo[d$set == "added"], hi_added = d$hi[d$set == "added"],
             n_removed = d$n[d$set == "removed"], ratio = attr(d, "ratio"),
             jaccard = attr(d, "jaccard"), nested = attr(d, "nested"))
}))
print(disc); write_tab(disc, "entrapment_discordant_package")
d1 <- disc[disc$alpha == ALPHA_MAIN, ]
stopifnot(d1$n_added == 6421, d1$E_added == 165, d1$n_removed == 0, abs(d1$ratio - 11.9) < 0.1)
say("Equation (5) discordant/core ratio at 1%%: %.1f (Jaccard %.2f)", d1$ratio, d1$jaccard)

## --- entrapment curve for Figure 3C/D ----------------------------------
curve <- rbind(cbind(universe = "global", diagFDR::dfdr_entrapment(U$global, ent_prec, alpha = ALPHA_GRID, r = r_pep)),
               cbind(universe = "minrunq", diagFDR::dfdr_entrapment(U$minrunq, ent_prec, alpha = ALPHA_GRID, r = r_pep)))
write_tab(curve, "entrapment_curve_package")
