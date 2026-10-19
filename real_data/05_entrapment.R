## =====================================================================
##  05 --- Entrapment: FDP, the size ratio r, the corrected estimator
##
##  Produces the entrapment results of Section 4.2 and Note S4.
##
##  ATTRIBUTION. Protein.Group in this export holds BARE UniProt
##  accessions (P55011, Q96JP5, ...) with no species tag, so a string
##  prefix such as "ARATH" matches nothing at all. An identification is
##  counted towards E_alpha when EVERY accession in its Protein.Group is
##  present in the A. thaliana FASTA ("pure entrapment").
##
##  Expected on the Zenodo data at alpha = 1%:
##    global Global.Q.Value   T = 81,643  E = 177  FDP = 0.217% [0.187, 0.251]
##    global min_r(Q.Value)   T = 88,064  E = 342  FDP = 0.388% [0.349, 0.432]
##    fold increase 1.79
## =====================================================================

source("00_config.R")
load_numbers()

## =====================================================================
##  PART A --- size ratio r and homology. Needs only the FASTA files.
## =====================================================================
message("\n=== A. entrapment database characterisation ===")
r_pep_unique <- NA_real_

if (file.exists(FASTA_TARGET) && file.exists(FASTA_ENTRAP)) {
  tgt <- read_fasta(FASTA_TARGET)
  ent <- read_fasta(FASTA_ENTRAP)
  con <- if (file.exists(FASTA_CONTAM)) read_fasta(FASTA_CONTAM)
         else list(ids = character(0), seqs = character(0))

  n_tgt_prot <- length(tgt$ids) + length(con$ids)
  n_ent_prot <- length(ent$ids)
  say("target proteins (incl. contaminants): %s", format(n_tgt_prot, big.mark = ","))
  say("entrapment proteins:                  %s", format(n_ent_prot, big.mark = ","))

  message("  digesting (a few minutes on the full proteomes)...")
  pep_tgt <- unique(unlist(lapply(c(tgt$seqs, con$seqs), digest)))
  pep_ent <- unique(unlist(lapply(ent$seqs, digest)))
  shared  <- intersect(pep_tgt, pep_ent)

  r_prot <- n_ent_prot / n_tgt_prot
  r_pep  <- length(pep_ent) / length(pep_tgt)
  r_pep_unique <- length(setdiff(pep_ent, shared)) / length(pep_tgt)

  say("tryptic peptides target/entrapment/shared: %s / %s / %s",
      format(length(pep_tgt), big.mark = ","), format(length(pep_ent), big.mark = ","),
      format(length(shared), big.mark = ","))
  say("r protein = %.3f | r peptide = %.3f | r peptide (shared removed) = %.3f",
      r_prot, r_pep, r_pep_unique)

  reg("rProtein", r_prot, "d2", "entrapment-to-target size ratio, protein level")
  reg("rPeptide", r_pep,  "d2", "peptide level -- the one that matters")
  reg("rPeptideUnique", r_pep_unique, "d2", "after removing shared peptides")
  reg("sharedPeptides", length(shared), "int", "homology overlap")
  reg("entrapProteins", n_ent_prot, "int", "entrapment database size")

  write_tab(data.frame(
    quantity = c("target proteins","entrapment proteins","target peptides",
                 "entrapment peptides","shared peptides","r protein","r peptide",
                 "r peptide unique"),
    value = c(n_tgt_prot, n_ent_prot, length(pep_tgt), length(pep_ent),
              length(shared), r_prot, r_pep, r_pep_unique)), "entrapment_composition")

  cat("\n--- database composition (Section 3.5) ---\n")
  cat(sprintf(
"so we report $r$ and the entrapment database composition ($r=%.2f$ at protein\nlevel, $%.2f$ at tryptic-peptide level after removing the $%s$ peptides shared\nwith the target proteomes; \\emph{A.~thaliana}, $%s$ entries) and give the\ncorrected estimator alongside the raw ratio.\n\n",
    r_prot, r_pep_unique, formatC(length(shared), format = "d", big.mark = "{,}"),
    formatC(n_ent_prot, format = "d", big.mark = "{,}")))
} else message("  FASTA files not found -- Part A skipped.")

## =====================================================================
##  PART B --- entrapment FDP
## =====================================================================
message("\n=== B. entrapment FDP ===")
U <- readRDS(file.path(OUT, "universes.rds"))

if (is.null(U$entrap)) {
  message("  No entrapment search imported; Section 4.2's entrapment paragraph")
  message("  cannot be regenerated. Point DIANN_ENTRAP at the report and rerun 02.")
} else {
  need_file(FASTA_ENTRAP, "entrapment FASTA")
  ACC <- entrap_accessions(FASTA_ENTRAP)
  say("%s entrapment accessions", format(length(ACC), big.mark = ","))

  pm  <- U$entrap$prot_map
  col <- if (ENTRAP_COL == "protein_ids") "protein_ids" else "protein"
  prot_of <- setNames(pm[[col]], pm$id)

  ## pure entrapment: every accession in the group is an entrapment accession
  is_pure <- function(ids) {
    s <- prot_of[ids]
    vapply(strsplit(ifelse(is.na(s), "", s), ";"), function(z) {
      z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) }, logical(1))
  }
  ## any entrapment: at least one accession is an entrapment accession
  is_any <- function(ids) {
    s <- prot_of[ids]
    vapply(strsplit(ifelse(is.na(s), "", s), ";"), function(z) any(z %in% ACC), logical(1))
  }

  ent_fdp <- function(tbl, alpha = ALPHA_MAIN) {
    acc <- tbl$id[!tbl$is_decoy & tbl$q <= alpha]
    ## run x precursor ids carry a "run||" prefix; strip for the lookup
    acc_base <- sub("^.*\\|\\|", "", acc)
    Ta <- length(acc_base)
    if (!Ta) return(data.frame(T_alpha = 0, E_alpha = 0, E_any = 0,
                               fdp = NA, lo = NA, hi = NA))
    Ea <- sum(is_pure(acc_base)); Ay <- sum(is_any(acc_base))
    w  <- wilson(Ea, Ta)
    data.frame(T_alpha = Ta, E_alpha = Ea, E_any = Ay,
               fdp = unname(w["est"]), lo = unname(w["lo"]), hi = unname(w["hi"]))
  }

  rows <- do.call(rbind, lapply(c("global","minrunq","runxprec"), function(u)
    cbind(universe = u, ent_fdp(U$entrap[[u]]))))
  print(rows); write_tab(rows, "entrapment_fdp")

  if (all(rows$E_alpha == 0))
    warning("Every E_alpha is zero. Check ENTRAP_COL and that the entrapment ",
            "search really used this FASTA. Run 01_preflight.R section 5.",
            call. = FALSE)

  ge <- function(u, c) { v <- rows[rows$universe == u, c]; if (!length(v)) NA else v[1] }
  reg("entrapGlobal",   ge("global","fdp"),  "d3", "Sec 4.2 entrapment FDP, correct global filtering")
  reg("entrapGlobalLo", ge("global","lo"),   "d3", "Wilson 95% CI")
  reg("entrapGlobalHi", ge("global","hi"),   "d3", "Wilson 95% CI")
  reg("entrapMinrun",   ge("minrunq","fdp"), "d3", "Sec 4.2, min_r aggregation")
  reg("entrapMinrunLo", ge("minrunq","lo"),  "d3", "Wilson 95% CI")
  reg("entrapMinrunHi", ge("minrunq","hi"),  "d3", "Wilson 95% CI")
  reg("entrapRunwise",  ge("runxprec","fdp"),"d3", "Sec 4.2, run-wise on its native universe")
  reg("entrapFold",     ge("minrunq","fdp") / ge("global","fdp"), "d2", "fold increase, Sec 4.2")
  reg("TentrapGlobal",  ge("global","T_alpha"), "int", "Sec 4.2")
  reg("EentrapGlobal",  ge("global","E_alpha"), "int", "Sec 4.2")
  reg("TentrapMinrun",  ge("minrunq","T_alpha"), "int", "Sec 4.2")
  reg("EentrapMinrun",  ge("minrunq","E_alpha"), "int", "Sec 4.2")

  ## --- Size-corrected estimator --------------------------------------
  ## Two estimators are in use and they estimate different quantities.
  ##   combined     E (1 + 1/r) / (T + E)   FDP of the target-plus-
  ##                                        entrapment list; Equation (4)
  ##                                        of the paper, reported there.
  ##   target-only  (E / T) / r             FDP of the target-only list.
  ## Both are monotone in E/T at fixed r, so comparisons between lists
  ## are the same under either; only the magnitudes differ (Note S4.2).
  if (is.finite(r_pep_unique)) {
    rows$fdp_corrected <- rows$E_alpha * (1 + 1 / r_pep_unique) /
                          (rows$T_alpha + rows$E_alpha)
    rows$fdp_corrected_target_only <-
                          (rows$E_alpha / r_pep_unique) / rows$T_alpha
    print(rows[, c("universe","fdp","fdp_corrected","fdp_corrected_target_only")])
    write_tab(rows, "entrapment_fdp_corrected")
    reg("entrapGlobalCorr", rows$fdp_corrected[rows$universe == "global"], "d3",
        "Equation (4), combined estimator, peptide-level r")
    reg("entrapMinrunCorr", rows$fdp_corrected[rows$universe == "minrunq"], "d3",
        "Equation (4), combined estimator")
    reg("entrapGlobalCorrTargetOnly",
        rows$fdp_corrected_target_only[rows$universe == "global"], "d3",
        "Note S4.2, target-only estimator")
  } else message("  r unavailable (Part A skipped) -- corrected estimator not computed.")

  ## --- protein-level entrapment --------------------------------
  if (!is.null(U$entrap$pg_global)) {
    pgt <- U$entrap$pg_global
    accp <- pgt$id[!pgt$is_decoy & pgt$q <= ALPHA_MAIN]
    Ep <- sum(vapply(strsplit(accp, ";"), function(z) {
      z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) }, logical(1)))
    wp <- wilson(Ep, length(accp))
    say("protein-level entrapment FDP = %.3f%% (T=%d, E=%d)",
        100 * wp["est"], length(accp), Ep)
    reg("entrapProteinLevel", unname(wp["est"]), "d3", "protein-level entrapment FDP")
    reg("TentrapProtein", length(accp), "int", "Sec 4.2")
    reg("EentrapProtein", Ep, "int", "Sec 4.2")
    write_tab(data.frame(universe = "pg_global", T_alpha = length(accp), E_alpha = Ep,
                         fdp = unname(wp["est"]), lo = unname(wp["lo"]), hi = unname(wp["hi"])),
              "entrapment_protein_level")
  } else {
    message("\n  No protein-group universe in the entrapment search.")
  }

  ## --- Figures 5C and 5D ----------------------------------------------
  p <- open_pdf("fig5C_entrapment", 6.2, 4.4)
  m <- rbind("decoy ratio" = 100 * vapply(c("global","minrunq","runxprec"),
               function(u) headline_of(U$entrap[[u]])$decoy_ratio, numeric(1)),
             "entrapment FDP" = 100 * rows$fdp)
  barplot(m, beside = TRUE, ylab = "%", names.arg = c("global","min_r","run-wise"),
          legend.text = rownames(m), args.legend = list(x = "topleft", bty = "n"))
  abline(h = 100 * ALPHA_MAIN, lty = 2); dev.off(); message("  -> ", p)

  gr <- do.call(rbind, lapply(ALPHA_GRID, function(a) {
    d <- ent_fdp(U$entrap$global, a); h <- headline_of(U$entrap$global, a)
    data.frame(alpha = a, entrap = d$fdp, ratio = h$decoy_ratio) }))
  write_tab(gr, "entrapment_by_alpha")
  yl <- 100 * c(gr$ratio, gr$entrap); yl <- yl[is.finite(yl) & yl > 0]
  p <- open_pdf("fig5D_entrapment_stability", 6.2, 4.4)
  plot(gr$alpha, pmax(100 * gr$ratio, min(yl)/2), type = "b", pch = 19, log = "xy",
       xlab = expression(alpha), ylab = "%", ylim = range(yl))
  lines(gr$alpha, pmax(100 * gr$entrap, min(yl)/2), type = "b", pch = 17, col = 2)
  legend("topleft", c("decoy ratio","entrapment FDP"), col = 1:2, pch = c(19,17), bty = "n")
  dev.off(); message("  -> ", p)
}

save_numbers()
message("\nNext: 06_msrescore.R")
