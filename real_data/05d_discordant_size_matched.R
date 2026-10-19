## =====================================================================
##  05d --- Is the scope-discordant set more than the margin of a
##          looser list?  (Note S4, main-text Section 4.2)
##
##  The discordant set of 05b is the set of precursors that the min_r
##  list adds to the Global.Q.Value list. Any list relaxed beyond the
##  global 1% list would also add a margin whose local error exceeds that
##  of the core, so the 11.9-fold ratio must be compared with a
##  size-matched relaxation of the correctly scoped list: the global
##  list cut at the size of the min_r list.
##
##  On these data the two columns rank target precursors almost
##  identically, and at alpha = 1% (and 5%) the min_r list IS the
##  Global.Q.Value list at a looser threshold:
##    entrapment search : 2.19% (10.1% at alpha = 5%)
##    basic search      : 2.23% (10.6%)
##    semi-tryptic      : 2.08% ( 9.8%)
##  so the size-matched control accepts the same precursors, and the
##  discordant-set rate measures the error of the margin added by the
##  mislabelled threshold.
##
##  Expected on the Zenodo data (entrapment search):
##    alpha = 0.01: |min_r| = 88,064, matched global threshold 0.02188,
##                  Jaccard(min_r, matched global) = 1.000,
##                  added E = 165 / 6,421 under both constructions
##    Spearman(Global.Q.Value, min_r Q.Value) over targets ~ 0.999
##
##  Uses out/universes.rds (02_import_diann.R) and the entrapment FASTA.
##  Output: tables/discordant_size_matched.csv
## =====================================================================

source("00_config.R")
load_numbers()
need_file(FASTA_ENTRAP, "entrapment FASTA")
U <- readRDS(file.path(OUT, "universes.rds"))
if (is.null(U$entrap)) stop("universes.rds has no entrapment universes; run 02 with the entrapment search.", call. = FALSE)

ACC <- entrap_accessions(FASTA_ENTRAP)
pm  <- U$entrap$prot_map
prot <- as.character(pm[[ENTRAP_COL]])
pure <- vapply(strsplit(ifelse(is.na(prot), "", prot), ";"),
               function(z) { z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) },
               logical(1))
names(pure) <- pm$id
if (!any(pure)) stop("No pure-entrapment precursor: attribution is wrong.", call. = FALSE)

## targets only, on a common identifier
g <- U$entrap$global;  g <- g[!g$is_decoy, c("id", "q")]
m <- U$entrap$minrunq; m <- m[!m$is_decoy, c("id", "q")]
qg <- setNames(g$q, g$id); qm <- setNames(m$q, m$id)
ids <- intersect(names(qg), names(qm))
rho_s <- suppressWarnings(cor(qg[ids], qm[ids], method = "spearman"))
say("Spearman(Global.Q.Value, min_r Q.Value) over target precursors = %.4f", rho_s)

rate <- function(v) {
  if (!length(v)) return(c(n = 0, E = 0, est = NA, lo = NA, hi = NA))
  E <- sum(pure[v], na.rm = TRUE); w <- wilson(E, length(v))
  c(n = length(v), E = E, est = unname(w["est"]), lo = unname(w["lo"]), hi = unname(w["hi"]))
}

og <- names(qg)[order(qg, names(qg))]
rows <- lapply(c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2), function(a) {
  G <- names(qg)[qg <= a]; M <- names(qm)[qm <= a]
  Gs <- og[seq_len(length(M))]                     # global list cut at |M|
  added_minr <- setdiff(M, G); added_glob <- setdiff(Gs, G)
  r1 <- rate(added_minr); r2 <- rate(added_glob); rc <- rate(intersect(G, M))
  data.frame(alpha = a, n_minr = length(M), q_global_matched = max(qg[Gs]),
             J_minr_vs_matched = length(intersect(M, Gs)) / length(union(M, Gs)),
             rate_core = rc["est"],
             n_added_minr = r1["n"], E_added_minr = r1["E"], rate_added_minr = r1["est"],
             n_added_global = r2["n"], E_added_global = r2["E"], rate_added_global = r2["est"],
             lo_added_global = r2["lo"], hi_added_global = r2["hi"], row.names = NULL)
})
res <- do.call(rbind, rows)
print(res, digits = 4)
write_tab(res, "discordant_size_matched")

main <- res[which.min(abs(res$alpha - ALPHA_MAIN)), ]
reg("spearmanGlobalMinr", rho_s, "d3", "Sec 4.2 / Note S4: rank agreement of the two columns")
reg("qGlobalMatched", main$q_global_matched, "pct2", "Sec 4.2: Global.Q.Value threshold equivalent to min_r at 1% (entrapment search)")
reg("JminrMatched", main$J_minr_vs_matched, "d3", "Note S4: Jaccard(min_r list, size-matched global list)")
save_numbers()
message("\nNext: 06_msrescore.R")
