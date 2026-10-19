## =====================================================================
##  05b --- The discordance-restricted entrapment rate
##
##  Equation (5) of the paper.
##
##  THE PROBLEM. J(run, global) = 0.93 at alpha = 1% -- high agreement on
##  any conventional reading -- while the entrapment FDP nearly doubles
##  between the same two lists. A Jaccard index is therefore not a proxy
##  for error impact.
##
##  THE MEASURE. Split the union of the two accepted lists into the CORE
##  (accepted by both) and the DISCORDANT set (accepted by one only), and
##  report the entrapment rate -- or, where entrapment is unavailable, the
##  decoy ratio -- separately in each. It uses the same two lists, costs
##  nothing extra, and unlike J it responds to the quantity at issue.
##
##  Expected on the Zenodo data at alpha = 1%:
##    J = 0.927
##    core     n = 81,643  E = 177  rate 0.217%
##    added    n =  6,421  E = 165  rate 2.570% [2.210, 2.986]
##    fold = 11.9
##  and 4.2x-11.9x across alpha in [1e-3, 5e-2] while J stays in
##  [0.91, 0.94].
##
##  NOTE ON THE COMPARISON. At alpha = 1% the global list is a strict
##  SUBSET of the min_r list, so "discordant" is entirely the set min_r
##  adds. The script does not assume this -- it reports both directions --
##  but the drop set is empty here and prints as such.
##
##  Outputs, in tables/:
##    discordant_by_alpha.csv   -> Section 4.2
##    fig5B_discordance.pdf     discordant-set rate against alpha
## =====================================================================

source("00_config.R")
load_numbers()

if (!file.exists(DIANN_ENTRAP))
  stop("This script needs the entrapment search. Point DIANN_ENTRAP at it.",
       call. = FALSE)
need_file(FASTA_ENTRAP, "entrapment FASTA")

suppressPackageStartupMessages(library(dplyr))

ACC <- entrap_accessions(FASTA_ENTRAP)
say("%s entrapment accessions", format(length(ACC), big.mark = ","))

rep_ent <- sanitise_report(diagFDR::read_diann_parquet(DIANN_ENTRAP),
                           COL$decoy, "entrapment")

## precursor -> Protein.Group map, and the pure-entrapment indicator.
## Computed ONCE: it is the same map for every alpha.
pm <- rep_ent[!duplicated(rep_ent[[COL$precursor]]), ]
prot_of <- setNames(as.character(pm[[COL$protein]]),
                    as.character(pm[[COL$precursor]]))
is_pure_vec <- vapply(strsplit(ifelse(is.na(prot_of), "", prot_of), ";"),
                      function(z) { z <- z[nzchar(z)]
                                    length(z) > 0 && all(z %in% ACC) },
                      logical(1))
names(is_pure_vec) <- names(prot_of)

## GUARD. Protein.Group holds BARE accessions with no species tag. If the
## attribution is wrong every count is silently zero, so assert before
## computing anything.
if (!any(is_pure_vec))
  stop("No precursor maps to a pure-entrapment protein group. The ",
       "attribution is wrong -- check ENTRAP_COL and that this search ",
       "really used ", basename(FASTA_ENTRAP), ".", call. = FALSE)

## Two experiment-wide universes on a common identifier (Precursor.Id):
##   global : Global.Q.Value, the correctly scoped list
##   minrun : min over runs of the run-wise Q.Value, the scope misuse
U <- rep_ent %>%
  group_by(id = as.character(.data[[COL$precursor]])) %>%
  summarise(qg = min(.data[[COL$q_global]], na.rm = TRUE),
            qm = min(.data[[COL$q_run]],    na.rm = TRUE),
            is_decoy = any(.data[[COL$decoy]] == 1L),
            .groups = "drop")

erate <- function(ids) {
  n <- length(ids)
  if (!n) return(c(n = 0, E = 0, rate = NA_real_, lo = NA_real_, hi = NA_real_))
  E <- sum(is_pure_vec[ids], na.rm = TRUE)
  w <- wilson(E, n)
  c(n = n, E = E, rate = unname(w["est"]), lo = unname(w["lo"]),
    hi = unname(w["hi"]))
}

rows <- lapply(ALPHA_GRID, function(a) {
  g <- U$id[U$qg <= a & !U$is_decoy]
  m <- U$id[U$qm <= a & !U$is_decoy]
  core <- intersect(g, m); added <- setdiff(m, g); dropped <- setdiff(g, m)
  J <- length(core) / length(union(g, m))
  c1 <- erate(core); c2 <- erate(added); c3 <- erate(dropped)
  data.frame(alpha = a, J = J,
             n_core = c1["n"], E_core = c1["E"], rate_core = c1["rate"],
             n_added = c2["n"], E_added = c2["E"], rate_added = c2["rate"],
             lo_added = c2["lo"], hi_added = c2["hi"],
             n_dropped = c3["n"], rate_dropped = c3["rate"],
             fold = unname(c2["rate"] / c1["rate"]),
             row.names = NULL)
})
res <- do.call(rbind, rows)
print(res[, c("alpha","J","n_core","rate_core","n_added","rate_added","fold")])
write_tab(res, "discordant_by_alpha")

main <- res[which.min(abs(res$alpha - ALPHA_MAIN)), ]
reg("discCore",    main$n_core,     "int",  "Sec 4.2: shared by both scopes")
reg("discCoreRate",main$rate_core,  "d5",   "Sec 4.2: entrapment rate, core")
reg("discAdded",   main$n_added,    "int",  "Sec 4.2: accepted only by min_r")
reg("discAddedRate",main$rate_added,"d5",   "Sec 4.2: entrapment rate, discordant")
reg("discAddedLo", main$lo_added,   "d5",   "Wilson 95% CI")
reg("discAddedHi", main$hi_added,   "d5",   "Wilson 95% CI")
reg("discFold",    main$fold,       "d1",   "Sec 4.2: fold difference")
reg("Jrunglobal",  main$J,          "d3",   "Sec 4.2: Jaccard overlap")

if (main$n_dropped == 0)
  say("the global list is a strict subset of the min_r list at alpha=%.3g", ALPHA_MAIN)

cat("\n--- Section 4.2 ---\n")
cat(sprintf(
"The global list is a strict subset of the $\\min_r$ list here, and the\n$%s$ identifications $\\min_r$ adds carry an entrapment rate of $%.3f\\%%$\n($[%.2f,%.2f]\\%%$) against $%.3f\\%%$ in the $%s$ shared by both --- a factor of\n$%.1f$, and $%.1f$ times nominal.\n\n",
  formatC(main$n_added, format = "d", big.mark = "{,}"),
  100 * main$rate_added, 100 * main$lo_added, 100 * main$hi_added,
  100 * main$rate_core,
  formatC(main$n_core, format = "d", big.mark = "{,}"),
  main$fold, main$rate_added / ALPHA_MAIN))

## --- figure --------------------------------------------------------------
ok <- is.finite(res$rate_core) & is.finite(res$rate_added) &
      res$n_added > 0 & res$alpha >= 1e-4
p <- open_pdf("fig5B_discordance", 6.6, 4.4)
op <- par(mar = c(4.5, 4.5, 2.5, 4.5))
ylim <- c(0, max(100 * res$rate_added[ok], na.rm = TRUE) * 1.15)
bp <- barplot(rbind(100 * res$rate_core[ok], 100 * res$rate_added[ok]),
              beside = TRUE, ylim = ylim,
              names.arg = formatC(res$alpha[ok], format = "g"),
              xlab = expression(alpha), ylab = "entrapment FDP (%)",
              legend.text = c("shared by both scopes",
                              "accepted only by min_r"),
              args.legend = list(x = "topleft", bty = "n", cex = 0.85))
abline(h = 100 * ALPHA_MAIN, lty = 2)
par(new = TRUE)
plot(colMeans(bp), res$J[ok], type = "b", pch = 19, col = "grey35",
     axes = FALSE, xlab = "", ylab = "", ylim = c(0.8, 1))
axis(4); mtext("Jaccard overlap", side = 4, line = 3, cex = 0.9)
par(op); dev.off(); message("  -> ", p)

save_numbers()
message("\nNext: 05c_entrapment_r_sensitivity.R")
