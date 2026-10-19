## =====================================================================
##  06b --- Before/after list overlap for MS2Rescore, on a COMMON KEY
##
##  WHY THIS SCRIPT EXISTS
##
##  06_msrescore.R builds the pre-rescoring universe with
##  diagFDR::read_dfdr_maxquant_msms(), which keys a PSM on the MaxQuant
##  integer `id`. The mokapot export keys a PSM on run + spectrum_id.
##  Those are different identifier spaces, so a Jaccard index computed
##  between the two accepted sets measures the key mismatch and not the
##  change in identifications -- it returns 0.05-0.10 whatever the data.
##
##  The mapping is in fact exact: MS2Rescore was run with
##      spectrum_id_pattern = ".*scan=(\\d+)$"
##  so mokapot's `spectrum_id` IS the MaxQuant `Scan number`, and
##  mokapot's `run` IS the MaxQuant `Raw file`. Verified on this dataset:
##  the two run names are identical and competing msms.txt to one winner
##  per (Raw file, Scan number) gives exactly 10,568 spectra, the same
##  number MS2Rescore reports after rank-1 filtering.
##
##  This script rebuilds the pre-rescoring universe on that key, competes
##  it the same way the post-rescoring universes are competed, and writes
##  the overlap curve and Figure 5E.
##
##  The before/after counts of Section 4.5 are given on the COMPETED
##  universes (10,568 spectrum-level winners on both sides):
##    before rescoring  T = 5,610, D = 55 at 1%;  T = 2,138, D = 0 at 5e-4
##    mokapot q-values  T = 7,070, D = 43 at 1%;  T = 5,958 at 5e-4
##    TDC on rescored   T = 7,164, D = 70 at 1%;  T = 6,013 at 5e-4
##  Figure 5D draws the pre-rescoring curve on the uncompeted PSM-level
##  universe of 06_msrescore.R (T = 5,481, D = 53 at 1%; T = 2,139, D = 0
##  at 5e-4); both are checked below.
##
##  Output: out/tables/msrescore_overlap_scankey.csv
##          out/figures/fig6E_overlap.pdf
## =====================================================================

source("00_config.R")

need_file(MAXQUANT_MSMS,        "MaxQuant msms.txt")
need_file(MOKAPOT_TARGETS, "mokapot target PSMs")
need_file(MOKAPOT_DECOYS,  "mokapot decoy PSMs")

## --- helpers ---------------------------------------------------------
## One winner per key, by highest score. Ties resolved by first occurrence,
## which matches the order MS2Rescore emits.
compete_by <- function(key, score) {
  o <- order(key, -score, method = "radix")
  key <- key[o]
  idx <- o[!duplicated(key)]
  sort(idx)
}

## --- (1) pre-rescoring universe from msms.txt ------------------------
mq <- readr::read_tsv(MAXQUANT_MSMS, show_col_types = FALSE, progress = FALSE,
                      na = c("", "NA"))
req <- c("Raw file", "Scan number", "Score", "Reverse")
if (length(setdiff(req, names(mq))))
  stop("msms.txt is missing: ", paste(setdiff(req, names(mq)), collapse = ", "))

## read_dfdr_maxquant_msms() drops potential contaminants when the column
## exists; do the same so the two arms are comparable.
if ("Potential contaminant" %in% names(mq))
  mq <- mq[is.na(mq[["Potential contaminant"]]) |
             mq[["Potential contaminant"]] != "+", ]

mq <- mq[is.finite(suppressWarnings(as.numeric(mq$Score))), ]
mq_key <- paste0(mq[["Raw file"]], "||", mq[["Scan number"]])
mq_dec <- !is.na(mq$Reverse) & mq$Reverse == "+"
mq_sc  <- as.numeric(mq$Score)

say("msms.txt: %s rows, %s decoys", format(nrow(mq), big.mark = ","),
    format(sum(mq_dec), big.mark = ","))

## (a) uncompeted, as 06_msrescore.R builds it -- reproduces Sec 4.5
q_unc <- tdc_qvalues(mq_sc, mq_dec, offset = 1)
## (b) competed to one winner per spectrum -- used for the overlap
i  <- compete_by(mq_key, mq_sc)
bk <- mq_key[i]; bd <- mq_dec[i]
q_cmp <- tdc_qvalues(mq_sc[i], bd, offset = 1)
say("competed pre-rescoring universe: %s spectra, %s decoys",
    format(length(bk), big.mark = ","), format(sum(bd), big.mark = ","))

## --- (2) post-rescoring universes ------------------------------------
rd <- function(p) readr::read_tsv(p, show_col_types = FALSE, progress = FALSE)
mt <- rd(MOKAPOT_TARGETS); md <- rd(MOKAPOT_DECOYS)
ak <- c(paste0(mt$run, "||", mt$spectrum_id), paste0(md$run, "||", md$spectrum_id))
as_ <- c(mt[["mokapot score"]],   md[["mokapot score"]])
aq <- c(mt[["mokapot q-value"]], md[["mokapot q-value"]])
ad <- c(rep(FALSE, nrow(mt)),    rep(TRUE, nrow(md)))

j  <- compete_by(ak, as_)
ak <- ak[j]; aq <- aq[j]; ad <- ad[j]; as_ <- as_[j]
q_res <- tdc_qvalues(as_, ad, offset = 1)   # TDC on the rescored score
say("competed post-rescoring universe: %s spectra, %s decoys",
    format(length(ak), big.mark = ","), format(sum(ad), big.mark = ","))

if (length(ak) != length(bk))
  warning("the two competed universes differ in size (", length(bk), " vs ",
          length(ak), "); the key mapping may be wrong for this dataset.",
          call. = FALSE)
ov <- length(intersect(bk, ak)) / length(union(bk, ak))
say("spectrum-key agreement between the two universes: %.4f", ov)
if (ov < 0.95)
  warning("the two universes barely share identifiers (Jaccard ", round(ov, 3),
          "). Check that mokapot's spectrum_id really is the MaxQuant scan ",
          "number for this run -- it depends on MS2Rescore's ",
          "spectrum_id_pattern.", call. = FALSE)

## --- (3) overlap across the threshold grid ---------------------------
jac <- function(a, b) if (!length(a) && !length(b)) NA_real_ else
  length(intersect(a, b)) / length(union(a, b))

res <- do.call(rbind, lapply(ALPHA_GRID, function(a) {
  B  <- bk[q_cmp <= a & !bd]
  Mk <- ak[aq    <= a & !ad]
  Tc <- ak[q_res <= a & !ad]
  data.frame(alpha = a,
             T_before_uncompeted = sum(q_unc <= a & !mq_dec),
             D_before_uncompeted = sum(q_unc <= a &  mq_dec),
             T_before_competed   = length(B),
             D_before_competed   = sum(q_cmp <= a & bd),
             T_mokapot           = length(Mk),
             D_mokapot           = sum(aq <= a & ad),
             T_rescoredTDC       = length(Tc),
             D_rescoredTDC       = sum(q_res <= a & ad),
             J_before_mokapot    = jac(B, Mk),
             J_before_rescoredTDC = jac(B, Tc))
}))
print(res, digits = 4)
write_tab(res, "msrescore_overlap_scankey")

## --- (4) consistency check against the paper --------------------
chk <- function(a, col, want, what) {
  got <- res[[col]][which.min(abs(res$alpha - a))]
  ok  <- isTRUE(all.equal(as.numeric(got), as.numeric(want)))
  message(sprintf("  [%s] %-34s %s (paper: %s)",
                  if (ok) "ok" else "MISMATCH", what, format(got, big.mark = ","),
                  format(want, big.mark = ",")))
  ok
}
message("\n--- reproducing the Section 4.5 counts ---")
all_ok <- c(
  chk(0.01,  "T_before_competed",   5610, "T before (winners), alpha = 1%"),
  chk(0.01,  "D_before_competed",     55, "D before (winners), alpha = 1%"),
  chk(5e-4,  "T_before_competed",   2138, "T before (winners), alpha = 5e-4"),
  chk(5e-4,  "D_before_competed",      0, "D before (winners), alpha = 5e-4"),
  chk(0.01,  "D_mokapot",             43, "D mokapot, alpha = 1%"),
  chk(0.01,  "D_rescoredTDC",         70, "D rescored TDC, alpha = 1%"),
  chk(0.01,  "T_before_uncompeted", 5481, "T before (all PSMs, Fig 5D), alpha = 1%"),
  chk(0.01,  "T_mokapot",           7070, "T mokapot, alpha = 1%"),
  chk(0.01,  "T_rescoredTDC",       7164, "T rescored TDC, alpha = 1%"),
    chk(5e-4,  "T_before_uncompeted", 2139, "T before (all PSMs, Fig 5D), alpha = 5e-4"),
  chk(5e-4,  "T_mokapot",           5958, "T mokapot, alpha = 5e-4"),
  chk(5e-4,  "T_rescoredTDC",       6013, "T rescored TDC, alpha = 5e-4"))
if (!all(all_ok))
    warning("at least one count of Section 4.5 did not reproduce.", call. = FALSE)

## --- (5) Figure 5E ---------------------------------------------------
pdf(file.path(FIGDIR, "fig6E_overlap.pdf"), width = 4.2, height = 3.6)
op <- par(mar = c(4.2, 4.2, 2.2, 1))
plot(res$alpha, res$J_before_mokapot, type = "b", pch = 16, log = "x",
     ylim = c(0, 1), xlab = expression(alpha),
     ylab = "Jaccard vs before rescoring", main = "E")
lines(res$alpha, res$J_before_rescoredTDC, type = "b", pch = 1, lty = 2)
abline(v = 0.01, col = "grey70", lty = 3)
legend("bottomright", bty = "n", cex = 0.8, pch = c(16, 1), lty = c(1, 2),
       legend = c("mokapot q-value", "TDC on rescored score"))
par(op); dev.off()
message("\n  -> ", file.path(FIGDIR, "fig6E_overlap.pdf"))

## --- (6) registered numbers ----------------------------------
g <- function(a, col) res[[col]][which.min(abs(res$alpha - a))]
reg("JrescoreStrict", g(5e-4, "J_before_mokapot"), "d2",
    "Sec 4.5 / Fig 5E: before-vs-after list overlap at 5e-4, common scan key")
reg("JrescoreOne",    g(0.01, "J_before_mokapot"), "d2",
    "Sec 4.5 / Fig 5E: before-vs-after list overlap at 1%, common scan key")
reg("JrescoreFive",   g(0.05, "J_before_mokapot"), "d2",
    "Sec 4.5 / Fig 5E: before-vs-after list overlap at 5%, common scan key")
reg("nSpectraRescore", length(bk), "int",
    "Sec 4.5 / Methods: spectra in each competed rescoring universe")
write_numbers_tex()
message("Next: 08_defaults_sensitivity.R")
