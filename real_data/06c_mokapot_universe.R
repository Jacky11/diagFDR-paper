## =====================================================================
##  06c --- On which universe are mokapot's q-values computed?
##          (main-text Section 4.5, Note S8)
##
##  The mokapot export carries every candidate PSM it received: 31,314
##  rows, up to ten candidates per spectrum (MS2Rescore input rank <= 10).
##  Its q-values are target-decoy counts, (D + 1) / T with a running
##  minimum, computed over THAT universe. Thresholding them on the
##  10,568 spectrum-level winners therefore applies a q-value defined on
##  one universe to another one.
##
##  Expected on the Zenodo data at alpha = 1%:
##    all candidate PSMs, mokapot q       : T = 7,987  D = 78  (0.98%)
##    max |mokapot q - (D+1)/T recount|    : < 6e-5
##    winners, mokapot q (as exported)    : T = 7,070  D = 43  (0.61%)
##    winners, (D+1)/T recount on score   : T = 7,164  D = 70  (0.98%)
##  The last count is the "7164 more PSMs ... at 0.01 FDR" of the
##  MS2Rescore log.
##
##  Output: tables/mokapot_universe.csv
## =====================================================================

source("00_config.R")
load_numbers()
need_file(MOKAPOT_TARGETS, "mokapot target PSMs")
need_file(MOKAPOT_DECOYS,  "mokapot decoy PSMs")

raw <- diagFDR::read_mokapot_psms(target_path = MOKAPOT_TARGETS,
                                  decoy_path  = MOKAPOT_DECOYS)
dec <- !diagFDR:::mokapot_is_target_to_logical(raw$is_target)
sc  <- as.numeric(raw[["mokapot score"]]); qm <- as.numeric(raw[["mokapot q-value"]])

q_all <- tdc_qvalues(sc, dec, offset = 1L)        # (D+1)/T on all candidates
maxdiff <- max(abs(q_all - qm))
say("max |mokapot q - (D+1)/T over all %s PSMs| = %.2g",
    format(nrow(raw), big.mark = ","), maxdiff)

## spectrum-level winners, as in 06_msrescore.R
win <- diagFDR::mokapot_competed_universe(
  raw, run_col = "run", spectrum_id_col = "spectrum_id",
  score_col = "mokapot score", q_col = "mokapot q-value", pep_col = "mokapot PEP",
  id_mode = "runxid", unit = "psm", scope = "global", q_source = "mokapot q-value")
q_win <- tdc_qvalues(win$score, win$is_decoy, offset = 1L)

cnt <- function(q, d, a = ALPHA_MAIN) c(T = sum(!d & q <= a), D = sum(d & q <= a))
res <- rbind(
  data.frame(universe = "all candidate PSMs", q = "mokapot q (exported)", t(cnt(qm, dec))),
  data.frame(universe = "spectrum-level winners", q = "mokapot q (exported)", t(cnt(win$q, win$is_decoy))),
  data.frame(universe = "spectrum-level winners", q = "(D+1)/T recount on mokapot score", t(cnt(q_win, win$is_decoy))))
res$decoy_ratio <- res$D / res$T
res$n_universe <- c(nrow(raw), nrow(win), nrow(win))
print(res)
write_tab(res, "mokapot_universe")

reg("mokapotQmaxdiff", maxdiff, NULL, "Sec 4.5: mokapot q reproduced by (D+1)/T on all candidates")
reg("mokapotAllT", res$T[1], "int", "Note S8"); reg("mokapotAllD", res$D[1], "int", "Note S8")
reg("mokapotNall", nrow(raw), "int", "Sec 4.5: candidate PSMs")
save_numbers()
