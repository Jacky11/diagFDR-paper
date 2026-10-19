## =====================================================================
##  05 --- Discrimination against decoy support under two operating rules
##
##  The separation delta is raised from 2 to 5 under scenario S0 with
##  Gaussian nulls. The number of accepted decoys is recorded
##    (a) at a fixed nominal level alpha = 1%, and
##    (b) at a fixed sensitivity: the score threshold that accepts 60% of
##        the correct targets.
##  Each point is the mean over B_SIGNFLIP replicates.
##
##  Writes out/tables/discrimination_vs_support.csv and
##         out/figures/figS1_signflip.pdf
## =====================================================================

source("parameters.R")

DELTAS      <- seq(2, 5, by = 0.25)
SENSITIVITY <- 0.6

one_replicate <- function(delta) {
  is_null <- runif(N_UNITS) < PI0
  X <- rnorm(N_UNITS) + ifelse(is_null, 0, delta)
  Y <- rnorm(N_UNITS)
  score <- pmax(X, Y); is_target <- X > Y
  q <- tdc_qvalues(score, is_target)
  threshold <- delta + qnorm(1 - SENSITIVITY)
  c(D_fixed_alpha = sum(q <= ALPHA_MAIN & !is_target),
    D_fixed_sensitivity = sum(score >= threshold & !is_target))
}

rows <- lapply(seq_along(DELTAS), function(j) {
  set.seed(SEED + 8e6 + j)
  m <- rowMeans(replicate(B_SIGNFLIP, one_replicate(DELTAS[j])))
  data.frame(delta = DELTAS[j], auroc = pnorm(DELTAS[j] / sqrt(2)),
             D_fixed_alpha = m[["D_fixed_alpha"]],
             D_fixed_sensitivity = m[["D_fixed_sensitivity"]])
})
signflip <- do.call(rbind, rows)
write_table(signflip, "discrimination_vs_support")

pdf(file.path(OUT_FIGURES, "figS1_signflip.pdf"), width = 9, height = 4.2)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.4, 1))
ylab <- expression("mean accepted decoys  " * D[alpha])
plot(signflip$auroc, signflip$D_fixed_alpha, type = "b", pch = 19, log = "y",
     xlab = "AUROC", ylab = ylab, main = "Fixed nominal level (Prop. S2)")
plot(signflip$auroc, pmax(signflip$D_fixed_sensitivity, 0.01), type = "b",
     pch = 19, log = "y", xlab = "AUROC", ylab = ylab,
     main = "Fixed sensitivity (Prop. S1)")
dev.off()
## embed the fonts (journals reject PDFs with non-embedded fonts); needs Ghostscript
if (nzchar(tools::find_gs_cmd()))
  grDevices::embedFonts(file.path(OUT_FIGURES, "figS1_signflip.pdf"),
                        options = "-dPDFSETTINGS=/prepress -dEmbedAllFonts=true") else
  message("  Ghostscript not found: fonts of figS1_signflip.pdf are not embedded")
message("  wrote ", file.path(OUT_FIGURES, "figS1_signflip.pdf"))
