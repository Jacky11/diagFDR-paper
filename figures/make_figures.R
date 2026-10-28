## =====================================================================
##  diagFDR paper --- Figures 3, 4 and 5, generated from the package
##
##  Every curve and every number drawn here is computed with the exported
##  functions of diagFDR (>= 0.2.1) from the universes written by
##  real_data/02_import_diann.R and real_data/06_msrescore.R.
##
##  Run from the figures/ folder after real_data/run_all.R:
##      Rscript make_figures.R
##  Writes, for each figure, Figure{3,4,5}.pdf (vector), Figure{3,4,5}.png
##  (300 dpi) and Figure{3,4,5}.pptx (one slide, every panel as native
##  PowerPoint shapes through the rvg package, so titles, labels and colours
##  can be edited in PowerPoint). Figures 1 and 2 are drawn by hand
##  (Figure1.pptx, Figure2.pptx); figure2B_report.R regenerates the report
##  captured in Figure 2B.
##
##  Needs: diagFDR >= 0.2.1, ggplot2, patchwork, officer, rvg, scales, ragg.
## =====================================================================

if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE))
  invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))      # Greek letters in the labels
suppressPackageStartupMessages({
  library(diagFDR); library(ggplot2); library(patchwork); library(officer); library(rvg); library(scales)
})
if (utils::packageVersion("diagFDR") < "0.2.1")
  stop("diagFDR >= 0.2.1 is required (found ", utils::packageVersion("diagFDR"), ")")

OUT    <- Sys.getenv("DIAGFDR_OUT", unset = "../real_data/out")
DATA   <- Sys.getenv("DIAGFDR_DATA", unset = "../real_data/data")
FASTA_ENTRAP <- file.path(DATA, "uniprotkb_taxonomy_id_3702_AND_reviewed_2026_02_06.fasta")
FIGOUT <- Sys.getenv("DIAGFDR_FIGOUT", unset = "out")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

ALPHA_GRID <- c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1)
ALPHA_MAIN <- 0.01
EPS <- 0.2

message("diagFDR ", utils::packageVersion("diagFDR"), " | universes from ", OUT)
U  <- readRDS(file.path(OUT, "universes.rds"))

## --- style -------------------------------------------------------------
col_search <- c("Basic" = "#1F77B4", "Semi-tryptic" = "#E69F00")
col_list   <- c("global" = "#1F77B4", "runxprec" = "#0E9F6E", "minrunq" = "#D35400")
lab_list   <- c(global   = "Global list (Global.Q.Value)\n(correct experiment-wide scope)",
                runxprec = "Run×precursor list (Q.Value)\n(correct run-wise scope)",
                minrunq  = "Global precursor list\nfrom Q.Value (min-q over runs)")
lt_pass <- scale_linetype_manual(values = c("First pass" = 2, "Second pass" = 1))
th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10), legend.title = element_blank(),
        legend.background = element_blank(), strip.background = element_rect(fill = "grey95"))
pct_lab <- function(x) paste0(format(100 * x, trim = TRUE, drop0trailing = TRUE), "%")
log_x <- scale_x_log10(breaks = ALPHA_GRID[c(1, 3, 6, 9)], labels = pct_lab)
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

save_figure <- function(p, name, width, height) {
  ggsave(file.path(FIGOUT, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIGOUT, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white",
         device = ragg::agg_png)
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Blank", master = "Office Theme")
  ## one slide the size of the figure; every panel becomes editable shapes
  doc <- ph_with(doc, dml(ggobj = p, fonts = list(sans = "Arial")),
                 location = ph_location(left = 0, top = 0, width = width, height = height))
  target <- file.path(FIGOUT, paste0(name, ".pptx"))
  print(doc, target = target)
  set_slide_size(target, width, height)
  message("  -> ", name, ".pdf / .png / .pptx")
}

## officer writes 10 x 7.5 in slides; resize the slide to the figure so that
## nothing hangs over the edge (edits <p:sldSz> in ppt/presentation.xml).
set_slide_size <- function(pptx, width, height) {
  pptx <- normalizePath(pptx, mustWork = TRUE)
  td <- tempfile("pptx_"); dir.create(td)
  utils::unzip(pptx, exdir = td)
  pres <- file.path(td, "ppt", "presentation.xml")
  x <- readLines(pres, warn = FALSE, encoding = "UTF-8")
  x <- sub('<p:sldSz[^>]*/>', sprintf('<p:sldSz cx="%d" cy="%d"/>', round(width * 914400), round(height * 914400)), x)
  writeLines(x, pres, useBytes = TRUE)
  owd <- setwd(td); on.exit(setwd(owd))
  unlink(pptx)
  utils::zip(pptx, files = list.files(".", recursive = TRUE, all.files = TRUE), flags = "-r9Xq")
  unlink(td, recursive = TRUE)
}

## =====================================================================
##  Figure 3 --- scope audit and external validation (DIA-NN)
## =====================================================================
message("Figure 3")
head_of <- function(x, a = ALPHA_MAIN) dfdr_headline(x, alpha = a)

## A. recomputed decoy ratio per list, basic and semi-tryptic
A <- do.call(rbind, lapply(c("Basic", "Semi-tryptic"), function(s) {
  u <- if (s == "Basic") U$basic else U$semi
  do.call(rbind, lapply(names(lab_list), function(l) {
    h <- head_of(u[[l]])
    data.frame(search = s, list = l, T = h$T_alpha, D = h$D_alpha, ratio = h$FDR_hat)
  }))
}))
A$list <- factor(A$list, levels = names(lab_list), labels = lab_list)
A$search <- factor(A$search, levels = names(col_search))
pA <- ggplot(A, aes(list, 100 * ratio, fill = list)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 100 * ALPHA_MAIN, linetype = 2) +
  geom_text(aes(label = sprintf("T=%s\nD=%s\n%.2f%%", fmt(T), fmt(D), 100 * ratio)),
            vjust = -0.15, size = 2.9, lineheight = 0.9) +
  facet_wrap(~search) +
  scale_fill_manual(values = unname(col_list), guide = "none") +
  scale_y_continuous(limits = c(0, 3.1), expand = expansion(mult = c(0, 0))) +
  labs(title = "A. Scope audit at α = 1%: recomputed decoy ratio",
       subtitle = "Bars show Dα/Tα; dashed line is the nominal 1% claim",
       x = NULL, y = "Decoy ratio (%)") + th +
  theme(axis.text.x = element_text(size = 7.5, angle = 20, hjust = 1), plot.margin = margin(5, 5, 5, 40))

## B. Jaccard overlap between the Global.Q.Value list and the min_r list
jac <- function(u) do.call(rbind, lapply(ALPHA_GRID, function(a) {
  g <- dfdr_accepted(u$global, a); m <- dfdr_accepted(u$minrunq, a)   # accepted target ids
  data.frame(alpha = a, J = length(intersect(g, m)) / length(union(g, m)))
}))
B <- rbind(cbind(search = "Basic", jac(U$basic)), cbind(search = "Semi-tryptic", jac(U$semi)))
B$search <- factor(B$search, levels = names(col_search))
B1 <- B[B$alpha == ALPHA_MAIN, ]
pB <- ggplot(B, aes(alpha, J, colour = search, shape = search)) +
  geom_vline(xintercept = ALPHA_MAIN, linetype = 3) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  geom_point(data = B1, size = 5, shape = 1, stroke = 1.1) +
  geom_text(data = B1, aes(label = sprintf("J=%.3f", J)), hjust = -0.25,
            vjust = c(-0.6, 1.8), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = col_search) + log_x +
  labs(title = "B. Scope-disagreement curve",
       subtitle = "Jaccard overlap: Global.Q.Value list vs min-over-runs Q.Value list",
       x = "α (log scale)", y = "Jaccard overlap") + th +
  theme(legend.position = c(0.82, 0.18))

## C. entrapment search: decoy ratio and raw entrapment FDP per list
ACC <- read_fasta_accessions(FASTA_ENTRAP)
pm  <- U$entrap$prot_map
ent_prec <- pm$id[vapply(strsplit(ifelse(is.na(pm$protein), "", pm$protein), ";"),
                         function(z) { z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) }, logical(1))]
C <- do.call(rbind, lapply(names(lab_list), function(l) {
  e <- dfdr_entrapment(U$entrap[[l]], ent_prec, alpha = ALPHA_MAIN)
  data.frame(list = l, T = e$T_alpha, D = e$D_alpha, E = e$E_alpha,
             `Decoy ratio (D/T)` = e$decoy_ratio, `Entrapment FDP (E/T)` = e$fdp_raw,
             lo = e$lo, hi = e$hi, check.names = FALSE)
}))
Cl <- rbind(data.frame(list = C$list, what = "Entrapment FDP (E/T)", value = C$`Entrapment FDP (E/T)`, lo = C$lo, hi = C$hi),
            data.frame(list = C$list, what = "Decoy ratio (D/T)",    value = C$`Decoy ratio (D/T)`, lo = NA, hi = NA))
Cl$list <- factor(Cl$list, levels = names(lab_list), labels = lab_list)
Cl$what <- factor(Cl$what, levels = c("Entrapment FDP (E/T)", "Decoy ratio (D/T)"))
Clab <- data.frame(list = factor(C$list, levels = names(lab_list), labels = lab_list),
                   y = pmax(C$`Decoy ratio (D/T)`, C$`Entrapment FDP (E/T)`),
                   label = sprintf("T=%s\nD=%s\nE=%s\n%.2f%%", fmt(C$T), fmt(C$D), fmt(C$E), 100 * C$`Decoy ratio (D/T)`))
pC <- ggplot(Cl, aes(list, 100 * value, fill = what)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = position_dodge(width = 0.7), width = 0.15, na.rm = TRUE) +
  geom_hline(yintercept = 100 * ALPHA_MAIN, linetype = 2) +
  geom_text(data = Clab, aes(list, 100 * y, label = label), inherit.aes = FALSE,
            vjust = -0.15, hjust = 0.1, size = 2.9, lineheight = 0.9) +
  geom_text(data = Cl[Cl$what == "Entrapment FDP (E/T)", ],
            aes(label = sprintf("%.2f%%", 100 * value)), position = position_dodge(width = 0.7),
            vjust = -0.4, hjust = 1.1, size = 2.9) +
  scale_fill_manual(values = c("Entrapment FDP (E/T)" = "#F5B183", "Decoy ratio (D/T)" = "#2E86C1")) +
  scale_y_continuous(limits = c(0, 3.3), expand = expansion(mult = c(0, 0))) +
  labs(title = "C. Scope audit at α = 1% with external entrapment",
       subtitle = "Blue: internal decoy ratio (D/T). Orange: external entrapment FDP (E/T) with 95% Wilson interval",
       x = NULL, y = "Rate (%)") + th +
  theme(axis.text.x = element_text(size = 7.5, angle = 20, hjust = 1), legend.position = c(0.3, 0.9),
        legend.direction = "horizontal")

## D. entrapment-based stability of the Global.Q.Value list
D <- do.call(rbind, lapply(ALPHA_GRID, function(a) {
  e0 <- dfdr_entrapment(U$entrap$global, ent_prec, alpha = a)
  e1 <- dfdr_entrapment(U$entrap$global, ent_prec, alpha = (1 + EPS) * a)
  data.frame(alpha = a, dE = abs(e1$E_alpha - e0$E_alpha),
             Dsi = e0$D_alpha^-0.5, Esi = e0$E_alpha^-0.5)
}))
pD1 <- ggplot(D, aes(alpha, dE)) + geom_vline(xintercept = ALPHA_MAIN, linetype = 3) +
  geom_line(colour = "#D35400", linewidth = 0.9) + log_x +
  labs(title = "D. Entrapment-based stability to cutoff perturbation",
       subtitle = sprintf("Change in accepted entrapment matches |E((1+ε)α) − E(α)|, ε = %.1f", EPS),
       x = NULL, y = NULL) + th
Dl <- rbind(data.frame(alpha = D$alpha, index = "Dα^-1/2", value = D$Dsi),
            data.frame(alpha = D$alpha, index = "Eα^-1/2", value = D$Esi))
pD2 <- ggplot(Dl, aes(alpha, value, colour = index, linetype = index)) +
  geom_vline(xintercept = ALPHA_MAIN, linetype = 3) + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("grey25", "#6C5CE7")) + scale_linetype_manual(values = c(1, 2)) + log_x +
  labs(subtitle = "Decoy-support index Dα^-1/2 and entrapment-support index Eα^-1/2",
       x = "α (log scale)", y = NULL) + th + theme(legend.position = c(0.75, 0.75))
pD <- pD1 / pD2

fig3 <- (pA | pC) / (pB | pD) + plot_layout(heights = c(1, 1.1))
save_figure(fig3, "Figure3", 15, 9.5)

## =====================================================================
##  Figure 4 --- calibration and stability diagnostics (DIA-NN)
## =====================================================================
message("Figure 4")
L <- list("Basic|Second pass" = U$basic$global, "Semi-tryptic|Second pass" = U$semi$global)
if (!is.null(U$basic_fp)) L[["Basic|First pass"]] <- U$basic_fp$global
if (!is.null(U$semi_fp))  L[["Semi-tryptic|First pass"]] <- U$semi_fp$global else
  message("  [note] no semi-tryptic first-pass universe (report_semi-first-pass.parquet absent): that series is omitted")
split_key <- function(k) { s <- strsplit(k, "|", fixed = TRUE)[[1]]; data.frame(search = s[1], pass = s[2]) }
series <- function(fun) do.call(rbind, lapply(names(L), function(k) cbind(split_key(k), fun(L[[k]]))))
fix <- function(d) { d$search <- factor(d$search, names(col_search)); d$pass <- factor(d$pass, c("First pass", "Second pass")); d }

st <- fix(series(function(x) { s <- dfdr_curve_stability(x, ALPHA_GRID); data.frame(alpha = s$alpha, D = s$D_alpha, si = s$support_index) }))
sc <- max(st$D) / max(st$si)
pA4 <- ggplot(st, aes(alpha, colour = search, linetype = pass, shape = search)) +
  geom_line(aes(y = D), linewidth = 0.9) + geom_point(aes(y = D), size = 2) +
  geom_line(aes(y = si * sc), linewidth = 0.9) + geom_point(aes(y = si * sc), size = 2) +
  scale_y_continuous(labels = comma, sec.axis = sec_axis(~ . / sc, name = "Dα^-1/2")) +
  scale_colour_manual(values = col_search) + log_x +
  labs(title = "A. Decoy support and support index", subtitle = "Dα (left axis) and Dα^-1/2 (right axis)",
       x = "α", y = "Dα") + lt_pass + th + theme(legend.position = c(0.75, 0.6))

lt <- fix(series(function(x) { w <- dfdr_local_tail(x, ALPHA_GRID, win_rel = EPS, truncation = "drop"); data.frame(alpha = w$alpha, Dwin = w$D_alpha_win) }))
pB4 <- ggplot(lt, aes(alpha, pmax(Dwin, 1), colour = search, linetype = pass, shape = search)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) + scale_y_log10(labels = comma) +
  scale_colour_manual(values = col_search) + log_x +
  labs(title = "B. Local boundary support Dα,win", subtitle = "Decoys with q in (α, 1.2α]",
       x = "α", y = "Dα,win") + lt_pass + th + theme(legend.position = c(0.75, 0.3))

el <- fix(series(function(x) { e <- dfdr_elasticity(x, ALPHA_GRID, eps = EPS); data.frame(alpha = e$alpha, S = e$jaccard) }))
pC4 <- ggplot(el, aes(alpha, S, colour = search, linetype = pass, shape = search)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = col_search) + log_x +
  labs(title = "C. Threshold elasticity Sα(0.2)", subtitle = "List overlap: A(α) vs A(1.2α)",
       x = "α", y = "Jaccard overlap") + lt_pass + th + theme(legend.position = c(0.75, 0.25))

## D. decoy-reference PEP concordance (second pass; first pass semi-transparent)
relfun <- function(x) { r <- dfdr_pep_reliability_tdc(x, breaks = seq(0, 1, by = 0.05)); data.frame(pep = r$pep_mean_targets, ehat = r$err_hat_target, n = r$n_target + r$n_decoy) }
pr <- cbind(list = "global list", fix(series(relfun)))
prx <- rbind(cbind(search = "Basic", pass = "Second pass", relfun(U$basic$runxprec)),
             cbind(search = "Semi-tryptic", pass = "Second pass", relfun(U$semi$runxprec)))
pr <- rbind(pr, cbind(list = "run \u00d7 precursor", fix(prx)))
pr <- pr[pr$n >= 200, ]   # the IPE itself uses bins with mean PEP <= 0.5; the panel shows the whole range
pr$n_capped <- pmin(pr$n, 2000)      # point area: bin size, capped at 2,000 entries
pD4 <- ggplot(pr, aes(pep, ehat, colour = search, shape = list, size = n_capped, alpha = pass)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey30") + geom_point() +
  scale_alpha_manual(values = c("First pass" = 0.35, "Second pass" = 0.9)) +
  scale_shape_manual(values = c(16, 4)) +
  scale_colour_manual(values = col_search) + scale_size_area(max_size = 8, name = "Bin size\n(capped at 2,000)", breaks = c(500, 1000, 2000)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "D. Decoy-reference PEP concordance", x = "Mean predicted PEP", y = "D_b / T_b (decoy-derived reference)") +
  th + theme(legend.position = "right", legend.text = element_text(size = 8))

## E. equal-chance check: decoy fraction by q band
eq <- fix(series(function(x) { b <- dfdr_equal_chance_qbands(x)$bands; data.frame(q = b$q_mean, frac = b$decoy_frac, n = b$n) }))
pE4 <- ggplot(eq, aes(q, frac, colour = search, linetype = pass, shape = search)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey30") +
  geom_line(linewidth = 0.9) + geom_point(aes(size = n), show.legend = FALSE) +
  scale_colour_manual(values = col_search) + scale_size_area(max_size = 7) +
  labs(title = "E. Equal-chance check", subtitle = "Decoy fraction by q-value band (global lists)",
       x = "Mean q-value in band", y = "Decoy fraction") + lt_pass + th + theme(legend.position = c(0.75, 0.3))

sp <- fix(series(function(x) { s <- dfdr_sumpep(x, ALPHA_GRID); data.frame(alpha = s$alpha, sumpep = s$sum_PEP) }))
pF4 <- ggplot(sp, aes(alpha, sumpep, colour = search, linetype = pass, shape = search)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) + scale_y_log10(labels = comma) +
  scale_colour_manual(values = col_search) + log_x +
  labs(title = "F. Expected false targets", x = "α", y = "Σ PEP (accepted targets)") + lt_pass + th +
  theme(legend.position = c(0.75, 0.3))

fig4 <- (pA4 | pB4 | pC4) / (pD4 | pE4 | pF4)
save_figure(fig4, "Figure4", 16, 8.5)

## =====================================================================
##  Figure 5 --- MS2Rescore case study
## =====================================================================
message("Figure 5")
TB <- readRDS(file.path(OUT, "msrescore_universes.rds"))
ALPHA5  <- ALPHA_GRID[ALPHA_GRID >= 5e-4]              # below 5e-4 no list holds a decoy
ALPHA5C <- ALPHA5[ALPHA5 <= 0.1]                        # D_win: above 10% the window (1.2 alpha) leaves the export
nm_map <- c("before rescoring (Andromeda, TDC)" = "Before (q from before score)",
            "mokapot q-value (exported, all candidates)" = "After (mokapot q-value)",
            "after rescoring (mokapot score, TDC)" = "After (q from mokapot score)")
col_cons <- c("After (mokapot q-value)" = "#F8766D", "After (q from mokapot score)" = "#00BA38", "Before (q from before score)" = "#619CFF")
cons_series <- function(fun) {
  d <- do.call(rbind, lapply(names(TB), function(k) cbind(construction = nm_map[[k]], fun(TB[[k]]))))
  d$construction <- factor(d$construction, names(col_cons)); d }

## A. score densities before and after
dens <- function(x, lab) { x <- x[is.finite(x$score), ]; data.frame(stage = lab, score = x$score, class = ifelse(x$is_decoy, "Decoy", "Target")) }
sd <- rbind(dens(TB[["before rescoring (Andromeda, TDC)"]], "before"), dens(TB[["after rescoring (mokapot score, TDC)"]], "after"))
sd$stage <- factor(sd$stage, c("before", "after"))
## densities on the full score range, then zoomed to [0, 300] (before) and [-4, 4] (after)
zoom <- list(before = c(0, 300), after = c(-4, 4))
dd <- do.call(rbind, lapply(split(sd, list(sd$stage, sd$class)), function(g) {
  d <- stats::density(g$score, n = 1024); z <- zoom[[as.character(g$stage[1])]]
  keep <- d$x >= z[1] & d$x <= z[2]
  data.frame(stage = g$stage[1], class = g$class[1], score = d$x[keep], dens = d$y[keep]) }))
blank <- data.frame(stage = factor(c("before", "before", "after", "after"), c("before", "after")),
                    score = c(0, 300, -4, 4), dens = 0, class = "Target")
pA5 <- ggplot(dd, aes(score, dens, colour = class)) + geom_line(linewidth = 0.9) +
  geom_blank(data = blank) +
  facet_wrap(~stage, ncol = 1, scales = "free") +
  scale_colour_manual(values = c(Decoy = "#F8766D", Target = "#00BFC4")) +
  labs(title = "A. Score separation (targets vs decoys)", x = "Score", y = "Density") + th +
  theme(legend.position = c(0.85, 0.3), strip.text = element_text(face = "italic", size = 11))

st5 <- cons_series(function(x) { s <- dfdr_curve_stability(x, ALPHA5); data.frame(alpha = s$alpha, T = s$T_alpha, si = s$support_index) })
log_x5 <- scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 2e-1), labels = c("0.001", "0.01", "0.10", "0.20"))
pB5 <- ggplot(st5, aes(alpha, si, colour = construction)) + geom_line(linewidth = 1) + log_x5 +
  scale_colour_manual(values = col_cons) +
  labs(title = "B. Decoy-support index Dα^-1/2 vs α", x = "α (log scale)", y = "Dα^-1/2") + th +
  theme(legend.position = c(0.65, 0.85))
lt5 <- cons_series(function(x) { w <- dfdr_local_tail(x, ALPHA5C, win_rel = EPS, truncation = "drop"); data.frame(alpha = w$alpha, Dwin = w$D_alpha_win) })
pC5 <- ggplot(lt5, aes(alpha, Dwin, colour = construction)) + geom_line(linewidth = 1) + log_x5 +
  scale_colour_manual(values = col_cons) +
  labs(title = "C. Local decoy support Dα,win vs α", x = "α (log scale)", y = "Dα,win") + th +
  theme(legend.position = c(0.35, 0.85))
pD5 <- ggplot(st5, aes(alpha, T, colour = construction)) + geom_line(linewidth = 1) + log_x5 +
  scale_colour_manual(values = col_cons) + scale_y_continuous(labels = comma) +
  labs(title = "D. Accepted targets vs α", x = "α (log scale)", y = "# accepted targets") + th +
  theme(legend.position = c(0.65, 0.2))

## E. overlap of accepted spectra before vs after (common run x scan key, 06b)
ov <- read.csv(file.path(OUT, "tables", "msrescore_overlap_scankey.csv"))
ovl <- ov[ov$alpha %in% ALPHA5, ]
ov_long <- rbind(
  data.frame(alpha = ovl$alpha, pair = "Before vs After (mokapot q-value)", J = ovl$J_before_mokapot),
  data.frame(alpha = ovl$alpha, pair = "Before vs After (q from mokapot score)", J = ovl$J_before_rescoredTDC))
ov_long <- ov_long[is.finite(ov_long$J), ]
pE5 <- ggplot(ov_long, aes(alpha, J, colour = pair)) + geom_line(linewidth = 1) + log_x5 +
  scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
  labs(title = "E. Accepted target-set overlap before vs after", x = "α (log scale)", y = "Jaccard") + th +
  theme(legend.position = c(0.6, 0.2))

eq5 <- cons_series(function(x) { b <- dfdr_equal_chance_qbands(x)$bands; data.frame(q = b$q_mean, frac = b$decoy_frac, n = b$n) })
pF5 <- ggplot(eq5, aes(q, frac, colour = construction)) + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_line(linewidth = 0.9) + geom_point(aes(size = n)) + scale_size_area(max_size = 8, name = "Bin size") +
  scale_colour_manual(values = col_cons) +
  labs(title = "F. Decoy fraction by q-band", x = "Mean q in band", y = "Decoy fraction") + th +
  theme(legend.position = "right", legend.text = element_text(size = 8))

fig5 <- (pA5 | pB5 | pC5) / (pD5 | pE5 | pF5)
save_figure(fig5, "Figure5", 16, 9)

## --- record the values behind every panel ------------------------------
write.csv(A, file.path(FIGOUT, "Figure3A_values.csv"), row.names = FALSE)
write.csv(B, file.path(FIGOUT, "Figure3B_values.csv"), row.names = FALSE)
write.csv(C, file.path(FIGOUT, "Figure3C_values.csv"), row.names = FALSE)
write.csv(D, file.path(FIGOUT, "Figure3D_values.csv"), row.names = FALSE)
write.csv(st, file.path(FIGOUT, "Figure4A_values.csv"), row.names = FALSE)
write.csv(lt, file.path(FIGOUT, "Figure4B_values.csv"), row.names = FALSE)
write.csv(el, file.path(FIGOUT, "Figure4C_values.csv"), row.names = FALSE)
write.csv(pr, file.path(FIGOUT, "Figure4D_values.csv"), row.names = FALSE)
write.csv(eq, file.path(FIGOUT, "Figure4E_values.csv"), row.names = FALSE)
write.csv(sp, file.path(FIGOUT, "Figure4F_values.csv"), row.names = FALSE)
write.csv(st5, file.path(FIGOUT, "Figure5BD_values.csv"), row.names = FALSE)
write.csv(lt5, file.path(FIGOUT, "Figure5C_values.csv"), row.names = FALSE)
write.csv(ov_long, file.path(FIGOUT, "Figure5E_values.csv"), row.names = FALSE)
write.csv(eq5, file.path(FIGOUT, "Figure5F_values.csv"), row.names = FALSE)
writeLines(c(paste("diagFDR", utils::packageVersion("diagFDR")), paste("R", R.version.string), format(Sys.time())),
           file.path(FIGOUT, "figures_provenance.txt"))
message("done")
