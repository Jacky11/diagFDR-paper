## =====================================================================
##  03 --- Scope audit: Table 2 Panel A and Panel B, Figure 4A
##
##  Produces: Table 2 and the decoy ratios of Section 4.2.
##
##  Expected Panel A values on the Zenodo data, at alpha = 1%:
##    Basic  Global.Q.Value        T = 83,149  D =   832  1.001%  CV 0.0347
##    Basic  min_r(Q.Value)        T = 89,718  D = 2,001  2.230%  CV 0.0224
##    Basic  Q.Value (run x prec)  T =375,129  D = 3,751  1.000%  CV 0.0163
##  If your run differs from these, something upstream changed --- stop and
##  find out what before using the results.
## =====================================================================

source("00_config.R")
load_numbers()
U <- readRDS(file.path(OUT, "universes.rds"))

UNIV_LABELS <- c(global    = "Global.Q.Value",
                 minrunq   = "Q.Value aggregated by min_r",
                 runxprec  = "Q.Value (run x precursor)",
                 pg_global = "Global.PG.Q.Value",
                 pg_run    = "PG.Q.Value")

## ---------------------------------------------------------------------
##  Panel A
## ---------------------------------------------------------------------
rowsA <- do.call(rbind, lapply(c("basic", "semi"), function(s)
  do.call(rbind, lapply(names(UNIV_LABELS), function(u) {
    if (is.null(U[[s]][[u]])) return(NULL)
    h <- headline_of(U[[s]][[u]])
    data.frame(search = if (s == "basic") "Basic" else "Semi",
               universe = UNIV_LABELS[[u]], level = if (grepl("^pg", u)) "protein" else "precursor",
               T_alpha = h$T_alpha, D_alpha = h$D_alpha,
               decoy_ratio = h$decoy_ratio, CV = h$CV,
               stringsAsFactors = FALSE)
  }))))
print(rowsA); write_tab(rowsA, "table2_panelA")

g <- function(sr, un, col) {
  v <- rowsA[rowsA$search == sr & rowsA$universe == un, col]
  if (!length(v)) NA else v[1]
}
reg("TalphaBasic",      g("Basic","Global.Q.Value","T_alpha"), "int", "Table 2; the Discussion's verifiable-claim example")
reg("DalphaBasic",      g("Basic","Global.Q.Value","D_alpha"), "int", "ditto")
reg("CVBasic",          g("Basic","Global.Q.Value","CV"), "d3", "ditto")
reg("ratioBasicGlobal", g("Basic","Global.Q.Value","decoy_ratio"), "pct2", "Sec 4.2")
reg("ratioBasicMinrun", g("Basic","Q.Value aggregated by min_r","decoy_ratio"), "pct2", "Sec 4.2")
reg("ratioSemiMinrun",  g("Semi","Q.Value aggregated by min_r","decoy_ratio"), "pct2", "Sec 4.2")
reg("DalphaRunPrec",    g("Basic","Q.Value (run x precursor)","D_alpha"), "int", "Sec 4.3 negative control")
reg("CVRunPrec",        g("Basic","Q.Value (run x precursor)","CV"), "d3", "Sec 4.3")
if (!is.na(g("Basic","Global.PG.Q.Value","T_alpha"))) {
  reg("TpgGlobal",  g("Basic","Global.PG.Q.Value","T_alpha"), "int", "Sec 4.2 protein level")
  reg("DpgGlobal",  g("Basic","Global.PG.Q.Value","D_alpha"), "int", "Sec 4.2")
  reg("ratioPgGlobal", g("Basic","Global.PG.Q.Value","decoy_ratio"), "pct2", "Sec 4.2")
  reg("TpgRun",     g("Basic","PG.Q.Value","T_alpha"), "int", "Sec 4.2")
  reg("DpgRun",     g("Basic","PG.Q.Value","D_alpha"), "int", "Sec 4.2")
  reg("ratioPgRun", g("Basic","PG.Q.Value","decoy_ratio"), "pct2", "Sec 4.2")
}

k_runs <- attr(U$basic, "n_runs"); if (!is.finite(k_runs) || k_runs < 2) k_runs <- 6
reg("minIndepPct", 1 - (1 - ALPHA_MAIN)^k_runs, "pct2",
    sprintf("1-(1-alpha)^k with k=%d runs, quoted in Sec 4.2 and Methods", k_runs))
reg("nRuns", k_runs, "int", "runs in the DIA-NN experiment")

## ---------------------------------------------------------------------
##  Panel B --- sensitivity to the label-propagation rule
##
##  Reported at BOTH levels, because the answer differs sharply:
##  at precursor level the label is effectively intrinsic and the rules
##  coincide; at protein-group level they do not, and decoy-if-any is
##  simply invalid.
## ---------------------------------------------------------------------
li_p <- attr(U$basic, "label_prec"); li_g <- attr(U$basic, "label_prot")

rowsB <- do.call(rbind, lapply(c("basic","semi"), function(s) rbind(
  do.call(rbind, lapply(names(U[[s]]$minrunq_by_rule), function(r) {
    h <- headline_of(U[[s]]$minrunq_by_rule[[r]])
    data.frame(search = if (s=="basic") "Basic" else "Semi", level = "precursor",
               universe = "min_r(Q.Value)", label_rule = r,
               T_alpha = h$T_alpha, D_alpha = h$D_alpha,
               decoy_ratio = h$decoy_ratio, CV = h$CV, stringsAsFactors = FALSE) })),
  do.call(rbind, lapply(names(U[[s]]$pg_by_rule), function(r) {
    h <- headline_of(U[[s]]$pg_by_rule[[r]])
    data.frame(search = if (s=="basic") "Basic" else "Semi", level = "protein",
               universe = "Global.PG.Q.Value", label_rule = r,
               T_alpha = h$T_alpha, D_alpha = h$D_alpha,
               decoy_ratio = h$decoy_ratio, CV = h$CV, stringsAsFactors = FALSE) })))))
print(rowsB); write_tab(rowsB, "table2_panelB")

pb <- function(lv, r, col) {
  v <- rowsB[rowsB$search=="Basic" & rowsB$level==lv & rowsB$label_rule==r, col]
  if (!length(v)) NA else v[1]
}
reg("precMixed", li_p$n_mixed, "int", "precursors carrying both labels")
reg("precKeys",  li_p$n_keys,  "int", "precursors in total")
reg("protMixed", li_g$n_mixed, "int", "protein groups carrying both labels")
reg("protKeys",  li_g$n_keys,  "int", "protein groups in total")
reg("rulePrecAny",  pb("precursor","any_decoy","decoy_ratio"),   "pct2", "Table 2B precursor")
reg("rulePrecWin",  pb("precursor","winning_row","decoy_ratio"), "pct2", "Table 2B precursor")
reg("rulePrecMaj",  pb("precursor","majority","decoy_ratio"),    "pct2", "Table 2B precursor")
reg("ruleProtAny",  pb("protein","any_decoy","decoy_ratio"),     "pct1", "Table 2B protein -- the failure case")
reg("ruleProtSplit",pb("protein","split","decoy_ratio"),         "pct2", "Table 2B protein -- correct")

prec_spread <- diff(range(rowsB$decoy_ratio[rowsB$search=="Basic" & rowsB$level=="precursor"], na.rm = TRUE))
reg("rulePrecSpread", prec_spread, "pct2", "spread across rules at precursor level")

## ---------------------------------------------------------------------
##  Scope disagreement (Jaccard overlap), reported descriptively
## ---------------------------------------------------------------------
J <- scope_disagreement_of(U$basic$runxprec, U$basic$global, ALPHA_GRID,
                           "run-wise", "global")
print(J); write_tab(J, "scope_jaccard")
reg("JrunGlobal", J$J[J$alpha == ALPHA_MAIN], "d2",
    "Jaccard between run-wise and global accepted lists at 1%")

p <- open_pdf("fig5B_jaccard", 5.5, 4.2)
plot(J$alpha, J$J, type = "b", pch = 19, log = "x", ylim = c(0, 1),
     xlab = expression(alpha), ylab = expression(J[run/global]),
     main = "Scope disagreement (descriptive)")
abline(h = 1, lty = 3); dev.off(); message("  -> ", p)

p <- open_pdf("fig5A_decoy_ratio", 6.5, 4.2)
bp <- rowsA[rowsA$universe %in% c("Global.Q.Value","Q.Value aggregated by min_r"), ]
m <- matrix(100 * bp$decoy_ratio, nrow = 2,
            dimnames = list(c("Global.Q.Value","min_r(Q.Value)"), c("Basic","Semi")))
barplot(m, beside = TRUE, ylab = "decoy ratio on the reported list (%)",
        legend.text = rownames(m), args.legend = list(x = "topleft", bty = "n"))
abline(h = 100 * ALPHA_MAIN, lty = 2); dev.off(); message("  -> ", p)

## ---------------------------------------------------------------------
##  Table 2 as LaTeX
## ---------------------------------------------------------------------
f <- function(x, d = 3) formatC(x, format = "f", digits = d)
i <- function(x) formatC(x, format = "d", big.mark = "{,}")
tex <- c("% Table 2 --- generated by 03_scope_audit.R",
         "\\multicolumn{6}{l}{\\textit{Global precursor list}} \\\\")
for (s in c("Basic","Semi")) for (u in c("Global.Q.Value","Q.Value aggregated by min_r")) {
  r <- rowsA[rowsA$search == s & rowsA$universe == u, ]; if (!nrow(r)) next
  un <- if (grepl("min_r", u)) "\\texttt{Q.Value} aggregated by $\\min_r$" else "\\texttt{Global.Q.Value}"
  tex <- c(tex, sprintf("%s & %s & %s & %s & %s\\%% & %s \\\\", s, un,
                        i(r$T_alpha), i(r$D_alpha), f(100*r$decoy_ratio), f(r$CV, 4)))
}
tex <- c(tex, "\\midrule", "\\multicolumn{6}{l}{\\textit{Run$\\times$precursor}} \\\\")
for (s in c("Basic","Semi")) {
  r <- rowsA[rowsA$search == s & rowsA$universe == "Q.Value (run x precursor)", ]; if (!nrow(r)) next
  tex <- c(tex, sprintf("%s & \\texttt{Q.Value} & %s & %s & %s\\%% & %s \\\\", s,
                        i(r$T_alpha), i(r$D_alpha), f(100*r$decoy_ratio), f(r$CV, 4)))
}
tex <- c(tex, "\\midrule", "\\multicolumn{6}{l}{\\textit{Protein group}} \\\\")
for (s in c("Basic","Semi")) for (u in c("Global.PG.Q.Value","PG.Q.Value")) {
  r <- rowsA[rowsA$search == s & rowsA$universe == u, ]; if (!nrow(r)) next
  tex <- c(tex, sprintf("%s & \\texttt{%s} & %s & %s & %s\\%% & %s \\\\", s, u,
                        i(r$T_alpha), i(r$D_alpha), f(100*r$decoy_ratio), f(r$CV, 4)))
}
tex <- c(tex, "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Panel B: sensitivity to the label-propagation rule}} \\\\")
for (k in seq_len(nrow(rowsB))) {
  r <- rowsB[k, ]
  tex <- c(tex, sprintf("%s & %s, %s (%s) & %s & %s & %s\\%% & %s \\\\",
    r$search, r$universe, gsub("_", " ", r$label_rule), r$level,
    i(r$T_alpha), i(r$D_alpha), f(100*r$decoy_ratio), f(r$CV, 4)))
}
writeLines(tex, file.path(TABDIR, "table2.tex"))
message("  -> ", file.path(TABDIR, "table2.tex"))

save_numbers()
message("\nNext: 04_calibration_stability.R")
