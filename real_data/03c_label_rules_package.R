## =====================================================================
##  03c --- Table 3 Panel B recomputed with diagFDR >= 0.2.0 only
##
##  Reproduces the label-rule audit (Note S9) and the picked-protein row
##  from the package functions, and ASSERTS the numbers of the paper so
##  that any divergence between package and scripts is caught.
##
##  Expected (basic search, Global.PG.Q.Value <= 1%):
##    split 8,991/90  winning_row 8,991/64  majority 8,987/68
##    any_decoy 8,965/90   picked 9,146/91   mixed groups at 1%: 26
##  Verified 2026-10-25 with diagFDR 0.2.0 on the deposited exports.
##  Output: out/tables/label_rules_package.csv, picked_protein_package.csv
## =====================================================================
source("00_config.R")
load_numbers()

rep <- sanitise_report(diagFDR::read_diann_parquet(DIANN_BASIC), COL$decoy, "basic")

say("label intrinsic at precursor level? %s",
    diagFDR::label_is_intrinsic(rep, COL$precursor)$intrinsic)
pg_mixed <- diagFDR::label_is_intrinsic(rep, COL$protein)
say("protein groups: %d identifiers, %d carrying both labels", pg_mixed$n_keys, pg_mixed$n_mixed)
reg("pgMixedGroups", pg_mixed$n_mixed, "int", "Sec 4.3: groups with both labels (full export)")

## Panel B is computed on the export PRE-FILTERED at Global.PG.Q.Value <= 1%,
## the export a laboratory would normally keep (Note S9). The unfiltered export
## is audited next: there decoy-if-any destroys the list (T = 963, D = 8,092).
sub <- rep[rep[[COL$q_pg_glob]] <= ALPHA_MAIN, , drop = FALSE]
aud <- diagFDR::dfdr_label_rule_audit(sub, COL$protein, COL$q_pg_glob, alpha = ALPHA_MAIN,
                                      q_max_export = ALPHA_MAIN, two_pass = TRUE)
print(aud); write_tab(aud, "label_rules_package")
aud_full <- diagFDR::dfdr_label_rule_audit(rep, COL$protein, COL$q_pg_glob, alpha = ALPHA_MAIN,
                                           q_max_export = QMAX$q_pg_glob, two_pass = TRUE)
print(aud_full); write_tab(aud_full, "label_rules_package_unfiltered")
stopifnot(aud_full$T_alpha[aud_full$rule == "any_decoy"] == 963,
          aud_full$D_alpha[aud_full$rule == "any_decoy"] == 8092,
          aud_full$D_alpha[aud_full$rule == "majority"] == 236)
reg("pgAnyDecoyFullT", 963, "int", "Sec 4.2: decoy-if-any on the unfiltered export")
reg("pgAnyDecoyFullD", 8092, "int", "Sec 4.2: decoy-if-any on the unfiltered export")

expected <- data.frame(rule = c("split", "winning_row", "majority", "any_decoy"),
                       T_alpha = c(8991, 8991, 8987, 8965), D_alpha = c(90, 64, 68, 90))
m <- merge(aud, expected, by = "rule", suffixes = c("", ".exp"))
stopifnot(all(m$T_alpha == m$T_alpha.exp), all(m$D_alpha == m$D_alpha.exp))
say("Panel B label rules reproduced by the package.")

for (r in aud$rule) {
  nm <- gsub("_", "", gsub("(^|_)([a-z])", "\\U\\2", r, perl = TRUE)); reg(paste0("pgT", nm), aud$T_alpha[aud$rule == r], "int", paste("Panel B", r))
  reg(paste0("pgD", nm), aud$D_alpha[aud$rule == r], "int", paste("Panel B", r))
  reg(paste0("pgRatio", nm), aud$decoy_ratio[aud$rule == r], "d3", paste("Panel B", r))
}

## picked-protein competition on the split universe (full export, not the 1% filter)
pg_split <- diagFDR::diann_global_proteingroup(rep, q_col = COL$q_pg_glob, rule = "split",
                                               q_max_export = QMAX$q_pg_glob, two_pass = TRUE)
pk <- diagFDR::dfdr_picked_protein(pg_split, alpha = ALPHA_MAIN, by = "q", add_decoy = 0)
say("picked: T=%d D=%d ratio=%.3f%% pairs=%d", pk$T_alpha, pk$D_alpha, 100 * pk$decoy_ratio, pk$n_pairs)
stopifnot(pk$T_alpha == 9146, pk$D_alpha == 91)
write_tab(data.frame(T_alpha = pk$T_alpha, D_alpha = pk$D_alpha, decoy_ratio = pk$decoy_ratio,
                     CV_benchmark = pk$CV_benchmark, n_pairs = pk$n_pairs), "picked_protein_package")
reg("pickedT", pk$T_alpha, "int", "Panel B picked"); reg("pickedD", pk$D_alpha, "int", "Panel B picked")
