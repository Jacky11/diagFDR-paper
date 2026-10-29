## Runs the real-data analyses in order. See README.md.
scripts <- c("01_preflight.R", "02_import_diann.R", "03_scope_audit.R",
             "03b_picked_protein.R", "03c_label_rules_package.R", "04_calibration_stability.R",
             "05_entrapment.R", "05b_entrapment_discordant.R",
             "05c_entrapment_r_sensitivity.R", "05d_discordant_size_matched.R", "05e_entrapment_package.R", "05f_entrapment_cluster_bootstrap.R",
             "06_msrescore.R", "06b_msrescore_overlap.R", "06c_mokapot_universe.R", "08_defaults_sensitivity.R", "09_double_filter.R")
for (s in scripts) {
  message("\n=== ", s, " ===")
  sys.source(s, envir = new.env())
}
