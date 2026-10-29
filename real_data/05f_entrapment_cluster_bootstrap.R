## =====================================================================
##  05f --- protein-group cluster bootstrap for the entrapment estimates
##
##  The Wilson intervals of 05_entrapment.R / 05b treat accepted precursors
##  as independent. This script resamples protein groups instead (Note S4,
##  "Clustered errors") for the shared / added precursors at alpha = 1% and
##  for the Global.Q.Value list at alpha = 1e-4.
##
##  Expected on the Zenodo data (diagFDR 0.2.1, seed 1):
##    core  0.217% [0.180, 0.255]   added 2.57% [2.17, 2.97]   ratio 11.9 [9.3, 14.9]
##    alpha = 1e-4: raw 5.4e-4, cluster bootstrap [2.7e-4, 9.0e-4]
## =====================================================================
source("00_config.R")
suppressPackageStartupMessages(library(diagFDR))
set.seed(1); B <- 2000
U <- readRDS(file.path(OUT, "universes.rds")); pm <- U$entrap$prot_map
ACC <- read_fasta_accessions(FASTA_ENTRAP)
is_ent <- function(p) vapply(strsplit(ifelse(is.na(p), "", p), ";"),
                             function(z) { z <- z[nzchar(z)]; length(z) > 0 && all(z %in% ACC) }, logical(1))
ent  <- pm$id[is_ent(pm$protein)]
prot <- setNames(pm$protein, pm$id)
boot <- function(ids) { e <- ids %in% ent; tab <- split(e, prot[ids]); k <- length(tab)
  replicate(B, mean(unlist(tab[sample(k, k, replace = TRUE)]))) }
ci <- function(v) stats::quantile(v, c(0.025, 0.975))

g <- dfdr_accepted(U$entrap$global, ALPHA_MAIN); m <- dfdr_accepted(U$entrap$minrunq, ALPHA_MAIN)
core <- intersect(g, m); added <- setdiff(m, g)
bc <- boot(core); ba <- boot(added)
res <- data.frame(set = c("core", "added", "ratio added/core", "global at 1e-4"),
                  estimate = c(mean(core %in% ent), mean(added %in% ent),
                               mean(added %in% ent) / mean(core %in% ent), NA),
                  lo = c(ci(bc)[1], ci(ba)[1], ci(ba / bc)[1], NA),
                  hi = c(ci(bc)[2], ci(ba)[2], ci(ba / bc)[2], NA),
                  n = c(length(core), length(added), NA, NA),
                  clusters = c(length(unique(prot[core])), length(unique(prot[added])), NA, NA))
g4 <- dfdr_accepted(U$entrap$global, 1e-4); b4 <- boot(g4)
res[4, c("estimate", "lo", "hi", "n", "clusters")] <- c(mean(g4 %in% ent), ci(b4), length(g4), length(unique(prot[g4])))
print(res); write.csv(res, file.path(TABDIR, "entrapment_cluster_bootstrap.csv"), row.names = FALSE)
