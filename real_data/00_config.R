## =====================================================================
##  diagFDR paper --- CONFIGURATION
##
##  This is the only file you normally edit.
##
##  DATA is the folder holding the files of the Zenodo record
##  https://doi.org/10.5281/zenodo.18621806 . Set it below, or through
##  the environment variable DIAGFDR_DATA. Results are written to ./out.
##  Column names and thresholds were checked against report_basic.parquet
##  (71 columns, 1,020,465 rows) and the mokapot exports.
## =====================================================================

## ---------------------------------------------------------------------
##  INPUT --- paths
## ---------------------------------------------------------------------

## Sys.setenv(DIAGFDR_DATA = "/path/to/zenodo/files")   # or export DIAGFDR_DATA
DATA <- Sys.getenv("DIAGFDR_DATA", unset = "data")

## --- DIA-NN ------------------------------------------------------------
DIANN_BASIC  <- file.path(DATA, "diann/Basic_search/report_basic.parquet")
DIANN_SEMI   <- file.path(DATA, "diann/Semi_Tryp_search/report_semi.parquet")
DIANN_ENTRAP <- file.path(DATA, "diann/Basic_search_wih_ArabidopsisT/report_entrapment.parquet")

## First-pass reports. Figure 3 compares first-pass with reanalysis mode,
## so both are needed for that panel. Set to NA to skip the comparison.
DIANN_BASIC_FP  <- file.path(DATA, "diann/Basic_search/report_basic-first-pass.parquet")
DIANN_SEMI_FP   <- file.path(DATA, "diann/Semi_Tryp_search/report_semi-first-pass.parquet")

## --- FASTA -------------------------------------------------------------
FASTA_TARGET <- file.path(DATA, "ProteoBenchFASTA_MixedSpecies_HYE.fasta")
FASTA_CONTAM <- file.path(DATA, "uniprot-proteome-human-contaminants.fasta")
FASTA_ENTRAP <- file.path(DATA, "uniprotkb_taxonomy_id_3702_AND_reviewed_2026_02_06.fasta")

## --- MS2Rescore --------------------------------------------------------
MOKAPOT_TARGETS <- file.path(DATA, "ms2rescore/maxquant/pyro_mq.ms2rescore.mokapot.psms.txt")
MOKAPOT_DECOYS  <- file.path(DATA, "ms2rescore/maxquant/pyro_mq.ms2rescore.mokapot.decoy.psms.txt")
MAXQUANT_MSMS   <- file.path(DATA, "ms2rescore/maxquant/msms.txt")

## --- Outputs -----------------------------------------------------------
OUT     <- Sys.getenv("DIAGFDR_OUT", unset = "out")
FIGDIR  <- file.path(OUT, "figures")
TABDIR  <- file.path(OUT, "tables")
NUMFILE <- file.path(OUT, "numbers.tex")

## ---------------------------------------------------------------------
##  INPUT --- analysis parameters
## ---------------------------------------------------------------------

ALPHA_MAIN <- 0.01
ALPHA_GRID <- c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1)
RHO        <- 0.20     # win_rel for D_{alpha,win}
EPSILON    <- 0.20     # eps for the elasticity S_alpha(eps)
NMIN_BIN   <- 200      # minimum bin occupancy entering the IPE
PEP_MAX    <- 0.5      # IPE uses bins with mean PEP <= this
PEP_BINWIDTH <- 0.05
EQ_BAND    <- c(0.2, 0.5)
EQ_BREAKS  <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
SEED       <- 20260811

## DIA-NN column names. VERIFIED against report_basic.parquet.
## Note there is NO CScore / score column in this export, so `score` stays
## NA throughout.
COL <- list(
  run       = "Run",
  precursor = "Precursor.Id",
  protein   = "Protein.Group",
  protein_ids = "Protein.Ids",
  q_run     = "Q.Value",
  q_global  = "Global.Q.Value",
  q_lib     = "Lib.Q.Value",
  q_pg      = "PG.Q.Value",
  q_pg_glob = "Global.PG.Q.Value",
  pep       = "PEP",
  decoy     = "Decoy"
)

## Export ceilings. THESE DIFFER BY COLUMN and dfdr_local_tail() /
## dfdr_equal_chance_qbands() both use them: a window or band above the
## ceiling is dropped. Verified: max(Q.Value) = 0.500 (the --qvalue 0.5
## export), max(Global.Q.Value) = 1.000.
QMAX <- list(q_run = 0.5, q_global = 1.0, q_pg = 1.0, q_pg_glob = 1.0)

## Entrapment attribution.
##
## IMPORTANT: Protein.Group holds BARE UniProt accessions (P55011, ...),
## with no species tag, so a string prefix such as "ARATH" matches
## nothing. Attribution is by membership of the accession in the
## A. thaliana FASTA. "pure" counts an identification towards E_alpha only
## when EVERY accession in its Protein.Group is an entrapment accession.
ENTRAP_RULE   <- "pure"
ENTRAP_COL    <- "protein"   # "protein" = Protein.Group (used in the paper)
                             # "protein_ids" = Protein.Ids (much more permissive)

## Protein-level label propagation.
## At protein-group level a group routinely contains BOTH target and decoy
## precursors (10,262 of 13,751 in the basic search), so decoy-if-any
## labels nearly every group a decoy and the ratio explodes. "split"
## treats target and decoy entries as distinct hypotheses, which is
## DIA-NN's own convention and reproduces the paper's numbers.
PROTEIN_LABEL_RULE <- "split"

## ---------------------------------------------------------------------
##  END OF INPUT
## ---------------------------------------------------------------------

set.seed(SEED)
for (d in c(OUT, FIGDIR, TABDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

CONFIG_DIR <- tryCatch({
  f <- sys.frame(1)$ofile
  if (is.null(f)) "." else dirname(f)
}, error = function(e) ".")
if (is.na(CONFIG_DIR) || !nzchar(CONFIG_DIR)) CONFIG_DIR <- "."

source(file.path(CONFIG_DIR, "R", "dfdr_api.R"))
source(file.path(CONFIG_DIR, "R", "utils.R"))

message("diagFDR ", DFDR_VERSION, " | output -> ", OUT)
