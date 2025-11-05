# R Corbett 2025
#
# Create SE splice event PSI matrix

# Load libraries
library(tidyverse)
library(data.table)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")

# Set file paths
hist_file <- file.path(root_dir, "analyses",
                       "cohort_summary", "results",
                       "histologies-plot-group.tsv")

indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")
  
# wrangle data
samples <- read_tsv(indep_file) %>%
  filter(cohort == "PBTA") %>%
  pull(Kids_First_Biospecimen_ID)

# Load rmats, filter for samples of interest and SE events, and select relevant columns
rmats <- data.table::fread(rmats_file) %>%
  dplyr::filter(sample_id %in% samples) %>%
  dplyr::select(sample_id, geneSymbol, chr,
                strand, exonStart_0base, exonEnd,
                upstreamES, upstreamEE,
                downstreamES, downstreamEE,
                IncLevel1)


# create splice_id column
rmats <- rmats %>%
  dplyr::mutate(exonStart_0base = exonStart_0base + 1) %>%
  dplyr::mutate(splice_id = glue::glue("{chr}:{geneSymbol}_{exonStart_0base}-{exonEnd}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}_{strand}")) %>%
  dplyr::select(sample_id, splice_id, IncLevel1)

# generate sample PSI matrix
psi_mat <- rmats %>%
  pivot_wider(
    names_from = splice_id,
    values_from = IncLevel1
  )

# write matrix to output
saveRDS(psi_mat,
        file.path(results_dir,
                  "pbta-splice-event-psis.RDS"))

# print session info
sessionInfo()


