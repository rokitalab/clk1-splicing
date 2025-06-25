# Author: Jo Lynne Rokita
# Function: Create geneset to use for pathway analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
data_dir <- file.path(root_dir, "data")

results_dir <- file.path(analysis_dir, "results")
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(analysis_dir, "input")

# output file
hallmark_kegg_splice_rds_output <- file.path(results_dir,"hallmark_kegg_splice_geneset_mrna.rds")

# read in genesets to keep
genesets_to_keep <- read_tsv(file.path(input_dir, "genesets.tsv")) %>%
  filter(Keep == "Yes")

# create pathways input file - we will use Hallmark + KEGG splice
hallmark_db  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique()

kegg_splice <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::filter(gs_name == "KEGG_SPLICEOSOME") %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique()

hallmark_kegg_splice <- hallmark_db %>%
  rbind(kegg_splice) %>%
  filter(gs_name %in% genesets_to_keep$Geneset) %>%
  unstack()

saveRDS(hallmark_kegg_splice, file = hallmark_kegg_splice_rds_output)
