# Author: Komal S. Rathi, Jo Lynne Rokita
# Function: 
# 1) Create KEGG geneset

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(vroom)
})


#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
data_dir <- file.path(root_dir, "data")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(analysis_dir, "input")


# read histologies file
clin_file <- read_tsv(file.path(data_dir, "histologies-plot-group.tsv")) %>% 
  filter(experimental_strategy == "RNA-Seq")

# read in genesets to keep
genesets_to_keep <- read_tsv(file.path(input_dir, "genesets.tsv")) %>%
  filter(Keep == "Yes")

# read splice dataset
splice_mat <- readRDS(file.path(input_dir, "pan_cancer_splicing_SE.gene.rds"))

# 2) read pbta dataset, filter to samples and genes in pbta splice matrix and save
pbta_subset <- readRDS(file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")) %>%
  dplyr::select(colnames(splice_mat)) %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% rownames(splice_mat)) %>%
  column_to_rownames("gene_symbol")

counts_rds_output <- file.path(input_dir,"raw_counts_pbta_subset.rds")
saveRDS(pbta_subset, file = counts_rds_output)

# 3) create pathways input file - we will use Hallmark + KEGG splice
hallmark_db  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique()

kegg_splice <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::filter(gs_name == "KEGG_SPLICEOSOME") %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique()

hallmark_splice <- hallmark_db %>%
  rbind(kegg_splice) %>%
  filter(gs_name %in% genesets_to_keep$Geneset) %>%
  unstack()

hallmark_splice_rds_output <- file.path(input_dir,"hallmark_splice_geneset_mrna.rds")

saveRDS(hallmark_splice, file = hallmark_splice_rds_output)
