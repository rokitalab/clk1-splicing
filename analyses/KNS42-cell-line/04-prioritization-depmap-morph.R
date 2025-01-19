################################################################################
# 04-prioritization-depmap-morph.R
# script that explores depmap scores in context of genes 
# with functional splicing, clk1 targets
# written by Jo Lynne Rokita
#
# usage: Rscript 04-prioritization-depmap-morph.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("tidyverse")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses")
kns_dir <- file.path(analysis_dir, "KNS42-cell-line")
sf_dir <- file.path(analysis_dir, "splicing-factor_dysregulation", "input")
func_dir <- file.path(analysis_dir, "splicing_events_functional_sites", "results")
morph_dir <- file.path(analysis_dir, "CLK1-splicing-impact-morpholino", "results")
input_dir <- file.path(kns_dir, "input")
plots_dir <- file.path(kns_dir, "plots")
res_dir <- file.path(kns_dir, "results")


if(!dir.exists(res_dir)){
  dir.create(res_dir, recursive=TRUE)
}

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir,"theme_for_plots.R"))

# read in splicing factor related genes
splicing_genes <- read_lines(file.path(sf_dir, "splicing_factors.txt"))

# read in genes with spliced exons annotated to functional sites
func_sites <- read_tsv(file.path(func_dir, "kinases-functional_sites.tsv")) %>%
  mutate(splice_gene = ifelse(gene %in% splicing_genes, "yes", "no"))

# which genes are involved in splicing?
func_sites %>%
  filter(splice_gene == "yes") %>%
  pull(gene) %>%
  unique()

#[1] "CLK1"  "CLK2"  "CLK3"  "CLK4"  "SRPK1"

func_sites_filt <- func_sites %>%
  select(gene, mean_dPSI, Average_TPM, splice_gene) %>%
  group_by(gene, Average_TPM, splice_gene) %>%
  summarise(max_dPSI = max(mean_dPSI)) %>%
  unique() 

func_sites_filt_splice %>% func_sites_filt
  filter(splice_gene == "yes") %>%
  write_tsv(file.path(res_dir, "splice_genes_functional.tsv"))
  
# explore whether any of these also overlap as dependencies
dep <- read_csv(file.path(input_dir, "CRISPRGeneEffect.csv"))
models <- read_csv(file.path(input_dir, "Model.csv")) %>%
# filter for pediatric, AYA
    filter(Age < 40,
           # filter brain only
         OncotreeLineage == "CNS/Brain",
         # glioma only
         grepl("Glioma", OncotreePrimaryDisease),
         # rm sarcoma
         OncotreeCOde != "GSARC")

crispr_df_long <- dep %>%
  filter(...1 %in% models$ModelID) %>%  # Filter for a specific row
  dplyr::rename(model_id=...1) %>%            # Rename the first column to "gene"
  pivot_longer(cols = -model_id,         # Convert all columns except "gene" to rows
               names_to = "gene",
               values_to = "CRISPR_score") %>%
  mutate(gene = str_remove_all(gene, "\\s*\\(.*\\)")) # Remove spaces and content within parentheses

mean_scores <- crispr_df_long %>%
  group_by(gene) %>%
  summarise(mean_CRISPR_score = mean(CRISPR_score, na.rm = TRUE)) %>%
  filter(mean_CRISPR_score != "NaN",
         mean_CRISPR_score < -0.1) %>%
  filter(gene %in% func_sites$gene) %>%
  right_join(func_sites_filt) %>%
  arrange(mean_CRISPR_score) %>%
  write_csv(file.path(res_dir, "mean_ped_glioma_crispr_scores_func_kinase_splice_events.csv"))

# investigate CLK1 targets from morph expt  
morph <- read_tsv(file.path(data_dir, "ctrl_vs_morpho.rsem.genes.results.tsv")) %>%
  mutate(gene_symbol = sub("^[^_]+_", "", gene)) %>%
  rowwise() %>% # Operate row-wise for calculations
  mutate(
    Mean_CTRL = mean(c_across(starts_with("CTRL"))),
    Mean_Treated = mean(c_across(starts_with("Treated"))),
    Fold_Change = Mean_Treated / Mean_CTRL # Fold change calculation
  ) %>%
  ungroup() # Remove row-wise operation

# pull DE genes post morph tx
de <- read_tsv(file.path(morph_dir, "ctrl_vs_treated.de.tsv"))
# intersect with known CLK1 targets from Tam, et al.
# Tam, B.Y., Chiu, K., Chung, H., Bossard, C., Nguyen, J.D., Creger, E., Eastman, B.W., Mak, C.C., Ibanez, M., Ghias, A., et al. (2020). 
# The CLK inhibitor SM08502 induces anti-tumor activity and reduces Wnt pathway gene expression in gastrointestinal cancer models. 
# Cancer Lett. 473, 186–197.
clk1_targets <- read_lines(file.path(input_dir, "tam_etal_clk1_targets.txt"))

de_clk <- de %>%
  dplyr::filter(Gene_Symbol %in% c(clk1_targets, "CLK1")) %>%
  select(Gene_Symbol, baseMean, log2FoldChange, padj, Significant) %>%
  arrange(log2FoldChange) %>%
  write_tsv(file.path(res_dir, "clk1_consensus_targets.tsv"))

