library(tidyverse)
source("~/Documents/GitHub/clk1-splicing/figures/theme_for_plots.R")
low_clk <- read_lines("~/Downloads/low_clk.txt")
low_bs <- v15 %>%
  filter(Kids_First_Biospecimen_ID %in% low_clk)



splicing_genes <- read_lines("~/Documents/GitHub/clk1-splicing/analyses/splicing-factor_dysregulation/input/splicing_factors.txt")

func_sites <- read_tsv("~/Documents/GitHub/clk1-splicing/analyses/splicing_events_functional_sites/results/kinases-functional_sites.tsv") %>%
  filter(Differential_freq >10,
         #mean_dPSI >= 0.2,
         Average_TPM >= 10) %>%
  mutate(splice_gene = ifelse(gene %in% splicing_genes, "yes", "no"))

func_sites_filt <- func_sites %>%
#  filter(splice_gene == "yes") %>%
  select(gene, mean_dPSI, Average_TPM, splice_gene) %>%
  group_by(gene, Average_TPM, splice_gene) %>%
  summarise(max_dPSI = max(mean_dPSI)) %>%
  unique() 

dep <- read_csv("~/Downloads/CRISPRGeneEffect.csv")
models <- read_csv("~/Downloads/Model.csv") %>%
  filter(Age <= 40,
         OncotreeLineage == "CNS/Brain",
         grepl("Glioma", OncotreePrimaryDisease))


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
  write_csv("~/Downloads/mean_ped_glioma_crispr_scores_func_kinase_splice_events.csv")
 # filter(!is.na(mean_CRISPR_score))

morph <- read_tsv("~/Documents/GitHub/clk1-splicing/data/ctrl_vs_morpho.rsem.genes.results.tsv") %>%
  mutate(gene_symbol = sub("^[^_]+_", "", gene)) %>%
  rowwise() %>% # Operate row-wise for calculations
  mutate(
    Mean_CTRL = mean(c_across(starts_with("CTRL"))),
    Mean_Treated = mean(c_across(starts_with("Treated"))),
    Fold_Change = Mean_Treated / Mean_CTRL # Fold change calculation
  ) %>%
  ungroup() # Remove row-wise operation

de <- read_tsv("~/Documents/GitHub/clk1-splicing/analyses/CLK1-splicing-impact-morpholino/results/ctrl_vs_treated.de.tsv")
clk1_targets <- read_lines("~/Downloads/tam_etal_clk1_targets.txt")
de_clk <- de %>%
  filter(Gene_Symbol %in% c(clk1_targets, "CLK1")) %>%
  select(Gene_Symbol, baseMean, log2FoldChange, padj, Significant) %>%
  arrange(log2FoldChange) %>%
  write_tsv("~/Downloads/clk1_consensus_targets.tsv")

