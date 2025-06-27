################################################################################
# 12-plot-clk1-ex4-transcripts-by-age.R
# Script that plots CLK1 exon 4 PSI variations across tumors by age
# written by Ammar Naqvi
#
# usage: Rscript 12-plot-clk1-ex4-transcripts-by-age.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(coin)  # For oneway_test
  library(vroom)
  library(data.table)
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(root_dir, "results")

# Specify file paths
clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

# Define file paths
histology_gtex_file <- file.path(input_dir,"gtex-samples-by-age.tsv")
gtex_trans_file <- file.path(data_dir,"gtex-harmonized-isoform-expression-rsem-tpm.rds")
expr_tpm_tumor_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

# Step 1: Load and filter clinical and independent data
indep_df <- read_tsv(indep_file)
cluster_df <- read_tsv(cluster_file) %>%
  rename(Kids_First_Biospecimen_ID = sample_id) %>%
  select(Kids_First_Biospecimen_ID, cluster)

hist_indep_rna_df <- read_tsv(clin_file) %>%
  filter(
    cohort == "PBTA",
    RNA_library == "stranded",
    Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID
  )

# Step 2: Process expression data
all_clk4_transcr_counts <- readRDS(expr_tpm_tumor_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("BS"), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = -transcript_id, names_to = "Kids_First_Biospecimen_ID", values_to = "TPM") %>%
  inner_join(hist_indep_rna_df, by = "Kids_First_Biospecimen_ID") %>%
  inner_join(cluster_df, by = "Kids_First_Biospecimen_ID") %>%
  mutate(
    age_years = age_at_diagnosis_days / 365,
    age_bin = cut(
      age_years,
      breaks = c(0, 7, 14, 18, 39),
      labels = c("[0-7]", "[8-14]", "[15-18]", "[19-39]"),
      right = TRUE
    )
  ) %>%
  select(Kids_First_Biospecimen_ID, cluster, plot_group, age_bin, transcript_id, TPM) %>%
  na.omit()


# Step 3: Calculate exon proportions
transcript_expr_CLK1_combined_df <- all_clk4_transcr_counts %>%
  group_by(Kids_First_Biospecimen_ID) %>%
  mutate(
    total_TPM = sum(TPM[transcript_id %in% c("Exon 4", "Other")], na.rm = TRUE),
    proportion = ifelse(transcript_id == "Exon 4", TPM, 0) / total_TPM
  ) %>%
  ungroup() %>%
  filter(transcript_id == "Exon 4")

# Step 4: Define color mappings
color_df <- hist_indep_rna_df %>%
  select(plot_group_hex, plot_group) %>%
  filter(!is.na(plot_group)) %>%
  distinct()

color_mapping <- setNames(color_df$plot_group_hex, color_df$plot_group)

transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  mutate(
    dot_color = case_when(
      plot_group %in% names(color_mapping) ~ color_mapping[plot_group],
      TRUE ~ "#b5b5b5"
    )
  ) %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  mutate(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
  ungroup()

# Step 5: Statistical tests
plot_p_values_wilcoxon <- data.frame(cluster = character(), p_value = numeric(), stringsAsFactors = FALSE)
plot_p_values_permutation <- data.frame(cluster = character(), p_value = numeric(), stringsAsFactors = FALSE)

kruskal_results <- list()

for (group in unique(transcript_expr_CLK1_combined_df$cluster)) {
  group_data <- filter(transcript_expr_CLK1_combined_df, cluster == group)
  
  if (length(unique(group_data$age_bin)) > 1) {
    perm_result <- oneway_test(proportion ~ factor(age_bin), data = group_data, distribution = approximate(nresample = 9999))
    wilcox_result <- pairwise.wilcox.test(group_data$proportion, group_data$age_bin, p.adjust.method = "bonferroni")
    
    kruskal_results[[group]] <- list(
      permutation_p_value = pvalue(perm_result),
      wilcoxon_p_values = wilcox_result$p.value
    )
    
    plot_p_values_wilcoxon <- rbind(plot_p_values_wilcoxon, data.frame(cluster = group, p_value = min(wilcox_result$p.value, na.rm = TRUE)))
    plot_p_values_permutation <- rbind(plot_p_values_permutation, data.frame(cluster = group, p_value = pvalue(perm_result)))
  } else {
    warning(sprintf("Only one group in age_bin for %s", group))
  }
}

# Step 6: Visualization
create_plot <- function(data, title, p_values, value_label) {
  ggplot(data, aes(x = age_bin, y = proportion)) +
    geom_jitter(aes(color = dot_color), width = 0.2, size = 2) +
    geom_boxplot(aes(group = age_bin), width = 0.6, color = "black", fill = "white", alpha = 0.2) +
    labs(title = title, x = "Age Group", y = "Exon 4 Isoform Fraction") +
    scale_color_identity(name = "Dot Color") +
    facet_wrap(~cluster, scales = "free_x", ncol = 4) +
    theme_Publication() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 75, hjust = 1)) +
    geom_text(
      data = p_values,
      aes(x = 1, y = max(data$proportion, na.rm = TRUE) * 1.05, label = paste0(value_label, " p = ", signif(p_value, 4))),
      inherit.aes = FALSE, size = 3, hjust = 0
    )
}

tpm_plot_wilcoxon <- create_plot(transcript_expr_CLK1_combined_df, "Relative CLK1 Exon 4 Transcript Expression", plot_p_values_wilcoxon, "")
tpm_plot_permutation <- create_plot(transcript_expr_CLK1_combined_df, "Relative CLK1 Exon 4 Transcript Expression", plot_p_values_permutation, "")

# Save plots
pdf(file.path(plots_dir, "clk1ex4-clusters-tpm-tumor-age-bin-wc-test.pdf"), height = 14, width = 12)
print(tpm_plot_wilcoxon)
dev.off()

pdf(file.path(plots_dir, "clk1ex4-clusters-tpm-tumor-age-bin-perm-test.pdf"), height = 14, width = 12)
print(tpm_plot_permutation)
dev.off()

transcript_expr_CLK1_prop_df <- transcript_expr_CLK1_combined_df %>% 
  dplyr::select(Kids_First_Biospecimen_ID, plot_group,proportion)
write_tsv(x = transcript_expr_CLK1_prop_df,file = file.path(results_dir,"clk1-exon4-proportion.tsv"))