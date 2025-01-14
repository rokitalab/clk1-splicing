################################################################################
# 11-plot-clk1-ex4-transcripts-by-age-gtex
# Script that plots CLK1 exon 4 PSI variations across gtex ctrls by age
# written by Ammar Naqvi
#
# usage: Rscript 11-plot-clk1-ex4-transcripts-by-age-gtex
################################################################################

library(dplyr)
library(ggplot2)
library(vroom)
library(coin)  # For oneway_test
library(tidyr)
library(stringr)
library(purrr)

## Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Define file paths
histology_gtex_file <- file.path(input_dir,"gtex-samples-by-age.tsv")
gtex_trans_file <- file.path(data_dir,"gtex-harmonized-isoform-expression-rsem-tpm.rds")

# Ensure the output directory exists
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

# Load histology data and filter by age and tissue
histology_gtex_df <- vroom(histology_gtex_file)

histology_brain_age_filtered_df <- histology_gtex_df %>%
  filter(gtex_group == "Brain") %>%
  mutate(gtex_subgroup = gsub("Brain - ", "", gtex_subgroup))

# Load GTEx TPM data and filter for CLK1 isoforms
gtex_clk1_transc_counts <- readRDS(gtex_trans_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("GTEX"), sum, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Kids_First_Biospecimen_ID",
    values_to = "TPM"
  ) %>%
  inner_join(histology_brain_age_filtered_df, by = "Kids_First_Biospecimen_ID") %>%
  rename(plot_group = gtex_subgroup) %>%
  select(transcript_id, Kids_First_Biospecimen_ID, plot_group, AGE, TPM) %>%
  rename(age_bin = AGE) %>%
  mutate(
    group = "Gtex",
    plot_group = str_replace_all(plot_group, "Brain - ", "")
  )

# Compute Exon 4 proportions
transcript_expr_CLK1_gtex_df <- gtex_clk1_transc_counts %>%
  group_by(Kids_First_Biospecimen_ID) %>%
  mutate(total_TPM = sum(TPM[transcript_id %in% c("Exon 4", "Other")], na.rm = TRUE)) %>%
  mutate(proportion = ifelse(transcript_id == "Exon 4", TPM, 0) / total_TPM) %>%
  ungroup() %>%
  filter(transcript_id == "Exon 4")

# Initialize results storage
plot_p_values_wilcoxon <- data.frame(
  plot_group = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
plot_p_values_permutation <- data.frame(
  plot_group = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
kruskal_results <- list()

# Perform Wilcoxon and Permutation tests for each plot_group
for (group2 in unique(transcript_expr_CLK1_gtex_df$plot_group)) {
  group_data <- subset(transcript_expr_CLK1_gtex_df, plot_group == group2)
  
  if (length(unique(group_data$age_bin)) > 1) {
    # Permutation Test
    perm_result <- oneway_test(
      proportion ~ factor(age_bin),
      data = group_data,
      distribution = approximate(nresample = 9999)
    )
    
    # Pairwise Wilcoxon Test
    wilcox_result <- pairwise.wilcox.test(
      x = group_data$proportion,
      g = group_data$age_bin,
      p.adjust.method = "bonferroni"
    )
    
    # Store results
    kruskal_results[[group2]] <- list(
      permutation_p_value = pvalue(perm_result),
      wilcoxon_p_values = wilcox_result$p.value,
      plot_group = group2
    )
    
    # Extract smallest p-value
    wilcoxon_p_val <- min(wilcox_result$p.value, na.rm = TRUE)
    plot_p_values_wilcoxon <- rbind(plot_p_values_wilcoxon, data.frame(
      plot_group = group2,
      p_value = wilcoxon_p_val
    ))
    permutation_p_val <- pvalue(perm_result)
    plot_p_values_permutation <- rbind(plot_p_values_permutation, data.frame(
      plot_group = group2,
      p_value = permutation_p_val
    ))
  } else {
    warning(paste("Only one group in age_bin for", group2))
  }
}

# Calculate max y-values for plot annotations
max_proportion <- transcript_expr_CLK1_gtex_df %>%
  group_by(plot_group) %>%
  summarise(y_max = max(proportion, na.rm = TRUE) * 1.05)

# Wilcoxon test plot
tpm_plot_wilcoxon <- ggplot(transcript_expr_CLK1_gtex_df, aes(x = age_bin, y = proportion)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(aes(group = age_bin), width = 0.6, color = "darkgrey", fill = "white", alpha = 0.2) +
  labs(
    title = "",
    x = "Age Group",
    y = "Exon 4 Isoform Fraction"
  ) +
  facet_wrap(~ plot_group, scales = "free_x", nrow = 3) +
  theme_Publication() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 75, hjust = 1),
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
  ) +
  geom_text(
    data = plot_p_values_wilcoxon %>% left_join(max_proportion, by = "plot_group"),
    aes(
      x = 1.5,
      y = y_max,
      label = paste0("p = ", signif(p_value, digits = 4))
    ),
    inherit.aes = FALSE,
    size = 3
  )

# Permutation test plot
tpm_plot_permutation <- ggplot(transcript_expr_CLK1_gtex_df, aes(x = age_bin, y = proportion)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(aes(group = age_bin), width = 0.6, color = "darkgrey", fill = "white", alpha = 0.2) +
  labs(
    title = "",
    x = "Age Group",
    y = "Exon 4 Isoform Fraction"
  ) +
  facet_wrap(~ plot_group, scales = "free_x", nrow = 3) +
  theme_Publication() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 75, hjust = 1),
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
  ) +
  geom_text(
    data = plot_p_values_permutation %>% left_join(max_proportion, by = "plot_group"),
    aes(
      x = 1.5,
      y = y_max,
      label = paste0("p = ", signif(p_value, digits = 4))
    ),
    inherit.aes = FALSE,
    size = 3
  )

# Save plots to files
pdf(file.path(plots_dir, "clk4-tpm-gtex-age-bin-wc-test.pdf"), height = 10, width = 12)
print(tpm_plot_wilcoxon)
dev.off()

pdf(file.path(plots_dir, "clk4-tpm-gtex-age-bin-perm-test.pdf"), height = 10, width = 12)
print(tpm_plot_permutation)
dev.off()

