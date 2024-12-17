# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
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

gtex_trans_file <- file.path("/Users/naqvia/d3b_coding/neoepitope-identification/data/gtex-harmonized-isoform-expression-rsem-tpm.rds")
ped_trans_file = "~/d3b_coding/neoepitope-identification/data/GSE243682_normal_rna-isoform-expression-rsem-tpm.rds"

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

# Step 1: Load and filter clinical and independent data
indep_df <- read_tsv(indep_file)

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
  mutate(
    age_years = age_at_diagnosis_days / 365,
    age_bin = cut(
      age_years,
      breaks = c(0, 14, 18, 39),
      labels = c("[1-14]", "[15-18]", "[19-39]"),
      right = TRUE
    )
  ) %>%
  select(Kids_First_Biospecimen_ID, plot_group, age_bin, transcript_id, TPM) %>%
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
  filter(!is.na(plot_group)) %>%
  group_by(plot_group) %>%
  mutate(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(plot_group = factor(
    plot_group,
    levels = transcript_expr_CLK1_combined_df %>%
      arrange(plot_group, desc(mean_proportion)) %>%
      pull(plot_group) %>%
      unique()
  ))

# Step 5: Statistical tests
plot_p_values_wilcoxon <- data.frame(plot_group = character(), p_value = numeric(), stringsAsFactors = FALSE)
plot_p_values_permutation <- data.frame(plot_group = character(), p_value = numeric(), stringsAsFactors = FALSE)

kruskal_results <- list()

for (group in unique(transcript_expr_CLK1_combined_df$plot_group)) {
  group_data <- filter(transcript_expr_CLK1_combined_df, plot_group == group)
  
  if (length(unique(group_data$age_bin)) > 1) {
    perm_result <- oneway_test(proportion ~ factor(age_bin), data = group_data, distribution = approximate(nresample = 9999))
    wilcox_result <- pairwise.wilcox.test(group_data$proportion, group_data$age_bin, p.adjust.method = "bonferroni")
    
    kruskal_results[[group]] <- list(
      permutation_p_value = pvalue(perm_result),
      wilcoxon_p_values = wilcox_result$p.value
    )
    
    plot_p_values_wilcoxon <- rbind(plot_p_values_wilcoxon, data.frame(plot_group = group, p_value = min(wilcox_result$p.value, na.rm = TRUE)))
    plot_p_values_permutation <- rbind(plot_p_values_permutation, data.frame(plot_group = group, p_value = pvalue(perm_result)))
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
    facet_wrap(~plot_group, scales = "free_x", ncol = 5) +
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
pdf(file.path(plots_dir, "clk4-tpm-tumor-age-bin-wc-test.pdf"), height = 14, width = 12)
print(tpm_plot_wilcoxon)
dev.off()

pdf(file.path(plots_dir, "clk4-tpm-tumor-age-bin-perm-test.pdf"), height = 14, width = 12)
print(tpm_plot_permutation)
dev.off()