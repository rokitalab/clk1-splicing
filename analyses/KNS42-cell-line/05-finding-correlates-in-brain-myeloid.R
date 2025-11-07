################################################################################
# 05-finding-correlates-in-brain-myeloid.R
#
# Identify and compare CLK1-correlated transcripts across cell line types
# - Calculate transcript-level correlations with CLK1 expression (|r| > 0.4)
# - Separate positive and negative correlations by cell lineage
# - Generate Venn diagrams showing overlap between CNS/Brain and Myeloid lines
# - Export comprehensive tables of correlated genes and transcripts
#
# Author: Ammar S Naqvi
################################################################################

# Load libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(ggpubr)
  library(ggthemes)
  library(ggforce)
  library(patchwork)
})

# Setup directories -----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "KNS42-cell-line")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

# Create output directories
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Source plotting theme
source(file.path(root_dir, "figures", "theme_for_plots.R"))

# Set seed
set.seed(2023)

# =============================================================================
# CLK1 CRISPR Dependency Analysis
# =============================================================================

# Load DepMap CRISPR data -----------------------------------------------------
depmap_file <- file.path(data_dir, "CLK1-CRISPR-DepMap-score.csv")
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>%
  rename(
    ModelID = `Depmap ID`,
    `Cell line` = `Cell Line Name`,
    Histology = `Lineage Subtype`,
    CRISPR_Score = `CRISPR (DepMap Public 23Q2+Score, Chronos)`
  )

# Filter for glioma cell lines
depmap_glioma <- depmap_data %>%
  filter(Histology == "Diffuse Glioma")

depmap_data_KNS42 <- depmap_glioma %>%
  filter(`Cell line` == "KNS42")

# Plot CLK1 dependency scores -------------------------------------------------
gene_score_plot <- ggplot(
  data = depmap_glioma,
  aes(x = reorder(`Cell line`, CRISPR_Score), y = CRISPR_Score)
) +
  geom_point(size = 3, colour = "gray89") +
  geom_point(size = 3, colour = "gray50", pch = 21) +
  geom_point(data = depmap_data_KNS42, colour = "red", size = 3) +
  geom_point(data = depmap_data_KNS42, colour = "black", size = 3, pch = 21) +
  labs(
    title = "CLK1 Perturbation",
    x = "Cell Line",
    y = "CRISPR Dependency Score"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8))

# Save plot
ggsave(
  file.path(plots_dir, "depmap_score_CLK1_vs_score_KNS42.pdf"),
  gene_score_plot,
  height = 4, width = 7.5
)

# =============================================================================
# Transcript Correlation Analysis
# =============================================================================

# Load omics data -------------------------------------------------------------
omics_mappings_file <- file.path(data_dir, "OmicsDefaultModelProfiles.csv")
tpm_file <- file.path(data_dir, "OmicsExpressionTranscriptsTPMLogp1Profile.csv")

# Function to calculate correlations ------------------------------------------
calculate_clk1_correlations <- function(lineage_name) {
  
  cat("\nProcessing", lineage_name, "cell lines...\n")
  
  # Map cell lines to omics profiles
  omics_mapping <- vroom(omics_mappings_file, show_col_types = FALSE) %>%
    inner_join(depmap_data, by = "ModelID") %>%
    filter(Lineage == lineage_name)
  
  cat("Found", nrow(omics_mapping), "cell lines\n")
  
  # Load expression data
  expr_data <- vroom(tpm_file, show_col_types = FALSE) %>%
    rename(ProfileID = `...1`) %>%
    inner_join(omics_mapping, by = "ProfileID")
  
  # Convert to long format
  expr_long <- expr_data %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "transcript",
      values_to = "expression"
    ) %>%
    mutate(gene = sapply(strsplit(transcript, " "), `[`, 1))
  
  # Convert back to wide format for correlation
  expr_wide <- expr_long %>%
    select(ProfileID, transcript, expression) %>%
    pivot_wider(names_from = transcript, values_from = expression)
  
  # Get CLK1 expression
  clk1_expr <- expr_wide[["CLK1 (ENST00000321356)"]]
  
  # Calculate correlations
  cor_results <- expr_wide %>%
    select(-ProfileID) %>%
    summarise(across(
      everything(),
      ~ cor(.x, clk1_expr, use = "pairwise.complete.obs")
    )) %>%
    pivot_longer(
      cols = everything(),
      names_to = "transcript",
      values_to = "corr_with_CLK1"
    )
  
  return(cor_results)
}

# Calculate correlations for both lineages ------------------------------------
cat("\n=== Calculating CLK1 Transcript Correlations ===\n")

cor_results_cns <- calculate_clk1_correlations("CNS/Brain")
cor_results_myeloid <- calculate_clk1_correlations("Myeloid")

# Filter and annotate results -------------------------------------------------
filter_and_annotate <- function(cor_results, threshold = 0.4) {
  cor_results %>%
    filter(abs(corr_with_CLK1) > threshold) %>%
    mutate(
      direction = case_when(
        corr_with_CLK1 > threshold ~ "positive",
        corr_with_CLK1 < -threshold ~ "negative",
        TRUE ~ NA_character_
      )
    )
}

cor_cns_filtered <- filter_and_annotate(cor_results_cns)
cor_myeloid_filtered <- filter_and_annotate(cor_results_myeloid)

# =============================================================================
# Venn Diagram Visualization
# =============================================================================

# Function to create Venn diagram ---------------------------------------------
create_venn_plot <- function(set1, set2, label1, label2, title, colors) {
  
  # Calculate overlaps
  overlap <- length(intersect(set1, set2))
  only1 <- length(setdiff(set1, set2))
  only2 <- length(setdiff(set2, set1))
  
  # Create circle data
  circle_data <- data.frame(
    x = c(-0.5, 0.5),
    y = c(0, 0),
    r = c(1, 1),
    label = c(label1, label2),
    fill = colors
  )
  
  # Create plot
  p <- ggplot(circle_data) +
    geom_circle(
      aes(x0 = x, y0 = y, r = r, fill = fill),
      alpha = 0.3, color = "black", size = 1.5
    ) +
    scale_fill_identity() +
    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(title = title)
  
  # Add count labels
  p <- p +
    annotate("text", x = -0.9, y = 0, label = only1,
             size = 10, fontface = "bold", color = colors[1]) +
    annotate("text", x = 0.9, y = 0, label = only2,
             size = 10, fontface = "bold", color = colors[2]) +
    annotate("text", x = 0, y = 0, label = overlap,
             size = 10, fontface = "bold", color = "#2C3E50")
  
  # Add category labels
  p <- p +
    annotate("text", x = -0.9, y = 1.5, label = label1,
             size = 5, fontface = "bold", color = colors[1]) +
    annotate("text", x = 0.9, y = 1.5, label = label2,
             size = 5, fontface = "bold", color = colors[2])
  
  return(p)
}

# Separate transcripts by direction -------------------------------------------
myeloid_pos <- cor_myeloid_filtered %>%
  filter(direction == "positive") %>%
  pull(transcript)

myeloid_neg <- cor_myeloid_filtered %>%
  filter(direction == "negative") %>%
  pull(transcript)

cns_pos <- cor_cns_filtered %>%
  filter(direction == "positive") %>%
  pull(transcript)

cns_neg <- cor_cns_filtered %>%
  filter(direction == "negative") %>%
  pull(transcript)

# Create Venn diagrams --------------------------------------------------------
venn_positive <- create_venn_plot(
  myeloid_pos, cns_pos,
  "Myeloid", "CNS/Brain",
  "Positive Correlations with CLK1 (r > 0.4)",
  colors = c("#E74C3C", "#E67E22")
)

venn_negative <- create_venn_plot(
  myeloid_neg, cns_neg,
  "Myeloid", "CNS/Brain",
  "Negative Correlations with CLK1 (r < -0.4)",
  colors = c("#3498DB", "#9B59B6")
)


## venn diagram plot
combined_venn_horizontal <- venn_positive | venn_negative

# Save Venn diagrams
ggsave(
  file.path(plots_dir, "CLK1_correlation_venn_horizontal.pdf"),
  combined_venn_horizontal,
  width = 14, height = 7
)

# =============================================================================
# Create Comprehensive Gene Table
# =============================================================================

# Extract gene names from transcript IDs
extract_gene <- function(transcript) {
  sapply(strsplit(transcript, " "), `[`, 1)
}

# Prepare datasets with gene names
cor_myeloid_with_gene <- cor_myeloid_filtered %>%
  mutate(gene = extract_gene(transcript))

cor_cns_with_gene <- cor_cns_filtered %>%
  mutate(gene = extract_gene(transcript))

# Get all unique transcripts from both datasets
all_transcripts <- unique(c(
  cor_myeloid_filtered$transcript,
  cor_cns_filtered$transcript
))

# Create comprehensive table
gene_table <- tibble(transcript = all_transcripts) %>%
  mutate(gene = extract_gene(transcript)) %>%
  left_join(
    cor_myeloid_filtered %>% select(transcript, myeloid_corr = corr_with_CLK1, myeloid_dir = direction),
    by = "transcript"
  ) %>%
  left_join(
    cor_cns_filtered %>% select(transcript, cns_corr = corr_with_CLK1, cns_dir = direction),
    by = "transcript"
  ) %>%
  mutate(
    # Determine which cell line(s) have this transcript
    cell_lines = case_when(
      !is.na(myeloid_corr) & !is.na(cns_corr) ~ "Both",
      !is.na(myeloid_corr) ~ "Myeloid",
      !is.na(cns_corr) ~ "CNS/Brain",
      TRUE ~ NA_character_
    ),
    
    # Overall direction (if in both, note if concordant)
    overall_direction = case_when(
      cell_lines == "Both" & myeloid_dir == cns_dir ~ paste0(myeloid_dir, " (concordant)"),
      cell_lines == "Both" & myeloid_dir != cns_dir ~ "discordant",
      cell_lines == "Myeloid" ~ myeloid_dir,
      cell_lines == "CNS/Brain" ~ cns_dir,
      TRUE ~ NA_character_
    ),
    
    # Specific directions per cell line
    myeloid_direction = ifelse(!is.na(myeloid_dir), myeloid_dir, "not_significant"),
    cns_direction = ifelse(!is.na(cns_dir), cns_dir, "not_significant")
  ) %>%
  # Reorder columns for clarity
  select(
    gene,
    transcript,
    cell_lines,
    overall_direction,
    myeloid_direction,
    myeloid_corr,
    cns_direction,
    cns_corr
  ) %>%
  arrange(cell_lines, gene)

# Save as tab-delimited file
write_tsv(
  gene_table,
  file.path(results_dir, "CLK1_correlated_genes_comprehensive.tsv")
)

