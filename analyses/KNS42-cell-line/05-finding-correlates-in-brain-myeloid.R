################################################################################
# 05-finding-correlates-in-brain-myeloid.R
#
# Identify and compare CLK1-correlated transcripts across cell line types
# - Calculate transcript-level correlations with CLK1 expression (|r| > 0.4)
# - Include p-value testing with multiple testing correction
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
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
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
depmap_model_file <- file.path(input_dir, "Model.csv")

depmap_models <- vroom(depmap_model_file, show_col_types = FALSE)
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>%   
  dplyr::rename(
    ModelID = `Depmap ID`,
    `Cell line` = `Cell Line Name`,
    Histology = `Lineage Subtype`,
    CRISPR_Score = `CRISPR (DepMap Public 23Q2+Score, Chronos)`
  ) %>%
  left_join(depmap_models) %>%
  filter(Age < 40) 
  
# Filter for glioma cell lines
depmap_glioma <- depmap_data %>%
  filter(Histology == "Diffuse Glioma")

depmap_data_KNS42 <- depmap_glioma %>%
  filter(`Cell line` == "KNS42")



# =============================================================================
# Transcript Correlation Analysis with P-values
# =============================================================================

# Load omics data -------------------------------------------------------------
omics_mappings_file <- file.path(data_dir, "OmicsDefaultModelProfiles.csv")
tpm_file <- file.path(data_dir, "OmicsExpressionTranscriptsTPMLogp1Profile.csv")

# Function to calculate correlations with p-values ---------------------------
calculate_clk1_correlations <- function(lineage_name) {
  
  cat("\nProcessing", lineage_name, "cell lines...\n")
  
  # Map cell lines to omics profiles
  omics_mapping <- vroom(omics_mappings_file, show_col_types = FALSE) %>%
    inner_join(depmap_data, by = "ModelID") %>%
    filter(Lineage == lineage_name)
  
  n_samples <- nrow(omics_mapping)
  cat("Found", n_samples, "cell lines\n")
  
  # Need at least 3 samples for meaningful correlation
  if(n_samples < 3) {
    cat("Too few samples. Skipping.\n")
    return(NULL)
  }
  
  # Load expression data
  expr_data <- vroom(tpm_file, show_col_types = FALSE) %>%
    dplyr::rename(ProfileID = `...1`) %>%
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
    dplyr::select(ProfileID, transcript, expression) %>%
    pivot_wider(names_from = transcript, values_from = expression)
  
  # Get CLK1 expression
  clk1_expr <- expr_wide[["CLK1 (ENST00000321356)"]]
  
  # Calculate correlations with p-values
  cat("Calculating correlations with p-values...\n")
  
  # Get all transcript columns
  transcript_cols <- setdiff(colnames(expr_wide), "ProfileID")
  
  # Calculate correlations one by one
  cor_results <- map_dfr(transcript_cols, function(transcript_name) {
    
    # Get expression values for this transcript
    transcript_expr <- expr_wide[[transcript_name]]
    
    # Skip if zero variance
    if(sd(transcript_expr, na.rm = TRUE) == 0) {
      return(NULL)
    }
    
    # Try to calculate correlation
    tryCatch({
      test_result <- cor.test(transcript_expr, clk1_expr, method = "pearson")
      
      tibble(
        transcript = transcript_name,
        corr_with_CLK1 = as.numeric(test_result$estimate),
        p_value = test_result$p.value,
        n_samples = n_samples
      )
    }, error = function(e) {
      return(NULL)
    })
  })
  
  cat("Calculated", nrow(cor_results), "correlations\n")
  
  return(cor_results)
}

# Calculate correlations for both lineages ------------------------------------
cat("\n=== Calculating CLK1 Transcript Correlations ===\n")

cor_results_cns <- calculate_clk1_correlations("CNS/Brain")
cor_results_myeloid <- calculate_clk1_correlations("Myeloid")

# Filter and annotate results -------------------------------------------------
filter_and_annotate <- function(cor_results, threshold = 0.3, p_threshold = 0.05) {
  
  if(is.null(cor_results) || nrow(cor_results) == 0) {
    return(tibble())
  }
  
  cor_results %>%
    filter(
      abs(corr_with_CLK1) > threshold,
      p_value < p_threshold  # Use unadjusted p-value
    ) %>%
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
    cor_myeloid_filtered %>% 
      dplyr::select(transcript, myeloid_corr = corr_with_CLK1, myeloid_dir = direction,
             myeloid_pval = p_value, myeloid_n = n_samples),
    by = "transcript"
  ) %>%
  left_join(
    cor_cns_filtered %>% 
      dplyr::select(transcript, cns_corr = corr_with_CLK1, cns_dir = direction,
             cns_pval = p_value, cns_n = n_samples),
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
  dplyr::select(
    gene,
    transcript,
    cell_lines,
    overall_direction,
    myeloid_direction,
    myeloid_corr,
    myeloid_pval,
    myeloid_n,
    cns_direction,
    cns_corr,
    cns_pval,
    cns_n
  ) %>%
  arrange(cell_lines, gene)

# Save as tab-delimited file
write_tsv(
  gene_table,
  file.path(results_dir, "CLK1_correlated_genes_comprehensive.tsv")
)

## GO Enrichment

# Subset genes with cell_lines == "Both"
# Define a helper function to get Entrez IDs
get_entrez_ids <- function(genes) {
  bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
}

# Background universe: all genes tested
background_genes <- unique(gene_table$gene)
background_ids <- get_entrez_ids(background_genes)

# Subset genes
genes_both    <- gene_table %>% filter(cell_lines == "Both") %>% distinct(gene) %>% pull(gene)
genes_cns     <- gene_table %>% filter(cell_lines == "CNS/Brain")  %>% distinct(gene) %>% pull(gene)
genes_myeloid <- gene_table %>% filter(cell_lines == "Myeloid") %>% distinct(gene) %>% pull(gene)

run_go_enrichment <- function(gene_list, universe_ids, ont="BP") {
  enrichGO(
    gene          = get_entrez_ids(gene_list),
    universe      = universe_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}

plot_go_dot <- function(ego_object, top_n = 15, title = "GO Enrichment", plot_path = "", height = 6, width = 8, label_char = 30) {
 
  options(enrichplot.colours = c("darkorange","blue"))
  plot <- enrichplot::dotplot(ego_object,
                              x = "GeneRatio",
                              size = "Count",
                              color = "p.adjust",
                              label_format = label_char,
                              showCategory = top_n) +   
    labs(y = "Pathway",
         x = "Gene Ratio",
         title = title) +
    theme_Publication() +
    scale_size(name = "Gene Count") +  
    scale_fill_gradient(low = "darkorange", high = "blue", name = "B-H p-value") +
    guides(
      fill = guide_colorbar(title = "B-H p-value", label.position = "right", barwidth = 1, barheight = 4)
    )
  
  ggplot2::ggsave(plot_path,
                  plot=plot,
                  width=width,
                  height=height,
                  device="pdf",
                  dpi=300)
}

ego_both    <- run_go_enrichment(genes_both, universe_ids = background_ids)
ego_cns     <- run_go_enrichment(genes_cns, universe_ids = background_ids)
ego_myeloid <- run_go_enrichment(genes_myeloid, universe_ids = background_ids)

# Generate plots
plot_go_dot(ego_both, title = "GO Enrichment — Shared (Both lineages)", plot_path = file.path(plots_dir,"GO_BP_dotplot_Both.pdf"), height = 6, width = 9, label_char = 50)
plot_go_dot(ego_cns, title = "GO Enrichment — CNS-only", plot_path = file.path(plots_dir,"GO_BP_dotplot_CNS.pdf"), height = 4)
plot_go_dot(ego_myeloid, title = "GO Enrichment — Myeloid-only", plot_path = file.path(plots_dir,"GO_BP_dotplot_Myeloid.pdf"))
