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
  library(ggVennDiagram)
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
#depmap_model_file <- file.path(input_dir, "Model.csv")

#depmap_models <- vroom(depmap_model_file, show_col_types = FALSE)
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>%   
  dplyr::rename(
    ModelID = `Depmap ID`,
    `Cell line` = `Cell Line Name`,
    Histology = `Lineage Subtype`,
    CRISPR_Score = `CRISPR (DepMap Public 23Q2+Score, Chronos)`
  ) #%>%
  #left_join(depmap_models) %>%
  #filter(Age < 40) 
  
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

# Extract gene names from transcript IDs
extract_gene <- function(transcript) {
  sapply(strsplit(transcript, " "), `[`, 1)
}

# Prepare datasets with gene names
cor_myeloid_with_gene <- cor_myeloid_filtered %>%
  mutate(gene = extract_gene(transcript))

cor_cns_with_gene <- cor_cns_filtered %>%
  mutate(gene = extract_gene(transcript))

# =============================================================================
# Venn Diagram Visualization
# =============================================================================

venn_diag <- ggVennDiagram(x = list(unique(cor_myeloid_with_gene$gene), unique(cor_cns_with_gene$gene)),
                           edge_lty = "dashed",
                           edge_size = 1,
                           label_size = 6,
                           set_size = 5,
                           category.names = c("Myeloid", "CNS/Brain"),
                           label_percent_digit = 1) +
  scale_fill_gradient(
    low = "#ffffff",
    high = "steelblue1",
    name = expression(bold("Gene count")),
    limits = c(0, NA)
  ) +
  labs(title = expression(bold("Correlations with CLK1 (|r| > 0.4)"))) +
  # keep normal left/right orientation; also prevent clipping of labels
  coord_flip(clip = "off") +
  theme(
    # extra left margin so label never gets cut off
    plot.margin = margin(t = 0, r = 6, b = 0, l = 18, unit = "mm"),
    # remove axes completely
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank()
  )
  
## remove gray boxes: convert geom_label -> geom_text
venn_diag$layers <- lapply(venn_diag$layers, function(l) {
  if (inherits(l$geom, "GeomLabel")) {
    l$geom <- ggplot2::GeomText
  }
  l
})
  
ggplot2::ggsave(
  file.path(plots_dir, "CLK1_correlation_venn.pdf"),
  plot = venn_diag,
  width = 6.5,
  height = 3.5,
  device = "pdf",
  dpi = 300
)

# =============================================================================
# Create Comprehensive Gene Table
# =============================================================================

# Get all unique transcripts from both datasets
all_transcripts <- unique(c(
  cor_results_myeloid$transcript,
  cor_results_cns$transcript
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
genes_both         <- gene_table %>% filter(cell_lines == "Both") %>% distinct(gene) %>% pull(gene)
genes_cns_only     <- gene_table %>% filter(cell_lines == "CNS/Brain")  %>% distinct(gene) %>% pull(gene)
genes_myeloid_only <- gene_table %>% filter(cell_lines == "Myeloid") %>% distinct(gene) %>% pull(gene)
genes_cns          <- gene_table %>% filter(cell_lines %in% c("CNS/Brain", "Both"))  %>% distinct(gene) %>% pull(gene)
genes_myeloid      <- gene_table %>% filter(cell_lines %in% c("Myeloid", "Both")) %>% distinct(gene) %>% pull(gene)

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

plot_go_dot <- function(ego_object, top_n = 15, title = "GO Enrichment", plot_path = "", height = 6, width = 9, label_char = 50) {
 
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

ego_both         <- run_go_enrichment(genes_both, universe_ids = background_ids)
ego_cns          <- run_go_enrichment(genes_cns, universe_ids = background_ids)
ego_cns_only     <- run_go_enrichment(genes_cns_only, universe_ids = background_ids)
ego_myeloid      <- run_go_enrichment(genes_myeloid, universe_ids = background_ids)
ego_myeloid_only <- run_go_enrichment(genes_myeloid_only, universe_ids = background_ids)

# Generate plots
plot_go_dot(ego_both, title = "GO Enrichment — Shared (Both lineages)", plot_path = file.path(plots_dir,"GO_BP_dotplot_Both.pdf"))
plot_go_dot(ego_cns, title = "GO Enrichment — CNS", plot_path = file.path(plots_dir,"GO_BP_dotplot_CNS.pdf"))
plot_go_dot(ego_cns_only, title = "GO Enrichment — CNS Only", plot_path = file.path(plots_dir,"GO_BP_dotplot_CNS_only.pdf"))
plot_go_dot(ego_myeloid, title = "GO Enrichment — Myeloid", plot_path = file.path(plots_dir,"GO_BP_dotplot_Myeloid.pdf"))
plot_go_dot(ego_myeloid_only, title = "GO Enrichment — Myeloid Only", plot_path = file.path(plots_dir,"GO_BP_dotplot_Myeloid_only.pdf"))

