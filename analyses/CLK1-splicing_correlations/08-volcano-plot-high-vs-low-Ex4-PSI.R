################################################################################
# 08-volcano-plot-high-vs-low-Ex4-PSI.R
# written by Ammar Naqvi 
#
# usage: Rscript 08-volcano-plot-high-vs-low-Ex4-PSI.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(vroom)
  library(data.table)
  library(DESeq2)
  library(EnhancedVolcano)
  library("msigdbr")
  library("org.Hs.eg.db")
  library("DOSE")
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
figures_dir <- file.path(root_dir, "figures")
input_dir   <- file.path(root_dir, "analyses", "splicing_index", "results")


# Specify file paths
sbi_file <-  file.path(input_dir,"splicing_index.SE.txt")
clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")
rmats_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
file_gene_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")  

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

# Path to save file as
hgg_plot_file <- file.path(plots_dir,"all_hgg_boxplot_SBI_high_vs_low_CLK1.pdf")
dmg_plot_file <- file.path(plots_dir,"dmg_boxplot_SBI_high_vs_low_CLK1.pdf")
other_hgg_plot_file <- file.path(plots_dir,"other_hgg_boxplot_SBI_high_vs_low_CLK1.pdf")

## Load clinical file
indep_df <- read_tsv(indep_file)

hist_rna_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

hgg_bs_id <- hist_rna_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

dmg_bs_id <- hist_rna_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group == "DIPG or DMG") %>%
  pull(Kids_First_Biospecimen_ID)

other_hgg_bs_id <- hist_rna_df %>%
  filter(plot_group == "Other high-grade glioma") %>%
  pull(Kids_First_Biospecimen_ID)

# List of sample ID groups and corresponding names for output files
sample_id_groups <- list(hgg_bs_id = hgg_bs_id, dmg_bs_id = dmg_bs_id, other_hgg_bs_id = other_hgg_bs_id)
output_names <- c("All", "DMGs", "Other HGGs")

# Function to process each sample ID group and save the plot
process_rmats_data <- function(sample_ids, name) {
  rmats_df <- vroom(rmats_file) %>%
    filter(geneSymbol == "CLK1", exonStart_0base == "200860124", exonEnd == "200860215") %>%
    dplyr::select(sample_id, geneSymbol, IncLevel1) %>%
    inner_join(hist_rna_df, by = c('sample_id' = 'Kids_First_Biospecimen_ID')) %>%
    filter(sample_id %in% sample_ids)
  
  # List of sample_ids to match for "low" low ones that are differentially spliced
  low_sample_ids <- c(
    "BS_0XHT9W4Q", "BS_17ZSMXFH", "BS_2EN3X6HB", "BS_7PVJ6Q3D",
    "BS_9Y1K14Q5", "BS_FYP2B298", "BS_MD0937XT", "BS_MDBT7S5Z",
    "BS_N6TQBQ8M", "BS_NPMEGFW8", "BS_TRPXB3AV", "BS_XM1AHBDJ"
  )
  
  # Step 1: Label entries with matching sample IDs as "low"
  rmats_labeled_df <- rmats_df %>%
    mutate(PSI = if_else(sample_id %in% low_sample_ids, "low", NA_character_))
  
  # Step 2: Count the number of "low" entries
  n_lows <- rmats_labeled_df %>%
    filter(PSI == "low") %>%
    nrow()
  
  # Step 3: Select the top `n_lows` entries with the highest IncLevel1
  top_high_samples <- rmats_labeled_df %>%
    arrange(desc(IncLevel1)) %>%
    slice_head(n = n_lows) %>%
    mutate(PSI = "high")
  
  # Step 4: Combine the "low" and "high" labels
  rmats_final_df <- rmats_labeled_df %>%
    left_join(top_high_samples %>% dplyr::select(sample_id, IncLevel1, PSI), 
              by = c("sample_id", "IncLevel1"),
              suffix = c("", "_high")) %>%
    mutate(PSI = coalesce(PSI_high, PSI)) %>%
    dplyr::select(-PSI_high) %>%
    filter(!is.na(PSI))
  
  #print(rmats_final_df)
  
  count_data <- readRDS(file_gene_counts) %>%
    dplyr::select(any_of(rmats_final_df$sample_id)) %>%
    mutate(gene = rownames(.)) %>%
    rowwise() %>%
    # Keep rows where every count in the selected samples is at least 10
    filter(min(c_across(where(is.numeric))) >= 2) %>%
    ungroup()
  
  design <- data.frame(row.names = rmats_final_df$sample_id, condition = rmats_final_df$PSI)
  cds <- DESeqDataSetFromMatrix(countData = round(dplyr::select(count_data, -gene)), colData = design, design = ~ condition)
  cds <- DESeq(cds)
  
  res <- results(cds) %>%
    as_tibble() %>%
    mutate(Significant = ifelse(padj < 0.05, "P-Adj < 0.05", "Not Sig")) %>%
    add_column(gene = count_data$gene)
  
  max_padj_log <- max(-log10(res$padj), na.rm = TRUE)
  
  volc_hgg_plot <- EnhancedVolcano(res,
                                   lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "", res$gene),
                                   x = 'log2FoldChange', y = 'padj',
                                   title = paste('High vs Low', name, 'Ex4 HGGs'),
                                   subtitle = NULL,
                                   caption = NULL,
                                   pCutoff = 0.05, 
                                   FCcutoff = 2,
                                   labSize = 2.5,
                                   ylim = c(0, max_padj_log)) # Set y-axis limits
  
  volc_hgg_plot <- volc_hgg_plot + labs(x = expression(bold(Log[2] * " Fold Change")),
                                        y = expression(bold("-Log"[10] * " p-value")))
  
  
  # Save significant genes to file
  res %>%
    filter(padj < 0.05) %>%
    write_tsv(file.path(results_dir, paste0(name, "_gene_sign_list.tsv")))
  
  # Save the plot
  ggsave(filename = file.path(plots_dir, paste0(name, "_volcano_plot.pdf")),
         plot = volc_hgg_plot, width = 8, height = 6)
  
  
  library("clusterProfiler")
  
  ## outplut file for plot
  ora_dotplot_func_path <- file.path(plots_dir, paste0("high_low_ex4_diff-genes-ora-dotplot-",name,".pdf"))
  
  ## get gene sets relevant to H. sapiens
  hs_msigdb_df <- msigdbr(species = "Homo sapiens")
  pathway_df <- hs_msigdb_df %>%
    dplyr::filter(gs_cat == "H")
  
  res_fc2 <- res %>% filter(abs(log2FoldChange) >= 2)
  
  ## run enrichR to compute and identify significant over-repr pathways
  ora_results <- enricher(
    gene = res_fc2$gene, # A vector of your genes of interest
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE = dplyr::select(
      pathway_df,
      gs_name,
      human_gene_symbol
    )
  )
  
  ora_result_df <- data.frame(ora_results@result)
  options(enrichplot.colours = c("darkorange","blue"))
  enrich_plot_func <- enrichplot::dotplot(ora_results,
                                          x = "geneRatio",
                                          size = "Count",
                                          color = "p.adjust",
                                          label_format = 30,
                                          showCategory = 20) +
    labs(y = "Pathway",
         x = "Gene Ratio") +
    theme_Publication() +
    scale_size(name = "Gene Count") +
    scale_fill_gradient(low = "darkorange", high = "blue", name = "B-H p-value") +
    guides(
      fill = guide_colorbar(title = "B-H p-value", label.position = "right", barwidth = 1, barheight = 4)
    )
  
  ggplot2::ggsave(ora_dotplot_func_path,
                  plot=enrich_plot_func,
                  width=7.5,
                  height=4,
                  device="pdf",
                  dpi=300)
  detach("package:clusterProfiler", unload = TRUE)
  
  return(volc_hgg_plot)
}

# Apply the function to each group
plots <- purrr::map2(sample_id_groups, output_names, process_rmats_data)