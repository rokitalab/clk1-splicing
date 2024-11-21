################################################################################
# 08-plot-Ex4PSI-correlations.R
# written by Ammar Naqvi 
#
# usage: Rscript 08-plot-Ex4PSI-correlations.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library("clusterProfiler")
  
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
    select(sample_id, geneSymbol, IncLevel1) %>%
    inner_join(hist_rna_df, by = c('sample_id' = 'Kids_First_Biospecimen_ID')) %>%
    filter(sample_id %in% sample_ids)
  
  quartiles_psi <- quantile(rmats_df$IncLevel1, probs = c(.25, .75), na.rm = TRUE)
  rmats_levels_df <- rmats_df %>%
    mutate(PSI = case_when(
      IncLevel1 > quartiles_psi[2] ~ "high",
      IncLevel1 < quartiles_psi[1] ~ "low",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(PSI))
  
  count_data <- readRDS(file_gene_counts) %>%
    select(any_of(rmats_levels_df$sample_id)) %>%
    mutate(gene = rownames(.)) %>%
    rowwise() %>%
    # Keep rows where every count in the selected samples is at least 10
    filter(min(c_across(where(is.numeric))) >= 2) %>%
    ungroup()
 
  design <- data.frame(row.names = rmats_levels_df$sample_id, condition = rmats_levels_df$PSI)
  cds <- DESeqDataSetFromMatrix(countData = round(select(count_data, -gene)), colData = design, design = ~ condition)
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
  return(volc_hgg_plot)
}

# Apply the function to each group
plots <- purrr::map2(sample_id_groups, output_names, process_rmats_data)


### CLK1 exon 4 correlates that are SRSF targets
rmats_clk1_df <- vroom(rmats_file) %>%
  filter(geneSymbol == "CLK1", exonStart_0base == "200860124", exonEnd == "200860215") %>%
  dplyr::select(sample_id, geneSymbol, IncLevel1) %>%
  inner_join(hist_rna_df, by = c('sample_id' = 'Kids_First_Biospecimen_ID')) %>%
  filter(sample_id %in% hgg_bs_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=sample_id) %>%
  dplyr::select(Kids_First_Biospecimen_ID,IncLevel1) %>%
  dplyr::rename("CLK1_PSI"=IncLevel1)


SRSF_targets <- c("BIN1", "RON", "BCL2L1", "MKNK2", "CASP9",
                  "CD45", "BCL2L1", "TP53", "FAS", "MCL1",
                  "TPM1", "SYN1", "TAU", "ERBB2", "SMN2",
                  "FAS", "MYC", "VEGF", "MCL1", "CASP2",
                  "FN1", "DMD", "CD44", "CASP8", "KLF6")



# View the filtered results
rmats_full_file = file.path(data_dir, "splice-events-rmats.tsv.gz")

rmats_full_df <- vroom(rmats_full_file) 
rmats_se_full_df <- vroom("~/d3b_coding/neoepitope-identification/data/pbta-splice-events-rmats.SE.tsv.gz")

rmats_clk1_targets <- rmats_se_full_df %>% 
  filter(geneSymbol %in% SRSF_targets) %>%  # Filter for genes in all_genes
  mutate(
    SpliceID = paste0(
      geneSymbol, ":",
      upstreamES, "-", upstreamEE, "_",
      exonStart_0base, "-", exonEnd, "_",
      downstreamES, "-", downstreamEE
    )  # Combine fields into SpliceID
  ) %>%
  dplyr::select(sample_id, geneSymbol, SpliceID, IncLevel1) %>%
  inner_join(hist_rna_df, by = c("sample_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(sample_id %in% hgg_bs_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id) %>%
  dplyr::select(Kids_First_Biospecimen_ID, SpliceID, IncLevel1)

# Check the resulting data frame
joined_df <- inner_join(rmats_clk1_targets,rmats_clk1_df, by='Kids_First_Biospecimen_ID')

# Directory to save plots
# Function to process each SpliceID
process_splice_id <- function(splice_id, data, output_dir) {
  # Filter data for the current SpliceID
  filtered_df <- data %>%
    filter(SpliceID == splice_id) %>%
    drop_na(CLK1_PSI, IncLevel1) # Remove rows with missing values
  
  # Check if there are enough data points for correlation analysis
  if (nrow(filtered_df) >= 2) {
    # Perform correlation test
    cor_result <- tryCatch({
      cor.test(filtered_df$CLK1_PSI, filtered_df$IncLevel1, use = "complete.obs")
    }, error = function(e) NULL)
    
    if (!is.null(cor_result) && 
        !is.na(cor_result$p.value) && 
        !is.na(cor_result$estimate)) {
      # Check if correlation meets criteria
      if (cor_result$p.value < 0.05 && abs(cor_result$estimate) > 0.4) {
        # Create scatter plot
        scatter_plot <- ggplot(filtered_df, aes(x = CLK1_PSI, y = IncLevel1)) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", color = "red", se = TRUE, formula = y ~ x) +
          labs(
            title = paste("Scatter Plot for", splice_id),
            subtitle = paste("r =", round(cor_result$estimate, 3), 
                             ", p-value =", signif(cor_result$p.value, 3)),
            x = "CLK1 PSI",
            y = "PSI"
          ) +
          theme_minimal() + 
          theme_Publication()
        
        # Save the plot
        file_name <- paste0(output_dir, "/", gsub("[:/]", "_", splice_id), ".pdf")
        ggsave(file_name, scatter_plot, width = 8, height = 6)
        
        return(scatter_plot)
      } else {
        message(paste("Correlation for", splice_id, "does not meet criteria: p-value >= 0.05 or |r| <= 0.4"))
        return(NULL)
      }
    } else {
      message(paste("Correlation test failed or produced NA for SpliceID:", splice_id))
      return(NULL)
    }
  } else {
    message(paste("Not enough data points for SpliceID:", splice_id))
    return(NULL)
  }
}

# Apply the function to each unique SpliceID
unique_splice_ids <- unique(joined_df$SpliceID)
plots <- map(unique_splice_ids, ~process_splice_id(.x, joined_df, plots_dir))
