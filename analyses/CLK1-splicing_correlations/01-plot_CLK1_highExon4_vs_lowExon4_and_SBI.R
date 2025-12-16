################################################################################
# 01-plot_CLK1_highExon4_vs_lowExon4_and_SBI.R
# written by Ammar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 01-plot_CLK1_highExon4_vs_lowExon4_and_SBI.R
################################################################################

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
figures_dir <- file.path(root_dir, "figures")
input_dir   <- file.path(root_dir, "analyses", "splicing_index", "results")

# Specify file paths
sbi_file <-  file.path(input_dir,"splicing_index.total.txt")
clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
rmats_file <- file.path(results_dir, "clk1-splice-events-rmats.tsv")

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

## Load clinical file
hist_rna_df  <-  read_tsv(clin_file) %>%
  # Select only "RNA-Seq" samples
  filter(experimental_strategy == "RNA-Seq")

hgg_bs_id <- hist_rna_df %>%
  filter(plot_group %in% c("Diffuse midline glioma", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

dmg_bs_id <- hist_rna_df %>%
  filter(plot_group == "Diffuse midline glioma") %>%
  pull(Kids_First_Biospecimen_ID)

other_hgg_bs_id <- hist_rna_df %>%
  filter(plot_group == "Other high-grade glioma") %>%
  pull(Kids_First_Biospecimen_ID)

all_bs_id <- hist_rna_df %>%
  pull(Kids_First_Biospecimen_ID)

## Load rmats file
rmats_df <-  vroom(rmats_file) %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  dplyr::select(sample_id, geneSymbol, IncLevel1) %>%
  # Join rmats data with clinical data
  inner_join(hist_rna_df, by=c('sample_id'='Kids_First_Biospecimen_ID')) 

## Load SBI file from previous module
sbi_vs_inclEx4_df <-  vroom(sbi_file) %>%
  inner_join(rmats_df, by=c("Sample"="sample_id"))

set.seed(45)

bs_list <- list("all_hgg" = hgg_bs_id, "dmg" = dmg_bs_id, "other_hgg" = other_hgg_bs_id, "all" = all_bs_id)
names <- names(bs_list)

# Loop through groups
for (each in names) {

  # Filter the DataFrame based on current group's IDs
  new_df <- sbi_vs_inclEx4_df %>%
    filter(Sample %in% bs_list[[each]])
  ## Compute quantiles to define high vs low Exon 4 PSI groups
  quartiles_psi <- quantile(new_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
  # Calculate IQR
  IQR_psi <- IQR(new_df$IncLevel1)
  
  # Get lower quantile (25%)
  lower_psi <- quartiles_psi[1] 
  # Get upper quantile (75%)
  upper_psi <- quartiles_psi[2]
  
  # Create df with high/low PSI
  sbi_vs_inclEx4_by_extremePSI_df <- new_df %>%
    mutate(PSI = case_when(new_df$IncLevel1 > upper_psi ~ "high",
                           new_df$IncLevel1 < lower_psi ~ "low",
                           TRUE ~ NA_character_)
           ) %>%
    filter(!is.na(PSI)) %>%
    mutate(PSI = factor(PSI, levels = c("high", "low"))) %>%
    dplyr::select(Sample,SI,PSI)
  
  # add a little extra room for wilcoxon test
  ylim_max <- max(sbi_vs_inclEx4_by_extremePSI_df$SI)+0.10*(max(sbi_vs_inclEx4_by_extremePSI_df$SI))
  
  ## Make box plot with stats
  boxplot_sbi_vs_incl <- ggboxplot(sbi_vs_inclEx4_by_extremePSI_df,x = "PSI", y = "SI", outlier.shape = NA) +
    xlab(expression(bold(bolditalic("CLK1")~"Exon 4 PSI Level"))) +
    ylab(expression(bold("Splicing Burden Index"))) +
    lims(y = c(0,ylim_max)) +
    ggforce::geom_sina(aes(color = PSI), alpha = 4, size = 1.5, fill = NA) +
    scale_color_manual(name = "PSI Level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
    stat_compare_means(position = "identity", label.x = 1, size = 4) +
    theme_Publication()

  # Save plot pdf
  pdf(file.path(plots_dir, paste0(each, "_SBI_high_vs_low_CLK1_stranded.pdf")), height = 4.5, width = 4.5)
  print(boxplot_sbi_vs_incl)
  dev.off()
}
