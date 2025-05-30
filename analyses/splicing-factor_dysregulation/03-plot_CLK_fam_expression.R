################################################################################
# 03-plot_CLK_fam_expression.R
# script that plots splicing factor related kinases (CLK1,CLK2,CLK3, CLK4, SRPK)
# that hit functional sites that were differential amongst HGGs
#
# written by Ammar Naqvi
# usage: Rscript 03-plot_CLK_fam_expression.R
################################################################################

# Libraries
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  library("ggplot2")
  library('tidyverse')
  library("rstatix")
  library("ggpubr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))


## get and setup input files
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- read_tsv(indep_file)

clin_file <- file.path(data_dir, "histologies-plot-group.tsv")
clin_tab  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         RNA_library == 'stranded',
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

clust_file <- file.path(root_dir, "analyses/sample-psi-clustering/results/sample-cluster-metadata-top-5000-events-stranded.tsv")
clust_df <- read_tsv(clust_file)

hgg_bs_id <- clin_tab %>%
  # Select only "RNA-Seq" samples
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

cluster6_bs_id <- clust_df %>%
  filter(cluster == 6) %>%
  pull(sample_id)

sf_list <- c("CLK1", "CLK2", "CLK3", "CLK4","SRPK1")

file_gene_counts = file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")

bs_list <- list("all_hgg" = hgg_bs_id, "cluster6" = cluster6_bs_id)
names <- names(bs_list)

for (each in names) {
  count_data <- readRDS(file_gene_counts) %>%
    #filter for HGG midline samples stranded and high sbi
    dplyr::select(any_of(bs_list[[each]])) 

  # Add gene names as a column to count_data
  count_data <- count_data %>%
    mutate(gene = rownames(.))

  count_data_sf <- count_data %>%
    dplyr::filter(gene %in% sf_list) %>%
    dplyr::select(gene, any_of(clin_tab$Kids_First_Biospecimen_ID)) %>%
    dplyr::rowwise() %>%  # Ensure you use parentheses here
    dplyr::filter(sum(c_across(where(is.numeric))) >= 1) %>%
    dplyr::ungroup()

  # Reshape data to long format
  data_long <- count_data_sf %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

  # Perform pairwise comparisons using Wilcoxon test and adjust p-values using Benjamini-Hochberg
  # Add a y.position column for each pairwise comparison to avoid the error
  # Perform pairwise comparisons (e.g., comparing CLK1 with all other genes)
  stat.test <- data_long %>%
    filter(gene %in% c("CLK1", "CLK2", "CLK3", "CLK4", "SRPK1")) %>%
    compare_means(expression ~ gene, data = ., 
                method = "wilcox.test", 
                p.adjust.method = "bonferroni")

  # Manually add a y.position for each comparison with triple the spacing
  stat.test$y.position <- seq(from = max(data_long$expression) + 1, 
                                 by = 9, length.out = nrow(stat.test))  # Triple the spacing

  stat.test_clk1 <- stat.test %>%
    # only retain CLK1 comparisons
    filter(group1 == "CLK1")

  # Create the boxplot with only CLK1 comparisons and p-values
  boxplot_tpm <- ggplot(data_long, aes(x = gene, y = expression)) +
    geom_boxplot(outlier.shape = NA, fill = "grey", color = "#2c3e50") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    geom_jitter(width = 0.2, size = 2, shape = 21, color = "black", fill="grey") +
    labs(title = "Kinase Expression", 
       x = "Gene", y = "TPM") +
    theme_Publication() + 
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 75, hjust = 1)) +
    scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 30))) + 
    stat_pvalue_manual(stat.test_clk1, 
                       label = "p.adj", 
                       hide.ns = TRUE, 
                       tip.length = 0.01)

  # Save plot as PDF
  pdf(file.path(plots_dir, paste0(each,"-CLK1-tpms-CLK-SPRK1-kinases.pdf")), 
      width = 6, height = 5)
  print(boxplot_tpm)
  dev.off()
}
