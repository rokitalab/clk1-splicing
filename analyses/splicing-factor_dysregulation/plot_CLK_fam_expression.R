## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  library(ggplot2)
  library('tidyverse')
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


hgg_bs_id <- clin_tab %>%
  # Select only "RNA-Seq" samples
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

sf_list <- c("CLK1", "CLK2", "CLK3", "CLK4")

file_gene_counts = "/Users/naqvia/d3b_coding/clk1-splicing/data/gene-expression-rsem-tpm-collapsed.rds"

count_data <- readRDS(file_gene_counts) %>% 
  #filter for HGG midline samples stranded and high sbi
  dplyr::select(any_of(hgg_bs_id)) 

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
stat.test <- stat.test %>%
  mutate(y.position = max(data_long$expression) * seq(1.05, 1.5, length.out = nrow(stat.test)))

# Plot with pairwise comparison results and mean labels
boxplot_tpm<- ggplot(data_long, aes(x = gene, y = expression)) +
  geom_boxplot(outlier.shape = NA, fill = "blue", color = "#2c3e50", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -0.5, color = "red", size = 4) + 
  stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = TRUE, tip.length = 0.01) + # Adjust if column is different
  labs(title = "CLK Family Expression", 
       x = "Gene", y = "TPM") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, color = "#2c3e50"),
    axis.text.x = element_text(color = "#2c3e50"),
    axis.text.y = element_text(color = "#2c3e50")
  ) + theme_Publication()