# Libraries
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  library("ggplot2")
  library('tidyverse')
  library("rstatix")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")

figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

# Get and setup input files
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- read_tsv(indep_file)
tpm_file <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")
clin_file <- file.path(data_dir, "histologies-plot-group.tsv")
clin_tab  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         RNA_library == 'stranded',
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

hgg_bs_id <- clin_tab %>%
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

sf_list <- c("CLK1", "CLK2", "CLK3", "CLK4")

# Read TPM data
count_data <- readRDS(tpm_file) %>% 
  dplyr::select(any_of(hgg_bs_id))

# Add gene names as a column to count_data
count_data <- count_data %>%
  mutate(gene = rownames(.))

# Filter for CLK family genes and samples with expression
count_data_sf <- count_data %>%
  dplyr::filter(gene %in% sf_list) %>%
  dplyr::select(gene, any_of(clin_tab$Kids_First_Biospecimen_ID)) %>%
  dplyr::rowwise() %>%
  dplyr::filter(sum(c_across(where(is.numeric))) >= 1) %>%
  dplyr::ungroup()

# Reshape data to long format
data_long <- count_data_sf %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

# Check for missing values in expression and remove rows with NAs
data_long_clean <- data_long %>%
  filter(!is.na(expression))

# Perform pairwise comparisons using Wilcoxon test for each gene
stat.test <- data_long_clean %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sample) %>%
  adjust_pvalue(method = "BH") %>%
  ungroup()

clk_comparisons <- list(c("CLK1", "CLK2"), c("CLK1", "CLK3"), c("CLK1", "CLK4"))


# Plot with pairwise comparison results and mean labels
boxplot_tpm <- ggplot(data_long_clean, aes(x = gene, y = expression)) +
  geom_boxplot(outlier.shape = NA, fill = "grey", color = "#2c3e50") +
  stat_compare_means(comparisons = clk_comparisons, 
                     label = "p.value", 
                     method = "wilcox.test", 
                     position = "identity", 
                     label.x = 1) +  # Set label.x to the first gene (CLK1)+
  geom_jitter(width = 0.2, size = 2, shape = 21, color = "black", fill = "grey") + # Add actual data points
  labs(title = "CLK Family Expression", 
       x = "Gene", y = "TPM") +
  theme_Publication() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 30))) # Wrap x-axis labels 

# Save plot as PDF
pdf(file.path(plots_dir, "CLK-tpms.pdf"), 
    width = 6, height = 5)
print(boxplot_tpm)
dev.off()
