################################################################################
# 08-plot-Ex4-PSI-all-tumors.R
# script that plots CLK1 exon 4 PSI variations
# written by Ammar Naqvi
#
# usage: Rscript 08-plot-Ex4-PSI-all-tumors.R
################################################################################

suppressPackageStartupMessages({
  library("reshape2")
  library("tidyverse")
  library("ggpubr")
  library("ggplot2")
  library("vroom")
  library("data.table")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")

## create plots dir if it doesn't exist
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## input files
rmats_file <- file.path(results_dir,"clk1-splice-events-rmats.tsv")
clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")
expr_file <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")
gtex_rmats <- file.path(data_dir,"gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz")
pedr_rmats <- file.path(data_dir,"GSE243682-normal-splice-events-rmats.tsv.gz")

## output files for final plots
hgg_plot_file <- file.path(plots_dir,"all_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")
dmg_plot_file <- file.path(plots_dir,"dmg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")
other_hgg_plot_file <- file.path(plots_dir,"other_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

indep_df <- vroom(indep_file)

# Function to add TPM values to a dataframe from a matrix
add_TPM_values <- function(df, matrix) {
  # Ensure gene_symbol and sample_id columns exist in df
  if (!("gene_symbol" %in% names(df)) | !("Kids_First_Biospecimen_ID" %in% names(df))) {
    stop("DataFrame must contain 'gene_symbol' and 'sample_id' as columns")
  }
  
  # Ensure rownames and colnames are properly set in the matrix
  if (is.null(rownames(matrix)) | is.null(colnames(matrix))) {
    stop("Matrix must have rownames and colnames set for gene_symbol and Sample_id, respectively")
  }
  
  
  # Extract the matching TPM values
  df$gene_tpm <- mapply(function(gene, sample) {
    if (gene %in% rownames(matrix) & sample %in% colnames(matrix)) {
      matrix[gene, sample]
    } else {
      NA
    }
  }, df$gene_symbol, df$Kids_First_Biospecimen_ID)
  
  return(df)
}

## load histologies info for HGG subty  
histologies_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         RNA_library=='stranded',
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

## load rmats input for CLK1
clk1_rmats <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4, HGGs
  dplyr::filter(
    geneSymbol=="CLK1",
    exonStart_0base=="200860124", 
    exonEnd=="200860215",
  ) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=sample_id) %>% 
  dplyr::select(Kids_First_Biospecimen_ID,IncLevel1) %>%
  # Join rmats data with clinical data
  inner_join(histologies_df, by='Kids_First_Biospecimen_ID') %>%
  mutate(gene_symbol="CLK1")

exp <- readRDS(expr_file) %>%
  dplyr::select(any_of(histologies_df$Kids_First_Biospecimen_ID))

var_exp_filt <- add_TPM_values(clk1_rmats, exp) 

# filter for very high expression 
high_expression <- quantile(var_exp_filt$gene_tpm, 0.90)
ex4_psi_filtered <- var_exp_filt %>%
   filter(gene_tpm > high_expression) %>%
  # Compute variance and add as a new column
  group_by(plot_group) %>%
  mutate(PSI_variance = sd(IncLevel1, na.rm = TRUE)) %>%
  filter(!is.na(PSI_variance))

# Plot with pairwise comparison results and mean labels
boxplot_tpm<- ggplot(var_exp_filt, aes(x = plot_group, y = IncLevel1)) +
  geom_boxplot(outlier.shape = NA, fill = "grey", color = "#2c3e50") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  geom_jitter(width = 0.2, size = 2, shape = 21, color = "black", fill="grey") + # Add actual data points
  labs(title = "CLK Exon 4 PSI", 
       x = "Histology", y = "PSI") +
  theme_Publication() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 30))) # Wrap x-axis labels 


var_plot<- ggplot(data=ex4_psi_filtered, 
       aes(reorder(plot_group, PSI_variance),PSI_variance,  
           group=1), color="black") +
  geom_point(aes(fill = 'black'), size = 3, pch = 21, color="black") +  # Map color inside aes()
  xlab("Histology") + 
  ylab("Standard Deviation") + 
  ggtitle("CLK1 Exon 4 variation") +
  theme_Publication() + 
  theme(axis.text.x=element_text(angle = 75, hjust = 1, size = 11),legend.position = "none") 

ex4_psi_range <- ex4_psi_filtered %>%
  group_by(plot_group) %>%
  mutate(PSI_range = max(IncLevel1, na.rm = TRUE) - min(IncLevel1, na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::select(plot_group,PSI_range) %>%
  unique()

psi_range_plot<- ggplot(data=ex4_psi_range, 
       aes(reorder(plot_group, PSI_range),PSI_range,  
           group=1), color="black") +
  geom_point(aes(fill = 'black'), size = 3, pch = 21, color="black") +  # Map color inside aes()
  xlab("Histology") + 
  ylab("PSI Range") + 
  ggtitle("CLK1 Exon 4 PSI Range") +
  theme_Publication() + 
  theme(axis.text.x=element_text(angle = 75, hjust = 1, size = 11),legend.position = "none") 


# Save plot as PDF
pdf(file.path(plots_dir, "CLK1-Ex4-sdev-across.pdf"), 
    width = 4, height = 6)
print(var_plot)
dev.off()

pdf(file.path(plots_dir, "CLK1-Ex4-range-across.pdf"), 
    width = 4, height = 6)
print(psi_range_plot)
dev.off()

## add ctrls to plots
gtex_psi_df <- vroom(gtex_rmats) %>% # Select CLK1 gene
  dplyr::filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  dplyr::select(sample_id, geneSymbol, IncLevel1) %>%
  dplyr::rename(gene_symbol=geneSymbol) %>% 
  dplyr::mutate(plot_group="gtex-brain")

ped_psi_df <- vroom(pedr_rmats) %>% # Select CLK1 gene
  dplyr::filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  dplyr::select(sample_id, geneSymbol, IncLevel1) %>%
  dplyr::rename(gene_symbol=geneSymbol) %>% 
  dplyr::mutate(plot_group="pediatric_ctrls")

ctrl_psi_df<- rbind(ped_psi_df,gtex_psi_df)

ctrl_psi_range <- ctrl_psi_df %>% group_by(plot_group) %>%
  mutate(PSI_range = max(IncLevel1, na.rm = TRUE) - min(IncLevel1, na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::select(plot_group,PSI_range) %>%
  unique()

combo_psi_range <- rbind(ctrl_psi_range,ex4_psi_range)

psi_range_plot<- ggplot(data=combo_psi_range, 
                        aes(reorder(plot_group, PSI_range),PSI_range,  
                            group=1), color="black") +
  geom_point(aes(fill = 'black'), size = 3, pch = 21, color="black") +  # Map color inside aes()
  xlab("Histology") + 
  ylab("PSI Range") + 
  ggtitle("CLK1 Exon 4 PSI Range") +
  theme_Publication() + 
  theme(axis.text.x=element_text(angle = 75, hjust = 1, size = 11),legend.position = "none") 

pdf(file.path(plots_dir, "CLK1-Ex4-range-across-ctrls.pdf"), 
    width = 4, height = 6)
print(psi_range_plot)
dev.off()
