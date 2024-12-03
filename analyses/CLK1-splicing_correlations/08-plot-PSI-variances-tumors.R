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
rmats_file <- file.path(data_dir,"splice-events-rmats.tsv.gz")
clk1_rmats_file <- file.path(results_dir,"splice-events-rmats.tsv")

clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")
expr_file <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")

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
rmats_df <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4, HGGs
  dplyr::filter(
    geneSymbol %in% c("CLK1", "CLK2", "CLK3", "CLK4", "SRPK1"),
    splicing_case=='SE'
  ) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=sample_id) %>% 
  dplyr::select(Kids_First_Biospecimen_ID,geneSymbol, exonStart_0base,exonEnd,IncLevel1) %>%
  # Join rmats data with clinical data
  inner_join(histologies_df, by='Kids_First_Biospecimen_ID') 

## create a splice_id 
rmats_modified <- rmats_df %>% mutate(SpliceID = paste0(geneSymbol, ":", exonStart_0base, "-", exonEnd))

## splice_ids of splice events from kinases of interest
filter_splice_ids <- c(
  "CLK1:200860124-200860215",
  "CLK2:155265859-155265954",
  "CLK2:155268707-155268792",
  "CLK2:155268707-155268795",
  "CLK3:74622119-74622216",
  "CLK3:74625801-74625968",
  "CLK4:178617343-178617434",
  "SRPK1:35870280-35870494"
)

rmats_modified <- rmats_modified %>%
  filter(SpliceID %in% filter_splice_ids) %>% 
  as_tibble() %>% 
  select(geneSymbol,SpliceID,IncLevel1, Kids_First_Biospecimen_ID,plot_group) %>%
  dplyr::rename(gene_symbol=geneSymbol)


exp <- readRDS(expr_file) %>%
  select(any_of(histologies_df$Kids_First_Biospecimen_ID))

var_exp_filt <- add_TPM_values(rmats_modified, exp) 

# filter for very high expression 
# Step 1: Compute the 90th percentile of gene_tpm for each gene_symbol
gene_tpm_90th <- var_exp_filt %>%
  group_by(gene_symbol) %>%
  summarize(high_expression_threshold = quantile(gene_tpm, 0.90, na.rm = TRUE), .groups = "drop")

# Step 2: Join the thresholds back to the original dataframe
var_exp_filt <- var_exp_filt %>%
  left_join(gene_tpm_90th, by = "gene_symbol")

# Step 3: Filter rows where gene_tpm exceeds the 90th percentile for the corresponding gene_symbol
ex4_psi_filtered <- var_exp_filt %>%
  filter(gene_tpm > high_expression_threshold) %>%  # Ensure column exists
  group_by(SpliceID, plot_group) %>%  # Group by both SpliceID and plot_group
  mutate(PSI_variance = sd(IncLevel1, na.rm = TRUE)) %>%  # Compute variance
  filter(!is.na(PSI_variance)) %>%  # Filter out rows where variance is NA
  ungroup()



var_plot<- ggplot(data = ex4_psi_filtered, 
                  aes(x = reorder(plot_group, PSI_variance), y = PSI_variance, group = 1)) +
  geom_point(size = 3, pch = 21, color = "black", fill = "black") +
  xlab("Histology") + 
  ylab("Standard Deviation") + 
  ggtitle("PSI Variation by SpliceID") +
  theme_Publication() + 
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 11),
    legend.position = "none"
  ) +
  facet_wrap(~SpliceID, scales = "free_y", ncol = 2)  # Facet by SpliceID, adjust columns


ex4_psi_range <- ex4_psi_filtered %>%
  group_by(SpliceID, plot_group) %>%  # Group by both SpliceID and plot_group
  mutate(PSI_range = max(IncLevel1, na.rm = TRUE) - min(IncLevel1, na.rm = TRUE)) %>%
  ungroup() %>% 
  select(plot_group,SpliceID,PSI_range) %>%
  unique()

psi_range_plot<- ggplot(data=ex4_psi_range, 
                        aes(reorder(plot_group, PSI_range),PSI_range,  
                            group=1), color="black") +
  geom_point(aes(fill = 'black'), size = 3, pch = 21, color="black") +  # Map color inside aes()
  xlab("Histology") + 
  ylab("PSI Range") + 
  ggtitle("CLK1 Exon 4 PSI Range") +
  theme_Publication() + 
  theme(axis.text.x=element_text(angle = 75, hjust = 1, size = 11),legend.position = "none") +
  facet_wrap(~SpliceID, scales = "free_y", ncol = 2)  # Facet by SpliceID, adjust columns



# Save plot as PDF
pdf(file.path(plots_dir, "CLK1-Ex4-sdev-across.pdf"), 
    width = 4, height = 6)
print(var_plot)
dev.off()