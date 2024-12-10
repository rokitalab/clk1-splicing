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
input_dir <- file.path(root_dir, "results")

# Specify file paths
clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
hist_file  <- file.path(data_dir,"histologies.tsv")
rmats_file <- file.path(results_dir, "clk1-splice-events-rmats.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

gtex_trans_file <- file.path("/Users/naqvia/d3b_coding/neoepitope-identification/data/gtex-harmonized-isoform-expression-rsem-tpm.rds")
ped_trans_file = "~/d3b_coding/neoepitope-identification/data/GSE243682_normal_rna-isoform-expression-rsem-tpm.rds"

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

## to get under 40 samples
gtex_rmats <- vroom("/Users/naqvia/d3b_coding/neoepitope-identification/data/gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz") %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample_id, geneSymbol, IncLevel1) %>%
  dplyr::rename(gene_symbol=geneSymbol)

indep_df <- read_tsv(indep_file)

hist_indep_rna_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         #grepl("poly", RNA_library),
         #RNA_library=='stranded',
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

all_clk4_transcr_counts <- readRDS(expr_tpm_tumor_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",  # Rename specified transcripts
      TRUE ~ "Other"  # All other transcripts are renamed to "Other"
    )
  ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("BS"), sum, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Kids_First_Biospecimen_ID",
    values_to = "TPM"
  ) %>%
  inner_join(hist_indep_rna_df, by="Kids_First_Biospecimen_ID") %>%
  #dplyr::select(transcript_id,Kids_First_Biospecimen_ID,plot_group, TPM) %>%
  dplyr::mutate(group="Tumors")

all_clk4_transcr_years_counts <- all_clk4_transcr_counts %>%  select(transcript_id,Kids_First_Biospecimen_ID,plot_group, TPM, age_at_diagnosis_days) %>% mutate(years=age_at_diagnosis_days/365) %>% select(transcript_id,Kids_First_Biospecimen_ID,plot_group, TPM, years) %>% 
  mutate(
    year_group = cut(
      years,
      breaks = c(0, 1, 2, 6, 12, 18, 21, 39),
      labels = c("[0-1)", "[1-2)", "[2-6)", "[6-12)", "[12-18)", "[18-21)", "[21-39)"),
      right = FALSE
    )) %>% 
  na.omit()



transcript_expr_CLK1_combined_df <- all_clk4_transcr_years_counts %>% 
  group_by(Kids_First_Biospecimen_ID) %>%
  mutate(total_TPM = sum(TPM[transcript_id %in% c("Exon 4", "Other")], na.rm = TRUE)) %>%
  mutate(proportion = ifelse(transcript_id == "Exon 4", TPM, 0) / total_TPM) %>%
  ungroup() %>% 
  filter(transcript_id=='Exon 4') 

color_df <- hist_indep_rna_df %>%
  dplyr::select(plot_group_hex, plot_group) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  unique()

cols <- as.character(color_df$plot_group_hex)
names(cols) <- as.character(color_df$plot_group)

# Ensure no duplicate rows
color_df <- color_df %>%
  dplyr::distinct()

# Create the color mapping
color_mapping <- setNames(color_df$plot_group_hex, color_df$plot_group)

# Add a color column to the data based only on plot_group
transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  mutate(
    dot_color = case_when(
      plot_group %in% names(color_mapping) ~ color_mapping[plot_group], # Use plot_group-specific colors
      TRUE ~ "#b5b5b5"                                                  # Default to grey
    )
  ) %>%
  filter(!is.na(plot_group)) # Remove rows with NA in plot_group

# Ensure plot_group is ordered as it appears in the dataframe
transcript_expr_CLK1_combined_df$plot_group <- factor(
  transcript_expr_CLK1_combined_df$plot_group,
  levels = unique(transcript_expr_CLK1_combined_df$plot_group)
)

# Calculate the mean of 'proportion' within each 'group' and 'plot_group'
transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  group_by(plot_group) %>%
  mutate(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
  ungroup()

# Reorder 'plot_group' within each 'group' based on 'mean_proportion' in descending order
transcript_expr_CLK1_combined_df$plot_group <- factor(
  transcript_expr_CLK1_combined_df$plot_group,
  levels = transcript_expr_CLK1_combined_df %>%
    arrange(plot_group, desc(mean_proportion)) %>%
    pull(plot_group) %>%
    unique()
)


## make plot for proportion
tpm_plot <-ggplot(transcript_expr_CLK1_combined_df, aes(x = year_group, y = proportion)) +
  geom_jitter(
    aes(color = dot_color),   # Use precomputed colors for jitter points
    width = 0.2, size = 2
  ) +
  geom_boxplot(
    aes(group = year_group),  # Create boxplots for each group
    width = 0.6,              # Adjust the width of the boxplots
    color = "black",          # Set the color of the boxplot borders
    fill = "white",           # Fill color for the boxplots
    alpha = 0.2
  ) +
  labs(
    title = "ENST00000321356 Expression (TPM)",
    x = "Age Group",
    y = "Isoform Fraction"
  ) +
  scale_color_identity(name = "Dot Color") +
  facet_wrap(~ plot_group, scales = "free_x", nrow = 1) +  # Facet with a single row
  theme_Publication() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 75, hjust = 1)
  )


pdf(file.path(plots_dir,"clk4-tpm-phgg-ctrls.age_years-v2.pdf"), height = 6, width = 26)
tpm_plot
dev.off()


