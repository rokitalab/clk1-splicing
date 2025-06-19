################################################################################
# 09-plot-clk1-ex4-transcripts-hgg-normals 
# Plot CLK1 exon 4 transcript expression with normals/ctrls
#
# written by Ammar S Naqvi
# Usage: Rscript 09-plot-clk1-ex4-transcripts-hgg-normals
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
input_dir <- file.path(analysis_dir, "input")

# Specify file paths
clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
hist_file  <- file.path(data_dir,"histologies.tsv")
rmats_file <- file.path(results_dir, "clk1-splice-events-rmats.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

## input files
gtex_trans_file <- file.path(data_dir,"gtex-harmonized-isoform-expression-rsem-tpm.rds")
ped_trans_file = file.path(data_dir,"GSE243682-normal-rna-isoform-expression-rsem-tpm.rds")
astro_trans_file <- file.path(data_dir,"GSE73721-normal-rna-isoform-expression-rsem-tpm.rds")
expr_tpm_tumor_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
gtex_rmats_file <- file.path(data_dir,"gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz")
expr_evodevo_file <-  file.path(data_dir,"evodevo_rna-isoform-expression-rsem-tpm.rds")
evodevo_hist_file <- file.path(data_dir,"evodevo-histologies.tsv")
astro_metadata_file<-file.path(data_dir,"GSE73721-normal-histologies.tsv")
metadata_pedr_file <-file.path(data_dir,"GSE243682-normal-brain-histologies.tsv")
pons_trans_file<- file.path(input_dir,"BA_KFWTGZPC_20250506.rsem.isoforms.results.gz")

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

## to get under 40 samples
gtex_rmats <- vroom(gtex_rmats_file) %>%
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
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

gtex_brain <- read_tsv(hist_file)  %>% 
  dplyr::filter(cohort == "GTEx",
                gtex_group == "Brain",
                Kids_First_Biospecimen_ID %in% gtex_rmats$sample_id) 


metadata_astrocytes = read_csv(astro_metadata_file) %>% 
  #filter(cell_type=="Astrocyte") %>%
  select(cell_type,Run) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=Run,
                plot_group=cell_type) %>%
    dplyr::mutate(plot_group=str_replace_all(plot_group, "Adult","Young")) %>%
  unique()

clk4_transcr_counts_astro <- readRDS(astro_trans_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",  # Rename specified transcripts
      TRUE ~ "Other"  # All other transcripts are renamed to "Other"
    )
  ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("SRR"), sum, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Kids_First_Biospecimen_ID",
    values_to = "TPM"
  ) %>%
  inner_join(metadata_astrocytes, by='Kids_First_Biospecimen_ID') %>%
  mutate(group="Cell Type Controls")


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
  dplyr::select(transcript_id,Kids_First_Biospecimen_ID,plot_group, TPM) %>%
  dplyr::mutate(group="PBTA")

gtex_clk1_transc_counts <- readRDS(gtex_trans_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",  # Rename specified transcripts
      TRUE ~ "Other"  # All other transcripts are renamed to "Other"
    )
  ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("GTEX"), sum, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Kids_First_Biospecimen_ID",
    values_to = "TPM"
  ) %>%
  inner_join(gtex_brain, by="Kids_First_Biospecimen_ID") %>%
  dplyr::rename(plot_group=gtex_subgroup) %>%
  
  dplyr::select(transcript_id,Kids_First_Biospecimen_ID,plot_group, TPM) %>%
  dplyr::mutate(group="GTEx",
                plot_group=str_replace_all(plot_group, "Brain - ", ""))
  
evo_devo_tpm <- readRDS(expr_evodevo_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%  # Filter for CLK1 genes
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",  # Rename specified transcripts
      TRUE ~ "Other" ) ) %>% # All other transcripts are renamed to "Other"  
  group_by(transcript_id) %>%
  summarise(across(starts_with("SAMEA"), sum, na.rm = TRUE)) %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  filter(Kids_First_Biospecimen_ID != "transcript_id") %>%  # Remove the "transcript_id" row
  pivot_longer(cols = starts_with("SAMEA"),  # Assuming your sample columns start with "SAMEA"
               names_to = "Sample_ID",
               values_to = "TPM") %>%
  select(transcript_id, Sample_ID, transcript_id, TPM) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=Sample_ID) # Select necessary columns


evodevo_histology_df <- vroom(evodevo_hist_file) 
evodevo_clk1_transc_counts <- inner_join(evodevo_histology_df,evo_devo_tpm, by="Kids_First_Biospecimen_ID") %>%
  dplyr::filter(primary_site=='Hindbrain') %>%
  dplyr::mutate(plot_group=pathology_free_text_diagnosis,
                group="Evo-Devo") %>%
  dplyr::select(Kids_First_Biospecimen_ID,TPM,plot_group,group,transcript_id) 

## pediatric metadata
metadata_ped <- vroom(metadata_pedr_file)

# Create a conversion table as a dataframe
sra_mapping <- data.frame(
  GSM = c("GSM7794191", "GSM7794190", "GSM7794185", "GSM7794180", "GSM7794176", "GSM7794214", "GSM7794155"),
  SRR = c("SRR26129063", "SRR26129064", "SRR26129066", "SRR26129078", "SRR26129084", "SRR26129123", "SRR26129178")
) %>% inner_join(metadata_ped,by=c("GSM" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(SRR, aliquot_id) %>%
  mutate(
    plot_group = case_when(
      grepl("CEREB", aliquot_id, ignore.case = TRUE) ~ "Cerebellum",
      grepl("PIT", aliquot_id, ignore.case = TRUE) ~ "Pituitary",
      grepl("PONS", aliquot_id, ignore.case = TRUE) ~ "Pons",
      grepl("frontal", aliquot_id, ignore.case = TRUE) ~ "Frontal_Cortex",
      TRUE ~ NA_character_  # Default to NA if no match
    )
  ) %>%
  dplyr::select(-aliquot_id) 
  
# ENST00000434813.3, ENST00000321356.9, ENST00000409403.6, ENST00000434813.3
ped_clk1_transcr_counts <- readRDS(ped_trans_file) %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  mutate(
    transcript_id = case_when(
      transcript_id %in% c("ENST00000321356.9","ENST00000434813.3", "ENST00000409403.6") ~ "Exon 4",  # Rename specified transcripts
      TRUE ~ "Other"  # All other transcripts are renamed to "Other"
    ) ) %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("SRR"), sum, na.rm = TRUE)) %>%
  #select(-gene_symbol) %>% # Remove 'gene_symbol'
  pivot_longer(
    cols = -transcript_id,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  dplyr::mutate(group="Pediatric normals") %>%
  dplyr::rename(Kids_First_Biospecimen_ID=Sample) %>%
  inner_join(sra_mapping,by=c("Kids_First_Biospecimen_ID" = "SRR"))

pons_counts <- vroom(pons_trans_file) %>%
  filter(grepl("_CLK1", transcript_id)) %>%
  mutate(transcript_id = ifelse(
    grepl("ENST00000321356.9|ENST00000434813.3|ENST00000409403.6", transcript_id),
    "Exon 4",  # Rename specified transcripts
    "Other"    # All other transcripts are renamed to "Other"
  )) %>%
  group_by(transcript_id) %>%
  dplyr::select(-gene_id) %>% # Remove 'gene_symbol'
  summarise(TPM = sum(TPM, na.rm = TRUE)) %>%
  dplyr::mutate(group="Pediatric normals",
                plot_group="Pons",
                Kids_First_Biospecimen_ID="BA_KFWTGZPC") %>%
  dplyr::select(Kids_First_Biospecimen_ID,transcript_id,TPM,group,plot_group)

ped_clk1_transcr_counts <- rbind(ped_clk1_transcr_counts, pons_counts) %>% 
  unique()

transcript_expr_CLK1_combined_df <- rbind(all_clk4_transcr_counts,ped_clk1_transcr_counts,gtex_clk1_transc_counts,evodevo_clk1_transc_counts,clk4_transcr_counts_astro) %>% 
  group_by(Kids_First_Biospecimen_ID) %>%
  mutate(total_TPM = sum(TPM[transcript_id %in% c("Exon 4", "Other")], na.rm = TRUE)) %>%
  mutate(proportion = ifelse(transcript_id == "Exon 4", TPM, 0) / total_TPM) %>%
  ungroup() %>% 
  filter(transcript_id=='Exon 4')#,
        # total_TPM > quantile(total_TPM, 0.25))

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




# Calculate the mean proportion for each 'plot_group' 
# Ensure plot_group is ordered as it appears in the dataframe
transcript_expr_CLK1_combined_df$plot_group <- factor(
  transcript_expr_CLK1_combined_df$plot_group,
  levels = unique(transcript_expr_CLK1_combined_df$plot_group)
)

# Calculate the mean of 'proportion' within each 'group' and 'plot_group'
transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  group_by(group, plot_group) %>%
  mutate(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
  ungroup()

# Reorder 'plot_group' within each 'group' based on 'mean_proportion' in descending order
transcript_expr_CLK1_combined_df$plot_group <- factor(
  transcript_expr_CLK1_combined_df$plot_group,
  levels = transcript_expr_CLK1_combined_df %>%
    arrange(group, desc(mean_proportion)) %>%
    pull(plot_group) %>%
    unique() )

# Define the desired order for groups
transcript_expr_CLK1_combined_df$group <- factor(
  transcript_expr_CLK1_combined_df$group, 
  levels = c("PBTA", "Pediatric normals", "Evo-Devo", "Cell Type Controls", "GTEx")
)

## order evodevo
evodevo_order <-  fct_relevel(
  "4 Week Post Conception",
  "5 Week Post Conception",
  "6 Week Post Conception",
  "7 Week Post Conception",
  "8 Week Post Conception",
  "9 Week Post Conception",
  "10 Week Post Conception",
  "11 Week Post Conception",
  "12 Week Post Conception",
  "13 Week Post Conception",
  "16 Week Post Conception",
  "Adolescent",
  "Elderly",
  "Infant",
  "Neonate",
  "School Age Child",
  "Toddler",
  "Young Adult",
  "Middle Adult")

# Reordering the plot group names for "Evo-Devo"
transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  mutate(plot_group = case_when(
    group == "Evo-Devo" ~ fct_relevel(plot_group,
                                      "4 Week Post Conception",
                                      "5 Week Post Conception",
                                      "6 Week Post Conception",
                                      "7 Week Post Conception",
                                      "8 Week Post Conception",
                                      "9 Week Post Conception",
                                      "10 Week Post Conception",
                                      "11 Week Post Conception",
                                      "12 Week Post Conception",
                                      "13 Week Post Conception",
                                      "16 Week Post Conception",
                                      "Adolescent",
                                      "Elderly",
                                      "Infant",
                                      "Neonate",
                                      "School Age Child",
                                      "Toddler",
                                      "Young Adult",
                                      "Middle Adult"
    ),
    TRUE ~ plot_group # Leave other groups as-is
  ))

## make plot for proportion
tpm_plot <- ggplot(transcript_expr_CLK1_combined_df, aes(x = plot_group, y = proportion)) +
  geom_jitter(
    aes(color = dot_color),   # Use precomputed colors for jitter points
    width = 0.2, size = 2) +
  geom_boxplot(
    aes(group = plot_group),  # Create boxplots for each group
    width = 0.6,              # Adjust the width of the boxplots
    color = "black",          # Set the color of the boxplot borders
    fill = "white",           # Fill color for the boxplots
    alpha = 0.2 ) +
  labs(
    title = "Relative CLK1 Exon 4 Transcript Expression",
    x = "Group",
    y = "Proportion CLK1 exon 4 inclusion transcript") +
  scale_color_identity(
    name = "Group"
  ) +
  facet_grid(
    ~ group,                  # Facet by 'group'
    scales = "free_x"         # Allow different x-axis scales for each facet
  )  +
  theme_Publication() +
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle = 75, hjust = 1)
  ) 

pdf(file.path(plots_dir,"clk1_ex4-tpm-phgg-ctrls.pdf"), height = 7, width = 22)
tpm_plot
dev.off()