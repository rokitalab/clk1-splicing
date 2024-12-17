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
astro_trans_file <- "~/downloads/normal-brain-isoform-expression-rsem-tpm.rds"
astro_metadata_file<- "~/downloads/_1_SraRun_astrocytes_under40_hot.csv"

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
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

gtex_brain <- read_tsv(hist_file)  %>% 
  dplyr::filter(cohort == "GTEx",
                gtex_group == "Brain",
                Kids_First_Biospecimen_ID %in% gtex_rmats$sample_id) 

expr_tpm_tumor_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")

metadata_astrocytes = read_csv(astro_metadata_file) %>% 
  filter(cell_type=="Astrocyte") %>%
  select(source_name,Run) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=Run,
                plot_group=source_name) %>%
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
  mutate(group="Astrocytes")


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
  dplyr::mutate(group="Tumors")

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
  dplyr::mutate(group="Gtex",
                plot_group=str_replace_all(plot_group, "Brain - ", ""))
  
evo_devo_tpm <- readRDS("~/d3b_coding/neoepitope-identification/data/evodevo_rna-isoform-expression-rsem-tpm.rds") %>%
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


evodevo_histology_df <- vroom("~/d3b_coding/neoepitope-identification/data/evodevo-histologies.tsv") 
evodevo_clk1_transc_counts <- inner_join(evodevo_histology_df,evo_devo_tpm, by="Kids_First_Biospecimen_ID") %>%
  dplyr::filter(primary_site=='Hindbrain') %>%
  dplyr::mutate(plot_group=pathology_free_text_diagnosis,
                group="Evodevo") %>%
  dplyr::select(Kids_First_Biospecimen_ID,TPM,plot_group,group,transcript_id) 

## pediatric metadata
metadata_ped <- vroom("~/d3b_coding/neoepitope-identification/data/ped-normal-brain-histologies.tsv")

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
  dplyr::mutate(group="Pediatric") %>%
  dplyr::rename(Kids_First_Biospecimen_ID=Sample) %>%
  inner_join(sra_mapping,by=c("Kids_First_Biospecimen_ID" = "SRR"))

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
    unique()
)


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
    y = "Isoform Fraction") +
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

pdf(file.path(plots_dir,"clk4-tpm-phgg-ctrls.pdf"), height = 6, width = 22)
tpm_plot
dev.off()





# Data for the inset pie chart
inset_data <- transcript_expr_CLK1_combined_df %>%
  filter(proportion > 0.75) %>%
  count(plot_group, name = "count") %>%
  mutate(percentage = round(count / sum(count) * 100, 1))  # Calculate percentages

# Main scatter plot with points colored by `plot_group`
scatter_plot <- ggplot(transcript_expr_CLK1_combined_df, aes(x = log(total_TPM, 2), y = proportion)) +
  geom_point(
    aes(color = group),  # Color points based on `plot_group`
    size = 2                  # Set point size
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25)  # Y-axis ticks at intervals of 0.25
  ) +
  labs(
    title = "Exon 4 Transcript Expression",
    x = "Log2(Total TPM)",
    y = "Isoform Proportion"
  ) +
  theme_Publication()

# Inset pie chart with numeric annotations
inset_piechart <-  ggplot(inset_data, aes(x = "", y = count, fill = dot_color)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +  # Create a pie chart
  geom_text(
    aes(label = paste0(percentage, "%")),  # Add percentage labels inside the pie
    position = position_stack(vjust = 0.5),
    size = 4, color = "white"
  ) +
  labs(
    title = "Plot Groups > 0.75 Proportion",
    x = NULL,
    y = NULL
  ) +
  scale_fill_identity(
    name = "Dot Color"  # Retain original `dot_color` scheme
  ) +
  theme_void() +              # Minimalistic theme for pie chart
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(size = 10)  # Adjust font size for readability
  )

# Combine the main plot with the inset pie chart
combined_plot <- ggdraw() +
  draw_plot(scatter_plot, 0, 0, 1, 1) +  # Main plot takes full canvas
  draw_plot(inset_piechart, x = 0.05, y = 0.6, width = 0.35, height = 0.35)  # Inset in top-left corner

# Display the combined plot
print(combined_plot)