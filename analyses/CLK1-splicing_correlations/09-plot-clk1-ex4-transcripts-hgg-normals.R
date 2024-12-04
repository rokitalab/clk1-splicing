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

gtex_trans_file <- file.path(data_dir,"GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")

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
hist_rna_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

gtex_brain <- read_tsv(hist_file)  %>% 
  dplyr::filter(cohort == "GTEx",
                gtex_group == "Brain") 

hgg_bs_id <- hist_rna_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) 

expr_tpm_tumor_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")

# transcript is ENST00000321356.9 in pbta
hgg_clk4_transcr_counts <- readRDS(expr_tpm_tumor_file) %>%
  filter(transcript_id == "ENST00000321356.9") %>%
  select(-gene_symbol) %>% # Remove 'gene_symbol'
  pivot_longer(
    cols = -transcript_id,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(group="pHGGs")

ped_trans_file = "~/d3b_coding/neoepitope-identification/data/GSE243682_normal_rna-isoform-expression-rsem-tpm.rds"
ped_clk1_transcr_counts <- readRDS(ped_trans_file) %>%
  filter(transcript_id == "ENST00000321356.9") %>%
  select(-gene_symbol) %>% # Remove 'gene_symbol'
  pivot_longer(
    cols = -transcript_id,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(group="pediatric_ctrls")

gtex_clk1_transcr_counts <- fread(gtex_trans_file, skip = 2) %>%
  select(c(transcript_id, gtex_brain$Kids_First_Biospecimen_ID)) %>%
  filter(transcript_id == "ENST00000321356.8") %>%
  select(-transcript_id) %>%
  {data.frame(t(.))} %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(ENST00000321356 = `t...`) %>%
  left_join(gtex_brain[,c("Kids_First_Biospecimen_ID", "gtex_subgroup")]) %>% 
  rename(group=gtex_subgroup,
         Sample=Kids_First_Biospecimen_ID,
         TPM=ENST00000321356) %>%
  mutate(transcript_id='ENST00000321356.9')


transcript_expr_CLK1_combined_df <- rbind(hgg_clk4_transcr_counts,ped_clk1_transcr_counts,gtex_clk1_transcr_counts)

## add pediatric metadata
metadata_ped <- vroom("~/d3b_coding/neoepitope-identification/data/ped-normal-brain-histologies.tsv")

# Create a conversion table as a dataframe
sra_mapping <- data.frame(
  GSM = c("GSM7794191", "GSM7794190", "GSM7794185", "GSM7794180", "GSM7794176", "GSM7794214", "GSM7794155"),
  SRR = c("SRR26129063", "SRR26129064", "SRR26129066", "SRR26129078", "SRR26129084", "SRR26129123", "SRR26129178")
) %>% inner_join(metadata_ped,by=c("GSM" = "Kids_First_Biospecimen_ID")) %>%
  inner_join(ped_clk1_transcr_counts,by=c("SRR" = "Sample")) %>%
  dplyr::select(SRR, aliquot_id) %>%
  mutate(
    region = case_when(
      grepl("CEREB", aliquot_id, ignore.case = TRUE) ~ "Cerebellum",
      grepl("PIT", aliquot_id, ignore.case = TRUE) ~ "Pituitary",
      grepl("PONS", aliquot_id, ignore.case = TRUE) ~ "Pons",
      grepl("frontal", aliquot_id, ignore.case = TRUE) ~ "Frontal_Cortex",
      TRUE ~ NA_character_  # Default to NA if no match
    )
  )
  

# View the conversion table
print(sra_mapping)

# Calculate mean and 2 standard deviations of pediatric_ctrls
ctrl_stats <- transcript_expr_CLK1_combined_df %>%
  filter(group == "pediatric_ctrls") %>%
  summarise(
    mean_ctrl = mean(TPM, na.rm = TRUE),
    sd_ctrl = sd(TPM, na.rm = TRUE)
  )

mean_ctrl <- ctrl_stats$mean_ctrl
sd_ctrl <- ctrl_stats$sd_ctrl

# Define region-specific colors
region_colors <- c(
  "Cerebellum" = "#E69F00",  # Orange
  "Pituitary" = "#56B4E9",  # Sky blue
  "Pons" = "#009E73",       # Green
  "Frontal_Cortex" = "#F0E442",    # Yellow
  "pHGGs" = "#0072B2",      # Blue
  "pediatric_ctrls" = "#D55E00"  # Vermillion
)

# Add region column to the combined dataframe for pediatric_ctrls
transcript_expr_CLK1_combined_df <- transcript_expr_CLK1_combined_df %>%
  left_join(sra_mapping %>% select(SRR, region), by = c("Sample" = "SRR"))


# Modify the plot
tpm_plot <- ggplot(transcript_expr_CLK1_combined_df, aes(x = group, y = TPM)) +
  geom_jitter(
    aes(color = ifelse(group == "pediatric_ctrls", region, group)),  # Use region for pediatric_ctrls
    width = 0.2, size = 2, alpha = 0.7
  ) +
  geom_hline(yintercept = mean_ctrl + 2 * sd_ctrl, linetype = "dashed", color = "black", size = 1) +  # +2 SD line
  labs(
    title = "ENST00000321356 Expression",
    x = "Group",
    y = "TPM"
  ) +
  scale_color_manual(
    values = region_colors,
    name = "Pediatric Samples"  # Set legend title
  ) +
  theme_Publication() +
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle = 75, hjust = 1)
  ) +
  scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 30))) + 
  scale_y_continuous(breaks = seq(0, max(transcript_expr_CLK1_combined_df$TPM), by = 50))  # Ticks every 50


pdf(file.path(plots_dir,"clk4-tpm-phgg-ctrls.pdf"), height = 8, width = 9)
tpm_plot
dev.off()


