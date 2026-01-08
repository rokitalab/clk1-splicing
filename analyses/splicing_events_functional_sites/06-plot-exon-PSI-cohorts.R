################################################################################
# 06-plot-exon-PSI-cohorts.R
# Plot exon transcript expression with normals/ctrls
#
# written by Ammar S Naqvi, Patricia Sullivan
# Usage: Rscript 06-plot-exon-PSI-cohorts.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(vroom)
  library(data.table)
  library(ggtext)
  library(ggpubr)
  library(rstatix)
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")

# Specify file paths
clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
hist_file  <- file.path(data_dir,"histologies.tsv")
rmats_file <- file.path(results_dir, "clk1-splice-events-rmats.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

## input files
gtex_trans_file <- file.path(data_dir,"gtex-harmonized-isoform-expression-rsem-tpm.rds")
expr_tpm_tumor_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
gtex_rmats_file <- file.path(data_dir,"gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz")
expr_evodevo_file <-  file.path(data_dir,"evodevo_rna-isoform-expression-rsem-tpm.rds")
evodevo_hist_file <- file.path(data_dir,"evodevo-histologies.tsv")
gtex_hist_file <- file.path(input_dir,"gtex-samples-by-age.tsv")
cluster_file <- file.path(root_dir, "analyses",
			"sample-psi-clustering", "results",
			"sample-cluster-metadata-top-5000-events-stranded.tsv")

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(root_dir, "figures/theme_for_plots.R"))

get_normals_psi <- function(gene, exon_no, transcript_list) {

  all_bs_id <- read_tsv(cluster_file) %>%
  	pull(sample_id)
  
  indep_df <- read_tsv(indep_file)
  hist_indep_rna_df  <-  read_tsv(hist_file) %>%
    filter(cohort == "PBTA",
  	       Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)
  
  gtex_brain <- read_tsv(gtex_hist_file)  %>% 
    dplyr::filter(gtex_group == "Brain")
  
  
  pbta_transcr_counts <- readRDS(expr_tpm_tumor_file) %>%
    filter(sub("-[^-]+$", "", gene_symbol) == gene) %>%
    mutate(
      transcript_id = case_when(
        transcript_id %in% transcript_list ~ "Target",  # Rename specified transcripts
        TRUE ~ "Other"  # All other transcripts are renamed to "Other"
      )
    ) %>%
    group_by(transcript_id) %>%
    summarise(across(starts_with("BS"), sum, na.rm = TRUE)) %>%
    pivot_longer(
      cols = -transcript_id,
      names_to = "Kids_First_Biospecimen_ID",
      values_to = "TPM"
    )
  
  pbta_stranded_transcr_counts <- pbta_transcr_counts %>%
    inner_join(hist_indep_rna_df, by="Kids_First_Biospecimen_ID") %>%
    filter(RNA_library == "stranded") %>%
    mutate(
      age_years = age_at_diagnosis_days / 365,
      group = "PBTA stranded",
      plot_group = cut(
        age_years,
        breaks = c(0, 14, 18, 39),
        labels = c("0-14", "15-18", "19-39"),
        right = TRUE
      ),
    ) %>%
    select(Kids_First_Biospecimen_ID, group, plot_group, transcript_id, TPM) %>%
    na.omit()
  
  pbta_polyA_transcr_counts <- pbta_transcr_counts %>%
    inner_join(hist_indep_rna_df, by="Kids_First_Biospecimen_ID") %>%
    filter(RNA_library %in% c("poly-A stranded", "poly-A")) %>%
    mutate(
      age_years = age_at_diagnosis_days / 365,
      group = "PBTA polyA",
      plot_group = cut(
        age_years,
        breaks = c(0, 14, 18, 39),
        labels = c("0-14", "15-18", "19-39"),
        right = TRUE
      ),
    ) %>%
    select(Kids_First_Biospecimen_ID, group, plot_group, transcript_id, TPM) %>%
    na.omit()
  
  
  gtex_transc_counts <- readRDS(gtex_trans_file) %>%
    filter(sub("-[^-]+$", "", gene_symbol) == gene) %>%
    mutate(
      transcript_id = case_when(
        transcript_id %in% transcript_list ~ "Target",  # Rename specified transcripts
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
    dplyr::mutate(group="GTEx",
                  plot_group=AGE) %>%
    dplyr::select(transcript_id,Kids_First_Biospecimen_ID, TPM, group, plot_group)
  
    
    
  evo_devo_tpm <- readRDS(expr_evodevo_file) %>%
    filter(sub("-[^-]+$", "", gene_symbol) == gene) %>%  # Filter for gene
    mutate(
      transcript_id = case_when(
        transcript_id %in% transcript_list ~ "Target",  # Rename specified transcripts
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
  
  ### Postnatally we sampled neonates, “infants” (6-9 months), “toddlers” (2-4 years), “school” (7-9 years), “teenagers” (13-19 years), and then adults from each decade until 63 years of age (Supplementary Tables 1–2).... ‘ya’ young adult (25-32 years) and ‘sen’ senior (58-63 years)....“yma” young middle age (39-41 years)
  ### https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/#S7
  
  evodevo_histology_df <- vroom(evodevo_hist_file) 
  evodevo_transc_counts <- inner_join(evodevo_histology_df,evo_devo_tpm, by="Kids_First_Biospecimen_ID") %>%
    dplyr::filter(primary_site=='Hindbrain') %>%
    dplyr::mutate(group="Evo-Devo",
                  plot_group = case_when(
                    pathology_free_text_diagnosis %in% c(
                      "4 Week Post Conception", "5 Week Post Conception", "6 Week Post Conception",
                      "7 Week Post Conception", "8 Week Post Conception", "9 Week Post Conception",
                      "10 Week Post Conception", "11 Week Post Conception", "12 Week Post Conception",
                      "13 Week Post Conception", "16 Week Post Conception") ~ "Fetal",
                pathology_free_text_diagnosis %in% c("Neonate", "Infant", "Toddler") ~ "0-4",
                pathology_free_text_diagnosis %in% c("School Age Child", "Adolescent") ~ "7-19",
                pathology_free_text_diagnosis == "Young Adult" ~ "25-32",
                pathology_free_text_diagnosis %in% c("Middle Adult", "Elderly") ~ "39-63",
                TRUE ~ NA)) %>%
    dplyr::select(Kids_First_Biospecimen_ID,TPM,plot_group,group,transcript_id) 
  
  transcript_expr_combined_df <- rbind(pbta_stranded_transcr_counts,pbta_polyA_transcr_counts,gtex_transc_counts,evodevo_transc_counts) %>% 
    group_by(Kids_First_Biospecimen_ID) %>%
    mutate(total_TPM = sum(TPM[transcript_id %in% c("Target", "Other")], na.rm = TRUE)) %>%
    mutate(proportion = ifelse(transcript_id == "Target", TPM, 0) / total_TPM) %>%
    ungroup() %>% 
    filter(transcript_id=='Target')#,
          # total_TPM > quantile(total_TPM, 0.25))
  
  # Define the desired order for groups
  transcript_expr_combined_df$plot_group <- factor(
    transcript_expr_combined_df$plot_group,
   # levels = c("0-14", "15-18", "19-39", "Fetal", "Early Childhood", "School Age", "Adult", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
    levels = c("0-14", "15-18", "19-39", "Fetal", "0-4", "7-19", "25-32", "39-63", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
  )
  
  color_pal <- c(
    "PBTA stranded" = "#FDBF6F",
    "PBTA polyA" = "coral",
    "Evo-Devo" = "mediumseagreen",
    "GTEx" = "#6ca6da")
  
  
  # Get stats across cohorts
  stat.tests.cohort <- transcript_expr_combined_df %>%
    wilcox_test(proportion ~ group, p.adjust.method = "BH")
  
  # Get stats for age bins across cohorts
  stat.tests.all <- transcript_expr_combined_df %>%
    mutate(full_group = paste0(group, ": ", plot_group)) %>%
    wilcox_test(proportion ~ full_group, p.adjust.method = "BH")
  
  # Save both
  bind_rows(stat.tests.cohort, stat.tests.all) %>%
    select(-`.y.`) %>%
    write_tsv(file = file.path(results_dir, paste0(gene, "-exon", exon_no, "-psi-normals-stats.tsv")))
  
  # Get stats within groups for plot
  stat.tests <- transcript_expr_combined_df %>%
    group_by(group) %>%
    wilcox_test(proportion ~ plot_group, p.adjust.method = "BH") %>%
    filter(p.adj < 0.05)
  
  # get position by each facet
  stat.tests <- stat.tests %>%
    group_by(group) %>%
    mutate(
      x    = as.numeric(factor(group1, levels = unique(transcript_expr_combined_df$plot_group))),
      xend = as.numeric(factor(group2, levels = unique(transcript_expr_combined_df$plot_group))),
      y_base   = max(transcript_expr_combined_df$proportion[transcript_expr_combined_df$group == first(group)]),
      
      # compute y-adjustment for overlapping comparisons within a facet
      step_h   = y_base * 0.05,
      comp_id  = row_number(),                           
      y.position = y_base * 1.05 + (comp_id - 1) * step_h
    ) %>%
    ungroup() %>%                                       
    select(-y_base, -step_h, -comp_id)
  
  ## make plot for proportion
  tpm_plot <- ggplot(transcript_expr_combined_df, aes(x = plot_group, y = proportion)) +
    geom_jitter(
      aes(color = group),   # Use precomputed colors for jitter points
      width = 0.2, size = 1.5) +
    geom_boxplot(
      aes(group = plot_group),  # Create boxplots for each group
      width = 0.6,              # Adjust the width of the boxplots
      color = "black",          # Set the color of the boxplot borders
      fill = "white",           # Fill color for the boxplots
      alpha = 0.2,
      outlier.shape = NA) +
    labs(
      x = "Age",
      y = paste0("Proportion <i>", gene, "</i> exon ", exon_no, "<br>inclusion in transcript")
     ) +
    theme_Publication() +
    theme(
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown()
    ) +
    scale_color_manual(values = color_pal) +
    facet_wrap(~group, scales = "free_x", nrow = 1) +
    stat_pvalue_manual(stat.tests,
                       tip.length   = 0.01) + 
    scale_y_continuous(
                         limits = c(0.0, 1.0),
                         breaks = seq(0.0, 1.0, by = 0.25)  # Adjust the step as needed
                       )
  
  pdf(file.path(plots_dir, paste0(gene, "-exon", exon_no, "-tpm-ctrls-summary.pdf")), height = 3.5, width = 9)
    print(tpm_plot)
  dev.off()
}

# Gene symbol, exon number, transcripts exon is included in
get_normals_psi("CLK1", "4", c("ENST00000321356.9", "ENST00000434813.3", "ENST00000409403.6"))
get_normals_psi("PRKDC", "80", c("ENST00000314191.7", "ENST00000521331.5"))

