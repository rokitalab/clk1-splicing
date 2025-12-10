################################################################################
# 07-plot-sbi-high-low.R
# script that takes in SBI tsv files and plots levels by histology
#
# written by Ammar Naqvi, Jo Lynne Rokita, Patricia Sullivan
#
# usage: Rscript 07-plot-sbi-high-low.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("ggforce")
  library("dplyr")
  library("tidyverse")
  library("viridis")
  library("RColorBrewer")
  })

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
map_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
palette_dir <- file.path(root_dir, "palettes")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

# read in SI files
splice_index_SE_file   <- file.path(results_dir, "splicing_index.SE.txt")
splice_index_RI_file   <- file.path(results_dir, "splicing_index.RI.txt")
splice_index_A5SS_file <- file.path(results_dir, "splicing_index.A5SS.txt")
splice_index_A3SS_file <- file.path(results_dir, "splicing_index.A3SS.txt")

splice_index_SE_df  <-  read_tsv(splice_index_SE_file) 
splice_index_RI_df  <-  read_tsv(splice_index_RI_file) 
splice_index_A5SS_df  <-  read_tsv(splice_index_A5SS_file) 
splice_index_A3SS_df  <-  read_tsv(splice_index_A3SS_file) 

# put them in a vector
si_files <- c(splice_index_SE_file,
              splice_index_RI_file,
              splice_index_A5SS_file,
              splice_index_A3SS_file)

# read and bind together
si_df <- si_files %>%
  map_dfr(~ read_tsv(.x, col_types = cols()))

# sum numeric columns by Sample
si_summary <- si_df %>%
  group_by(Sample, Histology) %>%
  summarise(
    Total    = sum(Total, na.rm = TRUE),
    AS_neg   = sum(AS_neg, na.rm = TRUE),
    AS_pos   = sum(AS_pos, na.rm = TRUE),
    AS_total = sum(AS_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SI = AS_total / Total) %>%
  write_tsv(file.path(results_dir, "splicing_index.total.txt"))

# read in color palette
palette_file <- file.path(map_dir, "histologies-plot-group.tsv")

palette_df <- read_tsv(palette_file) %>%
  dplyr::rename(Histology = plot_group) %>%
  select(Histology, plot_group_hex, Sample = Kids_First_Biospecimen_ID) %>%
  unique()

binary_palette_file <- file.path(palette_dir, "binary_color_palette.tsv")
binary_palette <- read_tsv(binary_palette_file)
binary_palette <- setNames(binary_palette$hex_codes, binary_palette$color_names)

# Function to calculate medians and ranks
prepare_data_for_plot <- function(df, grouping_variable = NULL, min_samples = 5) {
  df %>%
    # Group by specified column
    group_by({{grouping_variable}}) %>%
    # Only keep groups with the specified minimum number of samples
    filter(n() > min_samples) %>%
    # Calculate group median
    mutate(
      group_median = median(SI, na.rm = TRUE),
      group_rank = rank(SI, ties.method = "first") / n(),
      sample_size = paste0("n = ", n())
    ) %>%
    ungroup() 
}

# create filenames for plots
file_si_total_plot = "sbi-plot-total-boxplot.pdf"
file_si_SE_plot = "sbi-plot-SE-boxplot.pdf"
file_si_RI_plot = "sbi-plot-RI-boxplot.pdf"
file_si_A5SS_plot = "sbi-plot-A5SS-boxplot.pdf"
file_si_A3SS_plot = "sbi-plot-A3SS-boxplot.pdf"

plot_sbi <- function(sbi_df, plot_file,label) {
  
  si_cdf_plot <- sbi_df %>%
    as_tibble() %>%
    select(Sample, SI) %>%
    left_join(palette_df) %>%
    drop_na(Histology) %>%
    # Perform calculations needed for plot
    prepare_data_for_plot(grouping_variable = Histology) %>%
    # remove "Other" cancer group
   # filter(Histology != "Other") %>%
    # Order cancer groups by median TMB
    mutate(Histology = str_wrap(Histology, 22),
           Histology = fct_reorder(Histology, SI, .fun = median)
    ) 
  
  ## compute quantiles to define high vs low SBI tumors
  quartiles_sbi <- quantile(sbi_df$SI, probs=c(.25, .75), na.rm = FALSE)
  IQR_sbi <- IQR(sbi_df$SI)
  lower_sbi <- quartiles_sbi[1]
  upper_sbi <- quartiles_sbi[2]
  
  si_cdf_plot <- si_cdf_plot %>% 
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High SBI", 
                               SI < lower_sbi ~ "Low SBI",
			       TRUE ~ "Middle SBI"))

  plot_colors <- si_cdf_plot$plot_group_hex
  names(plot_colors) <- si_cdf_plot$Histology
  
  p <- ggplot(si_cdf_plot) + 
    aes(
      x = Histology,
      y = SI 
    ) 
  
  p <- p +   
    # Yellow shading from 0 to lower_sbi
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = lower_sbi,
            fill = "#fac50c", alpha = 0.3) +
    geom_hline(yintercept = lower_sbi, linetype = "dashed", color = "#fac50c", size = 1) +

    # Blue shading from upper_sbi to top of the plot
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = upper_sbi, ymax = Inf,
  	    fill = "#2a81e2", alpha = 0.1) +
    geom_hline(yintercept = upper_sbi, linetype = "dashed", color = "#2a81e2", size = 1) +

    # Data
    geom_boxplot(aes(fill = Histology), color = "black", alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
    geom_sina(aes(fill = Histology, group = Histology),
              shape = 21, size = 2, color = "black", alpha = 0.7) +
    
    # Labels
    labs(x = "Histology",
         y = bquote(bold(.(label) * " Splicing Burden Index"))
    ) +
    theme_Publication() + 
    theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle = 70, hjust = 1)) +

    # Histology fill
    scale_fill_manual(values = plot_colors)
  
  if(label == "Total") {
    # Scale to max value (for main figure)
    p2 <- p + coord_cartesian(ylim = c(0, max(si_cdf_plot$SI)+0.01)) +
	  scale_y_continuous(breaks = c(0, lower_sbi, upper_sbi, 0.1, 0.15, 0.2),
                       labels = c("0", "Low SBI", "High SBI", "0.1", "0.15", "0.2"))

    ggsave(filename = paste0("scale-",plot_file), path = plots_dir, plot = p2,
         height = 6, width = 8, useDingbats = FALSE)
  }

  # Scale to same max value and sqrt y-axis to see lower values (for supp figure)
  p <- p + coord_cartesian(ylim = c(0, 0.4)) +
	  scale_y_continuous(trans = "sqrt",
			                 breaks = c(0, lower_sbi, upper_sbi, 0.1, 0.2, 0.3, 0.4),
                       labels = c("0", "Low SBI", "High SBI", "0.1", "0.2", "0.3", "0.4"))

  # Save plots
  ggsave(filename = plot_file, path = plots_dir, plot = p,
         height = 6, width = 8, useDingbats = FALSE)
}

## plot SBI for each splicing case
plot_sbi(si_summary,file_si_total_plot,"Total")
plot_sbi(splice_index_SE_df,file_si_SE_plot,"SE")
plot_sbi(splice_index_RI_df,file_si_RI_plot,"RI")
plot_sbi(splice_index_A5SS_df,file_si_A5SS_plot,"A5SS")
plot_sbi(splice_index_A3SS_df,file_si_A3SS_plot,"A3SS")


