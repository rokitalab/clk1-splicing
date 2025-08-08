################################################################################
# 08-plot-sbi-clusters.R
# script that takes in SBI tsv files and plots levels by cluster
#
# written by Ammar Naqvi, Jo Lynne Rokita, Patricia Sullivan
#
# usage: Rscript 08-plot-sbi-clusters.R
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

# read in clusters
cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")
cluster_df <- read_tsv(cluster_file) %>%
  rename(Sample = sample_id)

# read in SI files
splice_index_SE_file   <- file.path(results_dir, "splicing_index.SE.txt") 
splice_index_RI_file   <- file.path(results_dir, "splicing_index.RI.txt")
splice_index_A5SS_file <- file.path(results_dir, "splicing_index.A5SS.txt")
splice_index_A3SS_file <- file.path(results_dir, "splicing_index.A3SS.txt")

splice_index_SE_df   <- readr::read_tsv(splice_index_SE_file)
splice_index_RI_df   <- readr::read_tsv(splice_index_RI_file)
splice_index_A5SS_df <- readr::read_tsv(splice_index_A5SS_file)
splice_index_A3SS_df <- readr::read_tsv(splice_index_A3SS_file)

# read in color palette
palette_file <- file.path(map_dir, "histologies-plot-group.tsv")

palette_df <- read_tsv(palette_file) %>%
  dplyr::rename(Histology = plot_group) %>%
  select(Histology, plot_group_hex) %>%
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
file_si_SE_plot = "clusters-sbi-plot-SE-boxplot.pdf"
file_si_RI_plot = "clusters-sbi-plot-RI-boxplot.pdf"
file_si_A5SS_plot = "clusters-sbi-plot-A5SS-boxplot.pdf"
file_si_A3SS_plot = "clusters-sbi-plot-A3SS-boxplot.pdf"

# define colors for clusters
cluster_cols <- c("#B2DF8A","#E31A1C","#33A02C","#A6CEE3","#FB9A99","#FDBF6F",
                  "#CAB2D6","#FFFF99","#1F78B4","#B15928","#6A3D9A","#FF7F00",
                  "#2ef4ca","#f4cced","#bd18ea")
names(cluster_cols) <- 1:length(cluster_cols)
cluster_cols <- cluster_cols[1:length(unique(cluster_df$cluster))]

plot_sbi <- function(sbi_df, plot_file,label) {
  
  si_cdf_plot <- sbi_df %>%
    as_tibble() %>%
    inner_join(cluster_df) %>%
    select(Sample, SI, cluster, plot_group, plot_group_hex, Histology) %>%
    # Perform calculations needed for plot
    prepare_data_for_plot(grouping_variable = cluster) %>%
    # remove "Other" cancer group
    # filter(Histology != "Other") %>%
    # Order cancer groups by median TMB
    mutate(cluster = factor(cluster)) 
  
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
  plot_colors <- plot_colors[sort(names(plot_colors))]
  
  p <- ggplot(si_cdf_plot) + 
    aes(
      x = cluster,
      y = SI 
    ) 
    
    # Data
  p <- p +   
    # Yellow shading from 0 to lower_sbi
    #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = lower_sbi,
    #        fill = "#fac50c", alpha = 0.3) +
    #geom_hline(yintercept = lower_sbi, linetype = "dashed", color = "#fac50c", size = 1) +
    
    # Blue shading from upper_sbi to top of the plot
    #annotate("rect", xmin = -Inf, xmax = Inf, ymin = upper_sbi, ymax = Inf,
    #    fill = "#2a81e2", alpha = 0.1) +
    #geom_hline(yintercept = upper_sbi, linetype = "dashed", color = "#2a81e2", size = 1) +
    
    # Data
    geom_boxplot(aes(fill = cluster), color = "black", alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
    geom_sina(aes(fill = plot_group, group = cluster),
              shape = 21, size = 2, color = "black", alpha = 0.7, show.legend = TRUE) +
    
    # Labels
    labs(x = "Cluster",
         y = bquote(bold(.(label) * " Splicing Burden Index"))
    ) +
    theme_Publication() + 
    theme(legend.position = "right") +
    labs(fill = "Histology") +
    
    # Histology fill
    scale_fill_manual(
      values = c(plot_colors, cluster_cols),
      breaks = names(plot_colors)  # ensures only plot_group shows
    )
  
  # Save plots
  ggsave(filename = plot_file, path = plots_dir, plot = p,
         height = 6, width = 10, useDingbats = FALSE)
}

## plot SBI for each splicing case
plot_sbi(splice_index_SE_df,file_si_SE_plot,"SE")
plot_sbi(splice_index_RI_df,file_si_RI_plot,"RI")
plot_sbi(splice_index_A5SS_df,file_si_A5SS_plot,"A5SS")
plot_sbi(splice_index_A3SS_df,file_si_A3SS_plot,"A3SS")
