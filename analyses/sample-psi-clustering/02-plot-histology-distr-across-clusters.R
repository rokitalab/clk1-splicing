################################################################################
# 02-plot-histology-distr-across-clusters.R
# Plot the distribution of histologies across clusters
#
# Author: Ammar Naqvi, Ryan Corbett
################################################################################

## libraries 
suppressPackageStartupMessages({
  library("tidyverse")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")

input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## filepaths 
stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")
polyA_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-poly-A stranded.tsv")

subtype_hex <- file.path(input_dir, "subtype_hex.tsv")

# Wrangle data
histologies <- read_tsv(file.path(data_dir, "v11", "histologies-plot-group.tsv"))
subtype_hex_codes <- read_tsv(subtype_hex)

# define list of cluster file paths
cluster_files <- c(stranded_cluster_file,
                   polyA_cluster_file)
names(cluster_files) <- c("stranded", "poly-A-stranded")

# loop through cluster files to plot histology and subtype distributions
for (type in names(cluster_files)){
  
  pdf(NULL)

  # read in palette and cluster file
  cluster_df <- read_tsv(cluster_files[type]) %>%
    dplyr::rename(Kids_First_Biospecimen_ID = sample_id)
  
  histologies_df <- cluster_df %>%
    mutate(molecular_subtype_display = case_when(grepl("KIAA", molecular_subtype) ~ "KIAA1549--BRAF",
                                                 grepl("LGG, BRAF V600|LGG, RTK, BRAF V600E", molecular_subtype) ~ "BRAF V600E",
                                                 grepl("LGG, IDH|LGG, other MAPK, IDH", molecular_subtype) ~ "IDH",
                                                 molecular_subtype %in% c("LGG, wildtype") ~ "LGG, Wildtype",
                                                 grepl("SEGA", molecular_subtype) ~ "SEGA",
                                                 grepl("To be classified", molecular_subtype) ~ "To be classified",
                                                 grepl("LGG", molecular_subtype) & grepl("other MAPK", molecular_subtype) ~ "Other MAPK",
                                                 grepl("LGG", molecular_subtype) ~ "Other alteration",
                                                 grepl("HGG, IDH", molecular_subtype) ~ "IDH",
                                                 grepl("HGG, H3 wildtype", molecular_subtype) ~ "H3 wildtype",
                                                 grepl("HGG, PXA", molecular_subtype) ~ "PXA",
                                                 grepl("IHG", molecular_subtype) ~ "IHG",
                                                 grepl("H3 G35", molecular_subtype) ~ "H3 G35",
                                                 grepl("H3 K28", molecular_subtype) ~ "H3 K28",
                                                 is.na(molecular_subtype) & plot_group == "Other high-grade glioma" ~ "Oligodendroglioma",
                                                 TRUE ~ molecular_subtype)) %>%
    left_join(histologies %>% dplyr::select(plot_group, plot_group_hex)) %>%
    distinct()
  
  color_df <- histologies_df %>%
    dplyr::select(plot_group_hex, plot_group) %>%
    dplyr::filter(!is.na(plot_group)) %>%
    unique()
  cols <- as.character(color_df$plot_group_hex)
  names(cols) <- as.character(color_df$plot_group)
  
  # create plot
  pdf(file.path(paste0(plots_dir, "/cluster_membership_", type, ".pdf")),
      height = 4, width = 7)
  
  hist_plot <- ggplot(histologies_df, aes(fill=plot_group, x= factor(cluster))) +
                geom_bar(stat="count", position="stack") + 
                xlab("Cluster") + ylab("Frequency") +
                scale_fill_manual("Histology", values = cols) + 
                theme_Publication()
  
  print(hist_plot)
  
  dev.off()
  
  ## stratify by molecular subtype
  # subset to only those plot groups: ATRT, MB, LGG, HGG, EPN
  hist_subset <- histologies_df %>%
    filter(plot_group %in% c("Medulloblastoma", 
                             "Atypical Teratoid Rhabdoid Tumor",
                             "Low-grade glioma",
                             "Other high-grade glioma",
                             "DIPG or DMG")) %>%
    mutate(plot_group = case_when(plot_group %in% c("Other high-grade glioma",
                                                    "DIPG or DMG") ~ "High-grade glioma",
                                  TRUE ~ plot_group))
  
  cols <- subtype_hex_codes$hex_code
  names(cols) <- subtype_hex_codes$molecular_subtype_display
  
  pdf(file.path(paste0(plots_dir, "/cluster_membership-subtypes_", type, ".pdf")),
      height = 6, width = 11)
  
  subtype_plot <- hist_subset %>% 
    ggplot(aes(fill=molecular_subtype_display, x= factor(cluster))) +
    geom_bar(stat="count", position="stack", color = "black", size = 0.2) + 
    facet_wrap(~plot_group, nrow = 2, scales = "free_y") +
    xlab("Cluster") + 
    ylab("Frequency") +
    theme_Publication() +
    scale_fill_manual(name = "Molecular Subtype", values = cols) 
  
  print(subtype_plot)
  
  dev.off()
  
}

# Session info
sessionInfo()