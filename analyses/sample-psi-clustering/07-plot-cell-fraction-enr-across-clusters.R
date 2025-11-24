################################################################################
# 07-plot-cell-fraction-enr-across-clusters.R
# Plot the distribution of high- and low-expressed cell type groups across splicing clusters
#
# Author: Ryan Corbett
################################################################################

## libraries 
suppressPackageStartupMessages({
  library("tidyverse")
  library("vroom")
  library("circlize")
  library("ComplexHeatmap")
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
source(file.path(analysis_dir, "util", "heatmap_function.R"))

# Set file paths
## filepaths 
stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")

cell_fraction_file <- file.path(results_dir,
                                "cell-proportion-estimate-stranded.tsv")


## Wrangle data
stranded_cluster_df <- read_tsv(stranded_cluster_file)

cell_fraction_df <- read_tsv(cell_fraction_file)

# write function to assign samples to quartile groups based on cell fraction estimates
quartile_label <- function(x, prefix) {
  qs <- summary(x)
  
  dplyr::case_when(
    x >= qs["3rd Qu."] ~ "4th Qu.",
    x >= qs["Median"]  ~ "3rd Qu.",
    x >= qs["1st Qu."] ~ "2nd Qu.",
    TRUE               ~ "1st Qu."
  )
}

# assign samples to cell fraction groups using newly-defined function
cell_fraction_df <- read_tsv(cell_fraction_file) %>%
  mutate(across(
    c(Astrocytes, Endothelial, Microglia, Neurons,
      Oligodendrocytes, OPCs),
    ~ quartile_label(.x, cur_column()),
    .names = "{.col}_group"
  ))

# join cell fraction group assignments to cluster df
stranded_cluster_df <- stranded_cluster_df %>%
  left_join(cell_fraction_df %>% 
              dplyr::select(contains(c("Kids_First_Biospecimen_ID", "group"))),
            by = c("sample_id" = "Kids_First_Biospecimen_ID"))

# Generate cluster enrichment plots for each cell type
cluster_astro_enr <- plot_enr(stranded_cluster_df, 
                            "cluster", "Astrocytes_group",
                            sort(unique(stranded_cluster_df$cluster)),
                            sort(unique(stranded_cluster_df$Astrocytes_group)),
                            padjust = TRUE)

cluster_endo_enr <- plot_enr(stranded_cluster_df, 
                              "cluster", "Endothelial_group",
                              sort(unique(stranded_cluster_df$cluster)),
                              sort(unique(stranded_cluster_df$Endothelial_group)),
                              padjust = TRUE)

cluster_micro_enr <- plot_enr(stranded_cluster_df, 
                             "cluster", "Microglia_group",
                             sort(unique(stranded_cluster_df$cluster)),
                             sort(unique(stranded_cluster_df$Microglia_group)),
                             padjust = TRUE)

cluster_neuro_enr <- plot_enr(stranded_cluster_df, 
                              "cluster", "Neurons_group",
                              sort(unique(stranded_cluster_df$cluster)),
                              sort(unique(stranded_cluster_df$Neurons_group)),
                              padjust = TRUE)

cluster_oligo_enr <- plot_enr(stranded_cluster_df, 
                              "cluster", "Oligodendrocytes_group",
                              sort(unique(stranded_cluster_df$cluster)),
                              sort(unique(stranded_cluster_df$Oligodendrocytes_group)),
                              padjust = TRUE)

cluster_opc_enr <- plot_enr(stranded_cluster_df, 
                              "cluster", "OPCs_group",
                              sort(unique(stranded_cluster_df$cluster)),
                              sort(unique(stranded_cluster_df$OPCs_group)),
                              padjust = TRUE)

# post-hoc add titles to each plot
cluster_astro_enr@column_title <- "Astrocytes"
cluster_endo_enr@column_title  <- "Endothelial"
cluster_micro_enr@column_title <- "Microglia"
cluster_neuro_enr@column_title <- "Neurons"
cluster_oligo_enr@column_title <- "Oligodendrocytes"
cluster_opc_enr@column_title   <- "OPCs"

# create object of merged complexheatmap objects 
ht_list <- cluster_astro_enr + cluster_endo_enr + cluster_micro_enr + cluster_neuro_enr + cluster_oligo_enr + cluster_opc_enr

# Draw & save plot
pdf(file.path(plots_dir, 
              "brain-cell-type-fraction-cluster-enr-heatmap-stranded.pdf"),
    width = 14, height = 6)

draw(ht_list)

dev.off()

hgg_cluster_df <- stranded_cluster_df %>%
  dplyr::filter(plot_group %in% c("Other high-grade glioma", "DIPG or DMG"))

# Generate cluster enrichment plots for each cell type
hgg_astro_enr <- plot_enr(hgg_cluster_df, 
                              "cluster", "Astrocytes_group",
                              sort(unique(hgg_cluster_df$cluster)),
                              sort(unique(hgg_cluster_df$Astrocytes_group)),
                              padjust = TRUE)

hgg_endo_enr <- plot_enr(hgg_cluster_df, 
                             "cluster", "Endothelial_group",
                             sort(unique(hgg_cluster_df$cluster)),
                             sort(unique(hgg_cluster_df$Endothelial_group)),
                             padjust = TRUE)

hgg_micro_enr <- plot_enr(hgg_cluster_df, 
                              "cluster", "Microglia_group",
                              sort(unique(hgg_cluster_df$cluster)),
                              sort(unique(hgg_cluster_df$Microglia_group)),
                              padjust = TRUE)

hgg_neuro_enr <- plot_enr(hgg_cluster_df, 
                              "cluster", "Neurons_group",
                              sort(unique(hgg_cluster_df$cluster)),
                              sort(unique(hgg_cluster_df$Neurons_group)),
                              padjust = TRUE)

hgg_oligo_enr <- plot_enr(hgg_cluster_df, 
                              "cluster", "Oligodendrocytes_group",
                              sort(unique(hgg_cluster_df$cluster)),
                              sort(unique(hgg_cluster_df$Oligodendrocytes_group)),
                              padjust = TRUE)

hgg_opc_enr <- plot_enr(hgg_cluster_df, 
                            "cluster", "OPCs_group",
                            sort(unique(hgg_cluster_df$cluster)),
                            sort(unique(hgg_cluster_df$OPCs_group)),
                            padjust = TRUE)


# post-hoc add titles to each plot
hgg_astro_enr@column_title <- "Astrocytes"
hgg_endo_enr@column_title  <- "Endothelial"
hgg_micro_enr@column_title <- "Microglia"
hgg_neuro_enr@column_title <- "Neurons"
hgg_oligo_enr@column_title <- "Oligodendrocytes"
hgg_opc_enr@column_title   <- "OPCs"

# create object of merged complexheatmap objects 
hgg_ht_list <- hgg_astro_enr + hgg_endo_enr + hgg_micro_enr + hgg_neuro_enr + hgg_oligo_enr + hgg_opc_enr

# Draw & save plot
pdf(file.path(plots_dir, 
              "brain-cell-type-fraction-hgg-cluster-enr-heatmap-stranded.pdf"),
    width = 14, height = 4)

draw(hgg_ht_list)

dev.off()

# print session info
sessionInfo()
