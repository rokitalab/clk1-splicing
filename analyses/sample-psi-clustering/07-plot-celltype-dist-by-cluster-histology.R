## Plot cell type distributions by cluster and histology
## Jo Lynne Rokita Sullivan 
## 2025
## 
## Usage: Rscript --vanilla 07-plot-celltype-dist-by-cluster-histology.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# read in files
source(file.path(root_dir, "figures", "theme_for_plots.R"))
types <- read_tsv(file.path(results_dir, "cell-proportion-estimate-stranded.tsv"))
plot_mapping <- read_tsv(file.path(root_dir, "analyses", "cohort_summary", "results", "histologies-plot-group.tsv"))
cluster_df <- read_tsv(file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")) %>%
  dplyr::select(Kids_First_Biospecimen_ID = sample_id, cluster)

# Identify cell type columns
cell_types <- c("Astrocytes", "Endothelial", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")

# Long format
df_long <- types %>%
  left_join(plot_mapping[,c("Kids_First_Biospecimen_ID", "plot_group", "plot_group_hex")]) %>%
  left_join(cluster_df) %>%
  pivot_longer(
    cols = all_of(cell_types), 
    names_to = "CellType", 
    values_to = "Score"
  ) %>%
  mutate(cluster = as.factor(cluster))

# set colors
group_colors <- plot_mapping %>%
  distinct(plot_group, plot_group_hex) %>%
  arrange(plot_group)

hist_colors <- setNames(group_colors$plot_group_hex, group_colors$plot_group)

pdf(file.path(plot_dir, "celltype_dist_by_cluster_histology.pdf"), height = 10, width = 16)
ggplot(df_long,
       aes(x = Score, y = CellType, fill = plot_group)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~ cluster, scales = "free_x", nrow = 4) +
  scale_fill_manual(values = hist_colors) +
  labs(fill = "Histology") +
  theme_Publication()
dev.off()

