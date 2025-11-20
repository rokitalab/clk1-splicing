## Estimate Brain cell type proportion for stranded samples
## Patricia Sullivan 
## 2025
## 
## Usage: Rscript --vanilla 06-estimate-cell-type.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(Biobase)
  library(BRETIGEA)
  library(lme4)
  library(performance)
  library(ggforce)
})


root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme_for_plots.R"))

## define file paths
expr_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")
stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")
hist_file <- file.path(root_dir, "analyses", "cohort_summary", "results", "histologies-plot-group.tsv")

## load expression data for stranded samples
hist <- read_tsv(hist_file) %>%
  select(Kids_First_Biospecimen_ID, plot_group, plot_group_hex)
stranded_clusters <- read_tsv(stranded_cluster_file)
expr <- readRDS(expr_file) %>%
  select(any_of(stranded_clusters$sample_id))

## run BRETIGEA
ct_res <- brainCells(expr, nMarker = 50, species = "human") %>%
  as.data.frame() %>%
  rename(Astrocytes = ast,
         Endothelial = end,
         Microglia = mic,
         Neurons = neu,
         Oligodendrocytes = oli,
         OPCs = opc) %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID") %>%
  left_join(hist) 

cell_types_long <- ct_res %>%
  pivot_longer(cols = c(Astrocytes, Endothelial, Microglia, Neurons, Oligodendrocytes, OPCs),
               names_to = "cell_type", values_to = "proportion_estimates")

## Plot cell types
plot_cols <- hist %>%
  distinct(plot_group, plot_group_hex) %>%
  deframe()

cell_prop_plot <- ggplot(cell_types_long, aes(x = plot_group, y = proportion_estimates, fill = plot_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina(alpha = 0.7, shape = 21) +
  theme_Publication() +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "none",
        panel.background = element_blank(),
  ) +
  scale_fill_manual(values = plot_cols) +
  xlab("Histology") +
  ylab("Cell type proportion estimate")

pdf(file = file.path(plot_dir, "cell-proportion-estimate-stranded.pdf"), height = 8, width = 12) 
  cell_prop_plot
dev.off()


