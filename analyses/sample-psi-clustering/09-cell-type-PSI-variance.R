## Calculate the PSI variance due to cell type
## Patricia Sullivan 
## 2025
## 
## Usage: Rscript --vanilla 09-cell-type-PSI-variance.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(performance)
  library(ggsci)
})

## define directory locations

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme_for_plots.R"))

# define file paths
psi_file <- file.path(results_dir, "psi-matrix-top-5000-events-stranded.rds")
cell_types_file <- file.path(results_dir, "cell-proportion-estimate-stranded.tsv")

# load files
psi_mat <- readRDS(psi_file)
cell_types <- read_tsv(cell_types_file) %>% 
  column_to_rownames("Kids_First_Biospecimen_ID")

# align PSI matrix with metadata
common <- intersect(colnames(psi_mat), rownames(cell_types))
psi_mat  <- psi_mat[, common]
cell_types <- cell_types[common, ]

# define the formula to model
form <- PSI ~ Astrocytes + Endothelial + Microglia + Neurons + Oligodendrocytes + OPCs

variance_df <- data.frame()

# Loop through events
for (event in rownames(psi_mat)) {
  y <- t(psi_mat[event, ])
  data_sub <- cbind(cell_types, PSI = y) %>%
    rename(PSI = event)
  
  # drop missing samples
  data_sub <- data_sub[!is.na(data_sub$PSI), ]
  # run model on event only if there are enough entries
  if (nrow(data_sub) < 5 || sd(data_sub$PSI) == 0) next
  
  # build a model per splicing event with PSI and cell type information
  m <- try(lm(form, data = data_sub), silent = TRUE)
  if (inherits(m, "try-error")) next
  
  # use ANOVA SS to get % variance explained
  a <- try(anova(m), silent = TRUE)
  if (inherits(a, "try-error")) next
  
  ss <- a[,"Sum Sq"]
  names(ss) <- rownames(a)
  # ensure Residuals included
  if (!"Residuals" %in% names(ss)) {
    ss["Residuals"] <- sum(residuals(m)^2)
  }
  tot <- sum(ss, na.rm = TRUE)
  frac <- 100 * ss / tot
  
  # Save event data
  all_cols <- c("Astrocytes", "Endothelial", "Microglia", "Neurons", "Oligodendrocytes", "OPCs", "Residuals")

  vals <- setNames(rep(0, length(all_cols)), all_cols)
  overlap <- intersect(names(frac), names(vals))
  vals[overlap] <- frac[overlap]
  
  out <- as.data.frame(t(vals))
  out$Event <- event
  out$n_samples_used <- nrow(data_sub)
  variance_df <- bind_rows(variance_df, out)
}

# mean variance 
vp_mean <- variance_df %>%
  select(-Event, -n_samples_used, -Residuals) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Mean Variance") %>%
  arrange(desc(`Mean Variance`)) 

pdf(file.path(plot_dir, "celltype_variance_contribution.pdf"), height = 4, width = 4)
ggplot(vp_mean, aes(x = Factor, y = `Mean Variance`, fill = Factor)) +
  geom_col(color = "black") +
  theme_Publication() +
  labs(
    x = "Cell type",
    y = "Mean variance (%)",
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_lancet()
dev.off()

# event-level cell type contribution
celltype_cols <- c("Astrocytes", "Endothelial", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")
variance_df$MaxCellVar <- apply(variance_df[, celltype_cols], 1, max, na.rm=TRUE)
variance_df$MaxCellType <- apply(variance_df[, celltype_cols, drop = FALSE], 1, function(x) names(x)[which.max(x)])
variance_df$TotalCellVar <- rowSums(variance_df[, celltype_cols])
variance_df$VarBinTotal <- cut(
  variance_df$TotalCellVar,
  breaks = c(-Inf, 10, 30, 50, Inf),
  labels = c("<10%\nMinimal", "10–30%\nModerate" ,"30-50%\nStrong","50%+\nDominant"),
  #labels = c("<10%", "10-30%","30-50%","50%+"),
  right = FALSE
)
variance_df$VarBinMax <- cut(
  variance_df$MaxCellVar,
  breaks = c(-Inf, 10, 30, 50, Inf),
  labels = c("<10%\nMinimal", "10–30%\nModerate" ,"30-50%\nStrong","50%+\nDominant"),
  #labels = c("<10%", "10-30%","30-50%","50%+"),
  right = FALSE
)

# save per-event variance df
write_tsv(variance_df, file = file.path(results_dir, "cell-type-variance-top-5000-events-stranded.tsv"))


# Total (sum of all cell-type variance)
plot_df <- variance_df %>%
  count(VarBinTotal, MaxCellType) %>%
  mutate(prop = n / sum(n))

pdf(file.path(plot_dir, "celltype_total_variance.pdf"), height = 4, width = 7)
ggplot(plot_df, aes(x = VarBinTotal, y = prop, fill = MaxCellType)) +
  geom_col(color = "black") +
  theme_Publication() +
  labs(
    x = "Total PSI variance explained by cell type",
    y = "Proportion of events",
    fill = "Max contributor"
  ) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_lancet()
dev.off()

# Max (maximum cell-type variance)
plot_df2 <- variance_df %>%
  count(VarBinMax, MaxCellType) %>% 
  mutate(prop = n / sum(n))

pdf(file.path(plot_dir, "celltype_max_variance.pdf"), height = 4, width = 7)
ggplot(plot_df2, aes(x = VarBinMax, y = prop, fill = MaxCellType)) +
  geom_col(color = "black") +
  theme_Publication() +
  labs(
    x = "Max PSI variance explained by cell type",
    y = "Proportion of events",
    fill = "Max contributor"
  ) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_lancet()
dev.off()
