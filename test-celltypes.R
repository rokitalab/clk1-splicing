source("clk1-splicing/figures/theme_for_plots.R")
types <- read_tsv("~/clk1-splicing/analyses/sample-psi-clustering/results/cell-proportion-estimate-stranded.tsv")
clk1_psi <- read_tsv("~/clk1-splicing/data/clk1-splice-events-rmats.tsv") %>%
  filter(exonEnd=="200860215") %>%
  select(sample_id, IncLevel1) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id,
                psi = IncLevel1)

hist <- read_tsv("~/clk1-splicing/analyses/cohort_summary/results/histologies-plot-group.tsv")
clusters <- read_tsv("clk1-splicing/analyses/sample-psi-clustering/results/sample-cluster-metadata-top-5000-events-stranded.tsv") %>%
  select(sample_id, cluster)  %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id)


comb <- types %>%
  left_join(clk1_psi) %>%
  left_join(hist) %>%
  left_join(clusters)


library(tidyverse)

# Identify cell type columns
cell_types <- c("Astrocytes", "Endothelial", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")   # replace with yours

# Calculate correlations
cor_df <- map_df(cell_types, ~{
  tibble(
    cell_type = .x,
    correlation = cor(comb$psi, comb[[.x]], use = "complete.obs", method = "pearson")
  )
})

cor_df



# 1. Long format
df_long <- comb %>%
  pivot_longer(
    cols = all_of(cell_types), 
    names_to = "CellType", 
    values_to = "Score"
  )

# 2. Safe color mapping from plot_group -> hex
group_colors <- comb %>%
  distinct(plot_group, plot_group_hex) %>%
  arrange(plot_group)   # optional, just for consistency

colors <- setNames(group_colors$plot_group_hex, group_colors$plot_group)

# 3. Basic plot (your version, with fixed colors)
ggplot(df_long, aes(x = Score, y = psi, color = plot_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_wrap(~ CellType, scales = "free_x") +
  scale_color_manual(values = colors) +
  theme_bw()

# subset HGG

# 3. Basic plot (your version, with fixed colors)
df_long_hgg <- df_long %>%
  filter(plot_group %in% c("Diffuse midline glioma", "Other high-grade glioma"))
ggplot(df_long_hgg, aes(x = Score, y = psi, color = plot_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_wrap(~ CellType, scales = "free_x") +
  scale_color_manual(values = colors) +
  theme_Publication()


df_long_cor <- df_long_hgg %>%
  group_by(CellType) %>%
  mutate(
    cor_label = paste0(
      "r = ", round(cor(psi, Score, use = "complete.obs"), 2)
    )
  ) %>%
  ungroup()



ggplot(df_long_cor, aes(x = Score, y = psi, color = plot_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_wrap(~ CellType, scales = "free_x") +
  scale_color_manual(values = colors) +
  geom_text(
    aes(label = cor_label),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.3,
    size = 3
  ) +
  theme_Publication()



pca <- prcomp(comb[, cell_types], scale. = TRUE)

df_pca <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  psi = comb$psi,
  histology = comb$plot_group
)

ggplot(df_pca, aes(x = PC1, y = PC2, color = psi)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  theme_Publication()

library(ggridges)


ggplot(df_long,
       aes(x = Score, y = CellType, fill = plot_group)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~cluster, scales = "free_x", nrow = 4) +
  scale_fill_manual(values = colors) +
  theme_Publication()
library(ggh4x)



library(ggridges)
library(ggh4x)

cluster_colors <-  c("1" = "#B2DF8A",
                     "2" = "#E31A1C",
                     "3" = "#33A02C",
                     "4" = "#A6CEE3",
                     "5" = "#FB9A99",
                     "6" = "#FDBF6F",
                     "7" = "#CAB2D6",
                     "8" = "#FFFF99",
                     "9" = "#1F78B4",
                     "10" = "#B15928",
                     "11" = "#6A3D9A")
pdf("clk1-splicing/analyses/celltypes_by_cluster.pdf", height = 10, width = 16)  
ggplot(df_long,
       aes(x = Score, y = CellType, fill = plot_group)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap2(
    ~ cluster,
    scales = "free_x",
    nrow = 4,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill   = unname(cluster_colors)      ),
      text_x = elem_list_text(
        colour = "black",
        face   = "bold"
      )
    )
  ) +
  scale_fill_manual(values = colors) +   # your plot_group palette
  theme_Publication()
dev.off()



