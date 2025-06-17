################################################################################
# 07-plot_sbi_vs_splicing_factor_expr.R
# script that calculates pearson correlation coefficients between SE SBI and splicing factor (SF) gene expression within identified clusters and tumor histologies
#
# written by Ryan Corbett
#
# usage: Rscript 07-plot_sbi_vs_splicing_factor_expr.R
################################################################################


suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("circlize")
  library("RColorBrewer")
  library("ComplexHeatmap")
  library("ggpubr")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

source(file.path(root_dir, "figures", "theme_for_plots.R"))

# set file paths
sbi_se_file <- file.path(results_dir,
                         "splicing_index.SE.txt")

expr_file <- file.path(data_dir,
                       "gene-expression-rsem-tpm-collapsed.rds")

cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", 
                          "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

sf_file <- file.path(root_dir, "analyses",
                     "splicing-factor_dysregulation",
                     "input", "splicing_factors.txt")

rmats_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")

# Wrangle data
sbi_se <- read_tsv(sbi_se_file)

expr <- readRDS(expr_file)

cluster_df <- read_tsv(cluster_file)

sfs <- read_lines(sf_file)
# some sfs have trailing "$" at end of gene name. Remove: 
sfs <- str_remove(sfs, "[$]")

# append cluster assignments to SBI df:
sbi_se <- sbi_se %>%
  left_join(cluster_df, by = c("Sample" = "sample_id"))

# subset expr matrix to only include SFs:
sf_expr <- expr[rownames(expr) %in% sfs,]

# define empty matrix to store SBI-SF expr correlation coefficients and p-values for each cluster
cor_mat <- matrix(0,
                  nrow(sf_expr),
                  length(unique(cluster_df$cluster)),
                  dimnames = list(rownames(sf_expr),
                                  sort(unique(cluster_df$cluster))))

p_mat <- cor_mat

# loop through clusters
for (clust in colnames(cor_mat)){
  
  # subset sbi df for cluster of interest
  sbi_sub <- sbi_se %>%
    dplyr::filter(cluster == clust)

  # loop through SFs
  for (sf in rownames(cor_mat)){
    
    # calculate pearson correlation coefficients and p-values
    cor_mat[sf,clust] <- cor.test(sbi_sub$SI, unlist(sf_expr[sf,sbi_sub$Sample]))$estimate
    p_mat[sf,clust] <- cor.test(sbi_sub$SI, unlist(sf_expr[sf,sbi_sub$Sample]))$p.value
    
  }
  
}

# calculate by-cluster false discovery rates
fdr_mat <- apply(p_mat, 2, function(x) p.adjust(x, "BH"))

# For plotting, we can limit to most significantly correlated SFs
plot_mat <- cor_mat[rowSums(fdr_mat < 1e-6 & !is.na(fdr_mat)) > 0,]
fdr_mat <- fdr_mat[rownames(plot_mat),]

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create the heatmap
cluster_heatmap <- Heatmap(plot_mat,
                          name = "Pearson\nCorrelation",
                          col = col_fun,
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          cluster_rows = TRUE,
                          cluster_columns = TRUE,
                          cell_fun = function(j, i, x, y, w, h, fill) {
                            if(fdr_mat[i, j] < 1e-5) {
                              grid.text("**", x, y)
                            } else if(fdr_mat[i, j] < 0.05) {
                              grid.text("*", x, y)
                            }
                          }
                          )

# save to output
pdf(NULL)
pdf(file.path(plots_dir, "sbi-sf-correlation-heatmap-byCluster.pdf"),
    width = 6, height = 10)

print(cluster_heatmap)

dev.off()

# repeat for by-histology correlation analyses
hist_cor_mat <- matrix(0,
                  nrow(sf_expr),
                  length(unique(cluster_df$plot_group)),
                  dimnames = list(rownames(sf_expr),
                                  sort(unique(cluster_df$plot_group))))

hist_p_mat <- hist_cor_mat

# loop through histologies
for (hist in colnames(hist_cor_mat)){
  
  # subset for hist of interest
  sbi_sub <- sbi_se %>%
    dplyr::filter(Histology == hist,
                  Sample %in% cluster_df$sample_id)
  
  # loop through SFs
  for (sf in rownames(hist_cor_mat)){
    
    # calculate pearson correlation coefficients and p-values
    hist_cor_mat[sf,hist] <- cor.test(sbi_sub$SI, unlist(sf_expr[sf,sbi_sub$Sample]))$estimate
    hist_p_mat[sf,hist] <- cor.test(sbi_sub$SI, unlist(sf_expr[sf,sbi_sub$Sample]))$p.value
    
  }
  
}

# calcualte FDRs
hist_fdr_mat <- apply(hist_p_mat, 2, function(x) p.adjust(x, "BH"))

hist_plot_mat <- hist_cor_mat[rowSums(hist_fdr_mat < 1e-7 & !is.na(hist_fdr_mat)) > 0,]
hist_fdr_mat <- hist_fdr_mat[rownames(hist_plot_mat),]

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Generate matrix of asterisks for significant p-values

# Create the heatmap
hist_heatmap <- Heatmap(t(hist_plot_mat),
                        name = "Pearson\nCorrelation",
                        col = col_fun,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(t(hist_fdr_mat)[i, j] < 1e-5) {
                            grid.text("**", x, y)
                          } else if(t(hist_fdr_mat)[i, j] < 0.05) {
                            grid.text("*", x, y)
                          }
                        }
                )

# Write to output
pdf(file.path(plots_dir, "sbi-sf-correlation-heatmap-byHist.pdf"),
    width = 10, height = 6)

print(hist_heatmap)

dev.off()

# extract CLK1 expression and append cluster, sbi info
clk1_expr <- expr %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  dplyr::filter(gene_name == "CLK1") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(CLK1_tpm = CLK1) %>%
  rownames_to_column("Sample") %>%
  left_join(cluster_df %>% dplyr::select(sample_id,
                                         cluster,
                                         plot_group,
                                         plot_group_hex),
            by = c("Sample" = "sample_id")) %>%
  dplyr::filter(!is.na(cluster)) %>%
  left_join(sbi_se %>% dplyr::select("Sample", "SI"))

# define plot group palette
plotgroup_palette <- unique(clk1_expr$plot_group_hex)
names(plotgroup_palette) <- unique(clk1_expr$plot_group)

# plot SBI against CLK1 expression in cluster 6
clk1_expr %>%
  dplyr::filter(cluster == 6) %>%
  ggplot(aes(x = log2(CLK1_tpm+1), y = log2(SI))) +
  geom_point(aes(color = plot_group)) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(Log[2] ~ "CLK1 TPM")),
       y = expression(bold(Log[2] ~ "SE SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 2, label.y = -3.5, size = 3) +
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-7.5,-3)) +
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-cluster6.pdf"),
       width = 7, height = 4)

# plot SBI against CLK1 expression in all other clusters
clk1_expr %>%
  dplyr::filter(cluster != 6) %>%
  ggplot(aes(x = log2(CLK1_tpm+1), y = log2(SI))) +
    geom_point(aes(color = plot_group),
               alpha = 0.7) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", 
                colour = "red",
                fill = "pink",
                linetype="dashed") +
    labs(x = expression(bold(Log[2] ~ "CLK1 TPM")),
         y = expression(bold(Log[2] ~ "SE SBI")),
         color = "Histology") + 
    stat_cor(method = "pearson",
             label.x = 0, label.y = -1, size = 3) +
    facet_wrap(~cluster, nrow = 2,
               labeller = labeller(cluster = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-8,0)) +
    theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-other-clusters.pdf"),
       width = 12, height = 5)

## Finally, assess correlations between CLK1 exon 4 PSI and SBI

# Load clk1 rmats and filter for exon 4
rmats_df <-  vroom(rmats_file) %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample_id, geneSymbol, IncLevel1) 

# append clk1 ex4 PSI to clk1 expr
clk1_expr <- clk1_expr %>%
  left_join(rmats_df %>% dplyr::select(sample_id, IncLevel1),
            by = c("Sample" = "sample_id"))

# Plot SBI against clk1 ex4 PSI in cluster 6
clk1_expr %>%
  dplyr::filter(cluster == 6) %>%
  ggplot(aes(x = IncLevel1, y = log2(SI))) +
  geom_point(aes(color = plot_group)) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = "CLK1 exon 4 PSI",
       y = expression(bold(Log[2] ~ "SE SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0.1, label.y = -3.5, size = 3) +
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-7.5,-3)) +
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-cluster6.pdf"),
       width = 7, height = 4)

# plot SBI vs. clk1 ex4 PSI in cluster 6
clk1_expr %>%
  dplyr::filter(cluster != 6) %>%
  ggplot(aes(x = IncLevel1, y = log2(SI))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = "CLK1 exon 4 PSI",
       y = expression(bold(Log[2] ~ "SE SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0.1, label.y = -2, size = 3) +
  facet_wrap(~cluster, nrow = 2,
             labeller = labeller(cluster = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-8,-1)) +
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-ex4-psi-other-clusters.pdf"),
       width = 12, height = 5)

# print session info
sessionInfo()
