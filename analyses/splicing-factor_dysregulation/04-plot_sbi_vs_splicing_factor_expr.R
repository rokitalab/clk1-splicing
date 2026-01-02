################################################################################
# 04-plot_sbi_vs_splicing_factor_expr.R
# script that calculates pearson correlation coefficients between Total SBI and splicing factor (SF) gene expression within identified clusters and tumor histologies
#
# written by Ryan Corbett, Patricia Sullivan
#
# usage: Rscript 04-plot_sbi_vs_splicing_factor_expr.R
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
  library("cowplot")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

source(file.path(root_dir, "figures", "theme_for_plots.R"))

# set file paths
sbi_total_file <- file.path(root_dir, "analyses", 
                            "splicing_index",
                            "results",
                            "splicing_index.total.txt")

expr_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")
expr_trans_file <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")

hist_file <- file.path(file.path(data_dir, "histologies-plot-group.tsv"))

cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", 
                          "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

sf_file <- file.path(input_dir, "splicing_factors.txt")

rmats_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")

# Wrangle data
sbi_df <- read_tsv(sbi_total_file)
hist_df <- read_tsv(hist_file) %>%
  select(sample_id = Kids_First_Biospecimen_ID, plot_group, plot_group_hex)

expr <- readRDS(expr_file)
expr_trans <- readRDS(expr_trans_file)

cluster_df <- read_tsv(cluster_file) %>%
  select(-c(plot_group, plot_group_hex)) %>%
  inner_join(hist_df)

sfs <- read_lines(sf_file)
# some sfs have trailing "$" at end of gene name. Remove: 
sfs <- str_remove(sfs, "[$]")

# append cluster assignments to SBI df:
sbi_df <- sbi_df %>%
  left_join(cluster_df, by = c("Sample" = "sample_id"))

# subset expr matrix to only include SFs:
sf_expr <- expr[rownames(expr) %in% sfs,]

vars <- list(clusters = unique(cluster_df$cluster),
             plot_groups = unique(cluster_df$plot_group))

for (group in names(vars)){
  
  groups <- vars[[group]]
  
  cor_mat <- matrix(0,
                    nrow(sf_expr),
                    length(groups),
                    dimnames = list(rownames(sf_expr),
                                    sort(groups)))
  
  p_mat <- cor_mat
  
  # loop through levels of group
  for (level in colnames(cor_mat)){
    
    if (group == "plot_groups"){
      
      sbi_sub <- sbi_df %>%
        dplyr::filter(plot_group == level)
      
    } else {
      
      sbi_sub <- sbi_df %>%
        dplyr::filter(cluster == level)
      
    }
    
    # loop through SFs
    for (sf in rownames(cor_mat)){
      
      # calculate pearson correlation coefficients and p-values
      cor_mat[sf,level] <- cor.test(log2(sbi_sub$SI),
                                    log2(unlist(sf_expr[sf,sbi_sub$Sample])),
                                    method = "pearson")$estimate
      p_mat[sf,level] <- cor.test(log2(sbi_sub$SI),
                                  log2(unlist(sf_expr[sf,sbi_sub$Sample])),
                                  method = "pearson")$p.value
      
    }
    
  }

  # calculate by-cluster false discovery rates
  fdr_mat <- apply(p_mat, 2, function(x) p.adjust(x, method = "BH"))
  fdr_mat <- as.matrix(fdr_mat)
  rownames(fdr_mat) <- rownames(p_mat)
  colnames(fdr_mat) <- colnames(p_mat)
  
  fdr_thresh <- 0.05
  r_cutoff   <- 0.50
  max_genes  <- 40
  
  # align fdr_mat and cor_mat (shared genes + shared clusters)
  common_genes <- intersect(rownames(fdr_mat), rownames(cor_mat))
  common_cols  <- intersect(colnames(fdr_mat), colnames(cor_mat))
  
  fdr2 <- fdr_mat[common_genes, common_cols, drop = FALSE]
  cor2 <- cor_mat[common_genes, common_cols, drop = FALSE]
  
  # genes that pass both thresholds in at least one cluster
  pass <- (fdr2 < fdr_thresh) & (abs(cor2) >= r_cutoff)
  keep <- rownames(fdr2)[rowSums(pass, na.rm = TRUE) > 0]
  
  # if still too many, rank by best (lowest) FDR among cells that also pass r cutoff
  if (length(keep) > max_genes) {
    best_fdr <- sapply(keep, function(g) {
      idx <- abs(cor2[g, ]) >= r_cutoff
      min(fdr2[g, idx], na.rm = TRUE)
    })
    keep <- names(sort(best_fdr))[1:max_genes]
  }
  
  gois_level <- keep
  
  plot_mat <- cor2[gois_level, , drop = FALSE]
  plot_mat <- plot_mat[rowSums(is.nan(plot_mat)) == 0, , drop = FALSE]
  
  fdr_plot_mat <- fdr2[rownames(plot_mat), colnames(plot_mat), drop = FALSE]
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # set some separate params
  if (group == "clusters") {
    set_width = 6
    set_height = 8
    set_rotation = 0
    set_title = "Cluster"
  }
  
  if (group == "plot_groups") {
    set_width = 7
    set_height = 11
    set_rotation = 70
    set_title = "Histology"
  }
  
  
  # Create the heatmap
  heatmap <- Heatmap(
    plot_mat,
    name = "Pearson\nCorrelation",
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontface = "italic"), 
    column_names_rot = set_rotation,
    column_title = set_title,
    column_title_side = "bottom",
    column_title_gp = gpar(fontface = "bold"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (fdr_plot_mat[i, j] < 1e-5) {
        grid.text("**", x, y)
      } else if (fdr_plot_mat[i, j] < 0.05) {
        grid.text("*", x, y)
      }
    }
  )
  
  # save to output
  pdf(NULL)
  pdf(file.path(plots_dir, glue::glue("sbi-sf-correlation-heatmap-by-{group}.pdf")),
      width = set_width, height = set_height)
  
  print(heatmap)
  
  dev.off()
  
  # Merge correlation stats and write to output
  sbi_sf_cor_df <- cor_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    pivot_longer(-gene_symbol,
                 names_to = group,
                 values_to = "pearson_r")
  
  sbi_sf_p_df <- p_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    pivot_longer(-gene_symbol,
                 names_to = group,
                 values_to = "pearson_pvalue")
  
  sbi_sf_fdr_df <- fdr_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    pivot_longer(-gene_symbol,
                 names_to = group,
                 values_to = "pearson_fdr")
  
  sbi_sf_cor_df <- sbi_sf_cor_df %>%
    left_join(sbi_sf_p_df) %>%
    left_join(sbi_sf_fdr_df) %>%
    filter(!is.na(pearson_r),
           pearson_fdr < 0.05) %>%
    arrange(desc(abs(pearson_r)))
  
  write_tsv(sbi_sf_cor_df,
            file.path(results_dir,
                      glue::glue("total-sbi-sf-expr-correlations-{group}.tsv")))
  
}

clk1_ex4_expr <- expr_trans %>%
  filter(grepl("^CLK1", gene_symbol)) %>%
  dplyr::filter(transcript_id == "ENST00000321356.9") %>%
  group_by(transcript_id) %>%
  summarise(across(starts_with("BS"), sum, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Sample",
    values_to = "CLK1_tpm"
  ) %>%
  dplyr::select(-transcript_id) %>%
  left_join(cluster_df %>% dplyr::select(sample_id,
                                         cluster,
                                         plot_group,
                                         plot_group_hex),
            by = c("Sample" = "sample_id")) %>%
  dplyr::filter(!is.na(cluster)) %>%
  left_join(sbi_df %>% dplyr::select("Sample", "SI"))

# extract CLK1 expression and append cluster, sbi info
clk1_expr <- sf_expr %>%
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
  left_join(sbi_df %>% dplyr::select("Sample", "SI"))

# define plot group palette
plotgroup_palette <- unique(clk1_expr$plot_group_hex)
names(plotgroup_palette) <- unique(clk1_expr$plot_group)

# plot SBI against CLK1 expression by cluster
clk1_expr %>%
  ggplot(aes(x = log2(CLK1_tpm), y = log2(SI))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(Log[2] ~ bolditalic("CLK1") ~ "TPM")),
       y = expression(bold(Log[2] ~ "SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0, label.y = -2.5, size = 3) +
  facet_wrap(~cluster, nrow = 2,
             labeller = labeller(cluster = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-by-cluster.pdf"),
       width = 13, height = 5)


clk1_expr %>%
  dplyr::filter(cluster == 7) %>%
  ggplot(aes(x = log2(CLK1_tpm), y = log2(SI))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(Log[2] ~ bolditalic("CLK1") ~ "TPM")),
       y = expression(bold(Log[2] ~ "SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 2, label.y = -2.5, size = 3) +
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-cluster7.pdf"),
       width = 7, height = 4)

# plot SBI against CLK1 expression by histology
clk1_expr %>%
  ggplot(aes(x = log2(CLK1_tpm+1), y = log2(SI))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7,
             show.legend = FALSE) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(Log[2] ~ bolditalic("CLK1") ~ "TPM")),
       y = expression(bold(Log[2] ~ "SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0, label.y = -2, size = 3) +
  facet_wrap(~plot_group, nrow = 3,
             labeller = labeller(plot_group = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-tpm-byHist.pdf"),
       width = 12, height = 7)

## Assess correlations between CLK1 exon 4 PSI and SBI

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

# append clk1 ex4 PSI to clk1 exon 4 expr
clk1_ex4_expr <- clk1_ex4_expr %>%
  left_join(rmats_df %>% dplyr::select(sample_id, IncLevel1),
            by = c("Sample" = "sample_id"))

# plot SBI vs. clk1 ex4 PSI by cluster
clk1_expr %>%
  ggplot(aes(x = IncLevel1, y = log2(SI))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(bolditalic("CLK1") ~ "exon 4 PSI")),
       y = expression(bold(Log[2] ~ "SBI")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0.1, label.y = -2.5, size = 3) +
  facet_wrap(~cluster, nrow = 2,
             labeller = labeller(cluster = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "sbi-vs-clk1-ex4-psi-by-cluster.pdf"),
       width = 12, height = 5)



## Assess correlations between CLK1 exon 4 PSI and CLK1 TPM

# plot CLK1 TPM vs. clk1 ex4 PSI by cluster
clk1_ex4_expr %>%
  ggplot(aes(x = IncLevel1, y = log2(CLK1_tpm + 1))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(bolditalic("CLK1") ~ "exon 4 PSI")),
       y = expression(bold(Log[2] ~ bolditalic("CLK1") ~ "ENST00000321356 TPM")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0, label.y = -.6, size = 3) +
  facet_wrap(~cluster, nrow = 2,
             labeller = labeller(cluster = label_wrap_gen(18))) + 
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-1,4)) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))

# save plot
ggsave(file.path(plots_dir,
                 "tpm-vs-clk1-ex4-psi-by-cluster.pdf"),
       width = 12, height = 5)


clk1_ex4_expr %>%
  dplyr::filter(cluster == 7) %>%
  ggplot(aes(x = IncLevel1, y = log2(CLK1_tpm + 1))) +
  geom_point(aes(color = plot_group),
             alpha = 0.7) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = expression(bold(bolditalic("CLK1") ~ "exon 4 PSI")),
       y = expression(bold(Log[2] ~ bolditalic("CLK1") ~ "ENST00000321356 TPM")),
       color = "Histology") + 
  stat_cor(method = "pearson",
           label.x = 0.1, label.y = 4.5, size = 3) +
  scale_color_manual(values = plotgroup_palette, breaks = names(plotgroup_palette)) + 
  ylim(c(-1,4.75)) +
  theme_Publication()

# save plot
ggsave(file.path(plots_dir,
                 "tpm-vs-clk1-ex4-psi-cluster7.pdf"),
       width = 7, height = 4)



# print session info
sessionInfo()
