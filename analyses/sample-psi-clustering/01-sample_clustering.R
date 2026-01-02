# R Corbett 2025
#
# Peform sample clustering of primary tumors based on splice event PSIs

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cluster)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(analysis_dir, 
                 "util", "heatmap_function.R"))

# set file paths

psi_mat_file <- file.path(results_dir,
                          "pbta-splice-event-psis.RDS")

histologies_file <- file.path(root_dir, "analyses", "cohort_summary", "results", "histologies-plot-group.tsv")

# wrangle data

psi_mat <- readRDS(psi_mat_file)

histologies <- read_tsv(histologies_file) %>%
  dplyr::mutate(molecular_subtype = case_when(
    short_histology == "Oligodendroglioma" ~ "OG",
    TRUE ~ molecular_subtype
  ))

plotgroup_palette <- unique(histologies$plot_group_hex)
names(plotgroup_palette) <- unique(histologies$plot_group)

# transpose psi matrix
psi_mat <- psi_mat %>%
  column_to_rownames("sample_id") %>%
  t()

# define vector of n variable events to test for clustering
n_events <- 5000

# extract stranded and poly-A stranded libraries to cluster separately
stranded_samples <- histologies %>%
  dplyr::filter(RNA_library == "stranded") %>%
  pull(Kids_First_Biospecimen_ID)

library_type <- list("stranded" = stranded_samples)

pdf(NULL)

# perform hierarchical clustering on different library groups and number of most variable splice events

# loop through library strategies
for (library in names(library_type)){
  
  # loop through number of events
  for (n in n_events){
    
    # filter for samples and splice events reported in at least 25% of samples
    cluster_mat <- psi_mat[,colnames(psi_mat) %in% library_type[[library]]]
    cluster_mat <- cluster_mat[rowSums(!is.na(cluster_mat)) > ncol(cluster_mat)*0.25,]

    # pull splice IDs with the highest variance
    psi_vars <- apply(cluster_mat, 1, var, na.rm = TRUE)
    var_splice_ids <- names(sort(psi_vars, decreasing = TRUE))[1:n]
    
    # filter matrix for most variable splice events
    cluster_mat <- cluster_mat[var_splice_ids,]
    
    # construct distance matrix using Euclidean distance measure
    dist_matrix <- as.matrix(dist(t(cluster_mat), method = "euclidean", diag = TRUE, upper = TRUE))
    
    dist_matrix <- as.dist(dist_matrix)
    
    # perform clustering using ward.D2 method
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    # choose optimal number of clusters, which has been previously assessed and are defined here for each library type and 
    clusters <- cutree(hc, k = 10)
    
    # create df of cluster assignment by sample
    cluster_assignment_df <- data.frame(sample_id = names(clusters),
                                        cluster = clusters)
    
    # create df with sample metadata and cluster assignment
    cluster_df <- data.frame(sample_id = colnames(cluster_mat)) %>%
      left_join(histologies %>% dplyr::select(Kids_First_Biospecimen_ID,
                                              plot_group,
                                              plot_group_hex,
                                              RNA_library,
                                              molecular_subtype),
                by = c("sample_id" = "Kids_First_Biospecimen_ID")) %>%
      left_join(cluster_assignment_df) %>%
      dplyr::mutate(cluster = factor(cluster,
                                     levels = sort(unique(cluster)))) %>%
      dplyr::arrange(cluster, plot_group)
    
    hist_ct <- cluster_df %>%
      dplyr::count(plot_group) %>%
      dplyr::rename(group_n = n)
    
    cluster_df <- cluster_df %>%
      left_join(hist_ct) %>%
      dplyr::mutate(plot_group_n = glue::glue("{plot_group} (n = {group_n})")) %>%
      dplyr::mutate(plot_group_n = str_replace(plot_group_n, "Atypical teratoid rhabdoid tumor",
                                               "ATRT"))
    
    # define plot group palette
    plotgroup_palette <- unique(cluster_df$plot_group_hex)
    names(plotgroup_palette) <- unique(cluster_df$plot_group_n)
    
    # define colors for clusters
    cluster_cols <- c("#B2DF8A","#E31A1C","#33A02C","#A6CEE3","#FB9A99","#FDBF6F",
                      "#CAB2D6","#FFFF99","#1F78B4","#B15928","#6A3D9A","#FF7F00",
                      "#2ef4ca","#f4cced","#bd18ea")
    names(cluster_cols) <- 1:length(cluster_cols)
    cluster_cols <- cluster_cols[1:length(unique(clusters))]
    
    # define sample annotation 
    ha = HeatmapAnnotation(
      Histology = cluster_df$plot_group_n,
      Cluster = cluster_df$cluster,
      `RNA library type` = cluster_df$RNA_library,
      col = list("Histology" = plotgroup_palette,
                 Cluster = cluster_cols,
                 `RNA library type` = c("exome_capture" = "black",
                                        "poly-A" = "red2",
                                        "poly-A stranded" = "yellow2",
                                        "stranded" = "green3")),
      annotation_name_gp= gpar(fontsize = 10))
    
    # define psi color scale
    col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    
    # reorder matrix to match cluster df
    cluster_mat <- cluster_mat[,cluster_df$sample_id]
    
    # plot cluster heatmap
    ht <- Heatmap(cluster_mat,
                  cluster_columns = TRUE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "ward.D2",
                  cluster_rows = TRUE,
                  clustering_distance_rows = "euclidean",
                  clustering_method_rows = "ward.D2",
                  col = col_fun,
                  na_col = "gray",
                  name = "PSI",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  top_annotation = ha,
                  use_raster = FALSE,
                  heatmap_legend_param = list(legend_gp = gpar(fontsize = 10),
                                              labels_gp = gpar(fontsize = 10)))
    
    ht_drawn <- draw(ht)
    col_order <- column_order(ht_drawn)
    ordered_colnames <- colnames(cluster_mat)[col_order] %>%
      as.data.frame() %>%
      write_tsv(file.path(results_dir, glue::glue("colorder-top-{n}-events-{library}.txt")))
    
    pdf(NULL)
    
    pdf(file.path(plot_dir,
                  glue::glue("sample-psi-heatmap-top-{n}-events-{library}.pdf")),
        width = 7, height = 7)
    
    print(ht)
    
    dev.off()
    
    
    
    # plot enrichment of tumor histologies within clusters
    histology_enr <- plot_enr(cluster_df, 
                              "plot_group", "cluster",
                              sort(unique(cluster_df$plot_group)),
                              unique(cluster_df$cluster),
                              padjust = TRUE)
    
    pdf(file.path(plot_dir, 
                  glue::glue("sample-cluster-histology-enr-top-{n}-events-{library}.pdf")),
        width = 11, height = 8)
    
    print(histology_enr)
    
    dev.off()
    
    # For HGG/DMG, LGG, MB, and ATRT, plot subtype distribution across identified clusters
    
    # HGG
    hgg_df <- cluster_df %>%
      dplyr::filter(plot_group %in% c("Other high-grade glioma",
                                      "Diffuse midline glioma"),
                    !is.na(molecular_subtype),
                    !grepl("To be classified", molecular_subtype)) %>%
      dplyr::mutate(molecular_subtype = str_replace(molecular_subtype, ", TP53", "")) %>%
      dplyr::mutate(mol_sub_group = case_when(
        grepl("IHG", molecular_subtype) ~ "IHG",
        grepl("DIPG", molecular_subtype) ~ str_replace(molecular_subtype, "DIPG", "DMG"),
        TRUE ~ molecular_subtype
      ))
    
    hgg_enr <- plot_enr(hgg_df, 
                        "mol_sub_group", "cluster",
                        sort(unique(hgg_df$mol_sub_group)),
                        unique(hgg_df$cluster),
                        padjust = TRUE)
    
    pdf(file.path(plot_dir, 
                  glue::glue("hgg-dmg-sample-cluster-subtype-enr-top-{n}-events-{library}.pdf")),
        width = 8, height = 6)
    
    print(hgg_enr)
    
    dev.off()
    
    # LGG
    lgg_df <- cluster_df %>%
      dplyr::filter(plot_group %in% c("Low-grade glioma"),
                    !is.na(molecular_subtype),
                    !grepl("To be classified", molecular_subtype)) %>%
      dplyr::mutate(mol_sub_group = case_when(
        grepl("V600E", molecular_subtype) ~ "LGG, BRAF V600E",
        grepl("-BRAF", molecular_subtype) ~ "LGG, BRAF fusion",
        grepl("SEGA", molecular_subtype) ~ "SEGA",
        grepl("wildtype", molecular_subtype) ~ "LGG, wildtype",
        !grepl("classified", molecular_subtype) ~ "LGG, Other alteration",
        TRUE ~ molecular_subtype
      ))
    
    lgg_enr <- plot_enr(lgg_df, 
                        "mol_sub_group", "cluster",
                        sort(unique(lgg_df$mol_sub_group)),
                        unique(lgg_df$cluster),
                        padjust = TRUE)
    
    pdf(file.path(plot_dir, 
                  glue::glue("lgg-sample-cluster-subtype-enr-top-{n}-events-{library}.pdf")),
        width = 7, height = 4)
    
    print(lgg_enr)
    
    dev.off()
    
    # Medullo
    mb_df <- cluster_df %>%
      dplyr::filter(plot_group %in% c("Medulloblastoma"),
                    !is.na(molecular_subtype))
    
    mb_enr <- plot_enr(mb_df, 
                       "molecular_subtype", "cluster",
                       sort(unique(mb_df$molecular_subtype)),
                       unique(mb_df$cluster),
                       padjust = TRUE)
    
    pdf(file.path(plot_dir, 
                  glue::glue("mb-sample-cluster-subtype-enr-top-{n}-events-{library}.pdf")),
        width = 5, height = 3)
    
    print(mb_enr)
    
    dev.off()
    
    # ATRT
    atrt_df <- cluster_df %>%
      dplyr::filter(plot_group %in% c("Atypical teratoid rhabdoid tumor"),
                    !is.na(molecular_subtype),
                    !grepl("To be classified", molecular_subtype))
    
    if (library != "poly-A_stranded"){
      
      atrt_enr <- plot_enr(atrt_df, 
                           "molecular_subtype", "cluster",
                           sort(unique(atrt_df$molecular_subtype)),
                           unique(atrt_df$cluster),
                           padjust = TRUE)
      
      pdf(file.path(plot_dir, 
                    glue::glue("atrt-sample-cluster-subtype-enr-top-{n}-events-{library}.pdf")),
          width = 6, height = 3)
      
      print(atrt_enr)
      
      dev.off()
      
    }
    
    # save psi matrix 
    
    saveRDS(cluster_mat,
            file.path(results_dir,
                      glue::glue("psi-matrix-top-{n}-events-{library}.rds")))
    write_tsv(cluster_df,
              file.path(results_dir,
                        glue::glue("sample-cluster-metadata-top-{n}-events-{library}.tsv")))
    
  }
  
}

# print session info
sessionInfo()