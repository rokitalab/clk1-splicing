# function to identify differential pathways using GSVA
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GSVA)
  library(pheatmap)
  library(Biobase)
})


root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# Set file paths
expr_file <- file.path(data_dir,
                       "gene-expression-rsem-tpm-collapsed.rds")

stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")
polyA_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-poly-A_stranded.tsv")

geneSet_file <- file.path(results_dir, "hallmark_kegg_splice_geneset_mrna.rds")

# Wrangle data
expr <- readRDS(expr_file)

pathways <- readRDS(geneSet_file)

# define list of cluster file paths
cluster_files <- c(stranded_cluster_file,
                   polyA_cluster_file)
names(cluster_files) <- c("stranded", "poly-A-stranded")

# loop through cluster files to run GSVA and plot differentially expressed pathways
for (libtype in names(cluster_files)){
  
  # read cluster df
  cluster_df <- read_tsv(cluster_files[libtype])

  # match up
  expr_mat <- expr %>%
    as.data.frame() %>%
    dplyr::select(cluster_df$sample_id)
  
  cluster_output <- cluster_df %>%
    mutate(tmp = sample_id) %>%
    column_to_rownames('tmp')
  
  # convert into an eset with cluster info
  eset <- Biobase::ExpressionSet(assayData = as.matrix(expr_mat), 
                                 phenoData = new("AnnotatedDataFrame", data = cluster_output))
    
  # GSVA analysis
  gsea_scores_param <- gsvaParam(eset,
                                 geneSets = pathways,
                                 kcdf = "Gaussian",
                                 assay = NA_character_,
                                 annotation = NA_character_,
                                 tau = 1,
                                 minSize = 1, 
                                 maxSize = 1500, ## Arguments from K. Rathi
                                 maxDiff = TRUE) ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
    gsva_eset <- gsva(gsea_scores_param, verbose = TRUE)
    
    # convert to long format
    gsva_scores <- gsva_eset@assayData$exprs %>%
      as.data.frame() %>%
      rownames_to_column("geneset") %>%
      gather("sample_id", "score", -c("geneset")) %>%
      dplyr::select(sample_id, geneset, score)
    
    # save gsva output
    readr::write_tsv(gsva_scores,
                     file = paste0(results_dir, "/gsva_output_", libtype, ".tsv"))
    
    # differential expression across clusters
    
    # per cluster analysis
    n_clusters <- unique(gsva_eset[["cluster"]])
    all_pathways <- c()
    
  for(i in 1:length(n_clusters)){
      
      # cluster of interest vs others
      cluster_of_interest <- n_clusters[i]
      type <- ifelse(gsva_eset[["cluster"]] == cluster_of_interest, "ref", "others")
      
      # create design and fit
      mod <- model.matrix(~ factor(type))
      colnames(mod) <- c("others", "ref_vs_others")
      rownames(mod) <- gsva_eset[["sample_id"]]
      fit <- limma::lmFit(object = gsva_eset, design = mod)
      fit <- limma::eBayes(fit)
      toptable_output <- limma::topTable(fit, coef = 2, n = Inf)
      toptable_output <- toptable_output %>%
        rownames_to_column("Geneset")
      
      # save filtered output
      out_file <- paste0("cluster_", i, "_pathway_", libtype, ".tsv")
      readr::write_tsv(x = toptable_output %>% 
                         filter(P.Value < 0.05), file = file.path(results_dir, out_file))
      
      #  pull N most significant pathways only
      DEpwys <- toptable_output %>%
        filter(adj.P.Val < 0.05) %>%
        arrange(adj.P.Val) %>%
        head(5) %>%
        pull(Geneset) 
      all_pathways <- c(all_pathways, DEpwys)
    }
    
      # plot top 5 pathways per cluster
    all_pathways <- unique(all_pathways)
  
    DEpwys_es <- Biobase::exprs(gsva_eset[all_pathways, ])
    rownames(DEpwys_es) <- gsub("KEGG_|HALLMARK_", "", rownames(DEpwys_es))
    rownames(DEpwys_es) <- gsub("_", " ", rownames(DEpwys_es))
    DEpwys_annot <- pData(gsva_eset[DEpwys,])
    DEpwys_annot['sample_id'] <- NULL
    DEpwys_annot$cluster <- as.character(DEpwys_annot$cluster)
    
    DEpwys_annot <- DEpwys_annot %>%
      rownames_to_column("Kids_First_Biospecimen_ID") %>%
      dplyr::select(Kids_First_Biospecimen_ID, cluster, plot_group, plot_group_hex) %>%
      column_to_rownames("Kids_First_Biospecimen_ID")
    
    # rename annotation columns
    DEpwys_annot <- DEpwys_annot %>%
      dplyr::rename("Histology" = "plot_group",
                    "Cluster" = "cluster")
    
    # colors for marking different clusters - to match first clustering heatmap
    thisPal <- c("#B2DF8A","#E31A1C","#33A02C","#A6CEE3","#FB9A99","#FDBF6F",
                 "#CAB2D6","#FFFF99","#1F78B4","#B15928","#6A3D9A","#FF7F00",
                 "#2ef4ca","#f4cced","#bd18ea","#f4cc03","#05188a","#e5a25a",
                 "#06f106", #bright green
                 "#85848f", #med gray
                 "#000000", #black
                 "#076f25", #dark green
                 "#93cd7f",#lime green
                 "#4d0776" #dark purple
    )
    
    l <- thisPal[1:length(unique(cluster_df$cluster))]
    names(l) <- as.character(1:length(unique(cluster_df$cluster)))
    mycolors <- list()
    mycolors[['Cluster']] <- l 
    
    # create annotation for plot group
    short_histology_palettes <- DEpwys_annot %>%
      dplyr::select(Histology, plot_group_hex) %>%
      unique()
    mycolors[['Histology']] <- short_histology_palettes$plot_group_hex
    names(mycolors[['Histology']]) <- short_histology_palettes$Histology
    
    # remove colors from annotation table
    DEpwys_annot$plot_group_hex <- NULL
    color_palette <- colorRampPalette(c("blue", "white", "darkorange"))
    heat_colors <- color_palette(6)
    
    # sort cluster df and differential expression matrix by cluster and histology
    cluster_df <- cluster_df %>%
      dplyr::arrange(cluster, plot_group)
    
    DEpwys_es <- DEpwys_es[,cluster_df$sample_id]
  
    # plot pathway heatmap
    pheatmap::pheatmap(DEpwys_es, scale = "row", 
                       fontsize = 8,
                       treeheight_row = 20, 
                       show_colnames = F, 
                       annotation = DEpwys_annot, 
                       annotation_colors = mycolors, 
                       cluster_cols = FALSE, 
                       color = heat_colors,
                       name = "GSVA score",
                       filename = file.path(paste0(plot_dir, "/top", 5, "_pathways_", libtype, ".pdf")), 
                       width = 12, height = 5.5)
  
}
