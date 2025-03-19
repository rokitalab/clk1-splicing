################################################################################
# 03-plot-sbi_with_cluster-mem.R
# script that takes in "splicing_index.total.txt" data file from previous module
# and intersects with cluster membership info from current module output file 
# and plots number of samples per cluster stratified by splicing burden
# 
# written by Ammar Naqvi, Jo Lynne Rokita, Ryan Corbett
#
# usage: Rscript 03-plot-sbi_with_cluster-mem.R
################################################################################
library("tidyverse")

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-psi-clustering")

plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## filepaths 
stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")
polyA_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-poly-A stranded.tsv")

sbi_file <- file.path(root_dir, "analyses", "clustering_analysis",
                      "input", "splicing_index.total.txt")

# Wrangle data
sbi_df <- read_tsv(sbi_file)

# define list of cluster file paths
cluster_files <- c(stranded_cluster_file,
                   polyA_cluster_file)
names(cluster_files) <- c("stranded", "poly-A-stranded")

# loop through cluster files to plot SBI group distributions
for (type in names(cluster_files)){

  # read cluster df
  cluster_df <- read_tsv(cluster_files[type])

  # define lower and upper quartile values
  SI_total_high     <- quantile(sbi_df$SI, probs=.75, names=FALSE)
  SI_total_low      <- quantile(sbi_df$SI, probs=.25, names=FALSE)

  ## add column with "High or "Low" for SBI info
  sbi_outliers <- sbi_df %>%
    filter(sbi_df$SI < SI_total_low | sbi_df$SI > SI_total_high) %>% 
    mutate(level = case_when(SI < SI_total_low ~ "Low", 
                           SI > SI_total_high  ~ "High" )) %>% 
    inner_join(cluster_df %>% 
                 dplyr::select(sample_id, cluster), 
               by = c("Sample" = "sample_id"))
  
  summ_sbi_cl <- sbi_outliers %>% group_by(cluster, level) %>%
    summarise(n = n()) %>%
    mutate(Freq = n/sum(n)) %>%
    ungroup()
  
  # are there any clusters without high or low SBI members? YES
  clusters <- 1:length(unique(cluster_df$cluster))
  
  missing_clusters <- setdiff(clusters, unique(as.character(summ_sbi_cl$cluster)))
  
  # Ensure missing_cluster_df has the correct structure with the specified columns
  missing_cluster_df <- data.frame(cluster = integer(),
                                   level = character(),
                                   n = integer(),
                                   Freq = numeric())
  
  # Loop through each missing cluster to add a new row for it
  for(missing in missing_clusters) {
    missing_cluster_df <- rbind(missing_cluster_df, 
                                data.frame(cluster = missing,
                                           level = "High",
                                           n = NA,
                                           Freq = NA))
  }
  
  summ_sbi_cl <- summ_sbi_cl %>%
    bind_rows(missing_cluster_df)
  
  # plot
  sbi_vs_cl_barplot <- ggplot(summ_sbi_cl, aes(x = as.factor(cluster), y = n, fill = level)) + 
    geom_bar(stat = "identity") +
    labs(x="Cluster", y="Number of Samples") +
    theme_Publication() + 
    scale_fill_manual(values= c("Low" = "#FFC20A", "High"="#0C7BDC"), name = "SBI Level")
  
  pdf(paste0(plot_dir, "/cluster_membership_sbi_group_", type, ".pdf"), 
      height = 3, width = 8)
  
  print(sbi_vs_cl_barplot)
  
  dev.off()
  
}

# session info
sessionInfo()

