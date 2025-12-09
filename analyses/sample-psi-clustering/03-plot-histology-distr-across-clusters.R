################################################################################
# 03-plot-histology-distr-across-clusters.R
# Plot the distribution of histologies across clusters
#
# Author: Ammar Naqvi, Jo Lynne Rokita, Ryan Corbett
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

## filepaths 
stranded_cluster_file <- file.path(results_dir, "sample-cluster-metadata-top-5000-events-stranded.tsv")

subtype_hex <- file.path(input_dir, "subtype_hex.tsv")

# Wrangle data
histologies <- read_tsv(file.path(root_dir, "analyses",
                                  "cohort_summary", "results",
                                  "histologies-plot-group.tsv"))
subtype_hex_codes <- read_tsv(subtype_hex)

# define list of cluster file paths
cluster_files <- c(stranded_cluster_file)
names(cluster_files) <- c("stranded")

# loop through cluster files to plot histology and subtype distributions
for (type in names(cluster_files)){
  
  pdf(NULL)

  # read in palette and cluster file
  cluster_df <- read_tsv(cluster_files[type]) %>%
    dplyr::rename(Kids_First_Biospecimen_ID = sample_id)
  
  histologies_df <- cluster_df %>%
    mutate(molecular_subtype_display = case_when(grepl("KIAA", molecular_subtype) ~ "KIAA1549::BRAF",
                                                 grepl("LGG, BRAF V600|LGG, RTK, BRAF V600E", molecular_subtype) ~ "BRAF V600E",
                                                 grepl("LGG, IDH|LGG, other MAPK, IDH", molecular_subtype) ~ "IDH",
                                                 molecular_subtype %in% c("LGG, wildtype") ~ "LGG, WT",
                                                 grepl("SEGA", molecular_subtype) ~ "SEGA",
                                                 grepl("To be classified", molecular_subtype) ~ "To be classified",
                                                 grepl("LGG", molecular_subtype) & grepl("other MAPK", molecular_subtype) ~ "Other MAPK",
                                                 grepl("LGG", molecular_subtype) ~ "Other alteration",
                                                 grepl("HGG, IDH", molecular_subtype) ~ "IDH",
                                                 grepl("HGG, H3 wildtype|DIPG, H3 wildtype", molecular_subtype) ~ "H3 wildtype",
                                                 grepl("HGG, PXA", molecular_subtype) ~ "PXA",
                                                 grepl("IHG", molecular_subtype) ~ "IHG",
                                                 grepl("H3 G35", molecular_subtype) ~ "H3 G35",
                                                 grepl("H3 K28", molecular_subtype) ~ "H3 K28",
                                                 is.na(molecular_subtype) & plot_group == "Other high-grade glioma" ~ "Oligodendroglioma",
                                                 TRUE ~ molecular_subtype)) %>%
    left_join(histologies %>% dplyr::select(plot_group, plot_group_hex)) %>%
    distinct()
  
  histologies_df <- histologies_df %>%
    mutate(cluster = as.numeric(cluster)) %>%
    add_count(cluster, name = "group_n") %>%
    mutate(
      cluster_label = paste0(cluster, " (n = ", group_n, ")"),
      cluster_label = factor(cluster_label, levels = unique(cluster_label[order(cluster)]))
    )
  
  color_df <- histologies_df %>%
    dplyr::select(plot_group_hex, plot_group) %>%
    dplyr::filter(!is.na(plot_group)) %>%
    unique()
  cols <- as.character(color_df$plot_group_hex)
  names(cols) <- as.character(color_df$plot_group)
  
  # create plot
  pdf(file.path(paste0(plots_dir, "/cluster_membership_", type, ".pdf")),
      height = 7, width = 5)
  
  hist_plot <- ggplot(histologies_df, aes(fill=plot_group, x=cluster_label)) +
                geom_bar(position="fill") + 
                xlab("Cluster") + ylab("Fraction") +
                scale_fill_manual("Histology", values = cols) + 
                theme_Publication() +
                theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(hist_plot)
  
  dev.off()
  
  ## stratify by molecular subtype
  # subset to only those plot groups: ATRT, MB, LGG, HGG, EPN
  hist_subset <- histologies_df %>%
    filter(plot_group %in% c("Medulloblastoma", 
                             "Atypical teratoid rhabdoid tumor",
                             "Low-grade glioma",
                             "Other high-grade glioma",
                             "Diffuse midline glioma")) %>%
    mutate(plot_group = case_when(plot_group %in% c("Other high-grade glioma",
                                                    "Diffuse midline glioma") ~ "High-grade glioma",
                                  TRUE ~ plot_group))
  
  cols <- subtype_hex_codes$hex_code
  names(cols) <- subtype_hex_codes$molecular_subtype_display
  
  pdf(file.path(paste0(plots_dir, "/cluster_membership-subtypes_", type, ".pdf")),
      height = 6, width = 9.5)
  
  subtype_plot <- hist_subset %>% 
    ggplot(aes(fill=molecular_subtype_display, x= factor(cluster))) +
    geom_bar(stat="count", position="stack", color = "black", size = 0.2) + 
    facet_wrap(~plot_group, nrow = 2, scales = "free_y") +
    xlab("Cluster") + 
    ylab("Frequency") +
    theme_Publication() +
    scale_fill_manual(name = "Molecular Subtype", values = cols) +
    guides(fill = guide_legend(ncol = 1))
  
  print(subtype_plot)
  
  dev.off()
  
  # Load BRAF fusion breakpoint df
  
  braf_fusions <- read_tsv(file.path(input_dir, "lgg-braf-fusion-breakpoint-annotation.tsv"))
  
  # subset cluster df for BRAF fusion positive LGG
  
  lgg_braf_cluster_df <- cluster_df %>%
    dplyr::filter(grepl("-BRAF", molecular_subtype),
                  plot_group == "Low-grade glioma") %>%
    left_join(braf_fusions %>%
                dplyr::select(Kids_First_Biospecimen_ID,
                              breakpoint_group)) %>%
    dplyr::filter(!is.na(breakpoint_group))
  
  # assess distribution of braf fusion breakpoint groups across clusters
  
  braf_fusion_enr <- plot_enr(lgg_braf_cluster_df, 
                      "breakpoint_group", "cluster",
                      sort(unique(lgg_braf_cluster_df$breakpoint_group)),
                      unique(lgg_braf_cluster_df$cluster),
                      padjust = TRUE)
  
  pdf(file.path(plots_dir, 
                glue::glue("lgg-braf-fusion-breakpoint-cluster-enr-top-5000-events-{type}.pdf")),
      width = 5, height = 3)
  
  print(braf_fusion_enr)
  
  dev.off()
  
  # Load MB-SHH molecular subgroup table
  
  mbshh_subtypes <- read_tsv(file.path(input_dir, "mb_shh_molecular_subtypes.tsv")) %>%
    dplyr::rename(shh_subtype = molecular_subtype)
  
  # define molecualr subgroups using methylation subtyping for Group3/4 MB, and SHH subtype table for MB-SHH
  mb_cluster_df <- cluster_df %>%
    dplyr::filter(plot_group == "Medulloblastoma") %>%
    # append match IDs
    left_join(histologies %>% 
                dplyr::select(Kids_First_Biospecimen_ID,
                              "match_id")) %>%
    # append high-confidence methylation subtypes
    left_join(histologies %>% 
                dplyr::select(match_id, 
                              dkfz_v12_methylation_subclass,
                              dkfz_v12_methylation_subclass_score) %>%
                dplyr::filter(dkfz_v12_methylation_subclass_score >= 0.7)) %>%
    # retain unique subtype per ID
    dplyr::arrange(desc(dkfz_v12_methylation_subclass_score)) %>%
    distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
    # append SHH subtypes
    left_join(mbshh_subtypes %>% 
                dplyr::select(match_id, shh_subtype) %>%
                distinct(match_id, .keep_all = TRUE)) %>%
    # define subgroup column that takes SHH subtypes and Group3/4 methylation subtypes
    dplyr::mutate(molecular_subgroup = case_when(
      molecular_subtype == "MB, SHH" ~ shh_subtype,
      TRUE ~ dkfz_v12_methylation_subclass
    )) %>%
    # filter for MB methylation subtypes and remove undefined SHH subtypes
    dplyr::filter(grepl("MB", molecular_subgroup),
                  molecular_subgroup != "MB, SHH")
  
  # plot SHH subtype distribution
  
  mbshh_cluster_df <- mb_cluster_df %>%
    dplyr::filter(molecular_subtype == "MB, SHH")
  
  if (length(unique(mbshh_cluster_df$cluster)) > 1){
  
    mbshh_subgroup_enr <- plot_enr(mbshh_cluster_df, 
                                "molecular_subgroup", "cluster",
                                sort(unique(mbshh_cluster_df$molecular_subgroup)),
                                sort(unique(mbshh_cluster_df$cluster)),
                                padjust = TRUE)
    
    pdf(file.path(plots_dir, 
                  glue::glue("mb-shh-subgroup-cluster-enr-top-5000-events-{type}.pdf")),
        width = 4, height = 3)
    
    print(mbshh_subgroup_enr)
    
    dev.off()
    
  }
  
  # plot Group 3/4 methylation subtype distribution
  
  mbg34_cluster_df <- mb_cluster_df %>%
    dplyr::filter(molecular_subtype %in% c("MB, Group3",
                                           "MB, Group4"))
  
  mbg34_subgroup_enr <- plot_enr(mbg34_cluster_df, 
                                 "molecular_subgroup", "cluster",
                                 sort(unique(mbg34_cluster_df$molecular_subgroup)),
                                 sort(unique(mbg34_cluster_df$cluster)),
                                 padjust = TRUE)
  
  pdf(file.path(plots_dir, 
                glue::glue("mb-group34-subgroup-cluster-enr-top-5000-events-{type}.pdf")),
      width = 3, height = 4)
  
  print(mbg34_subgroup_enr)
  
  dev.off()
  
  # Save BRAF fusion LGG breakpoint group and MB subgroup cluster membership dfs
  
  write_tsv(lgg_braf_cluster_df,
            file.path(results_dir, 
                      glue::glue("lgg-braf-fusion-cluster-membership-{type}.tsv")))
  
  write_tsv(mb_cluster_df,
            file.path(results_dir, 
                      glue::glue("mb-subgroup-cluster-membership-{type}.tsv")))
  
}

# Session info
sessionInfo()

