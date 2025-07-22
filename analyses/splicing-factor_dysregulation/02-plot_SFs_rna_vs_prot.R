################################################################################
# 02-plot_SFs_rna_vs_prot.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# This script usesCPTAC data to generate and plot heatmap of RNA vs proteomics
# of select splicing factors that was previously identified in 
# 01-volcano_plot_mRNA.R script
#
# usage: Rscript 02-plot_SFs_rna_vs_prot.R
################################################################################

## load packages
suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(circlize)
  library(ComplexHeatmap)
  library(msigdbr)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")
results_dir   <- file.path(analysis_dir, "results")

## check and create plots dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# get histology info
hist_file <- file.path(root_dir, "data", "histologies-plot-group.tsv")
cluster_file <- file.path(root_dir, "analyses", 
                          "sample-psi-clustering", "results", 
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

hist_df <- read_tsv(cluster_file) %>%
  rename(Kids_First_Biospecimen_ID = sample_id) %>%
  select(Kids_First_Biospecimen_ID, cluster) %>%
  inner_join(read_tsv(hist_file))

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTAC3-pbt.xls") 
hgg_de_file <- file.path(results_dir, "all_hgg-diffSFs_sig_genes.txt")
cluster6_de_file <- file.path(results_dir, "cluster6-diffSFs_sig_genes.txt")
sf_file <- file.path(input_dir, "splicing_factors.txt")
hugo_file <- file.path(root_dir, "analyses", "oncoprint", "input", "hgnc-symbol-check.csv")

# Load datasets
hgg_de_genes <- read_tsv(hgg_de_file) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  pull(gene)

cluster6_de_genes <- read_tsv(cluster6_de_file) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  pull(gene)

sf_genes <- read_lines(sf_file)

hugo_genes <- read_csv(hugo_file, skip = 1) %>%
  pull(`Approved symbol`)

kegg_splice <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::filter(gs_name == "KEGG_SPLICEOSOME") %>%
  pull(gene_symbol) %>%
  unique()

de_genes_list <- list("all_hgg" = hgg_de_genes, "cluster6" = cluster6_de_genes, "sf" = sf_genes, "hugo" = hugo_genes, "kegg_splice" = kegg_splice)
names <- names(de_genes_list)

for (each in names) {
  cptac_data <- readxl::read_excel(cptac_output_file) %>%
    # select only rows with CLK1 or SRFs, remove muts
    filter(grepl("\\srna|\\spro", idx), 
                 !`Data type` %in% c("mut", "cnv")) %>%
    # remove extra info from cols
    rename_with(~ gsub("X7316.", "7316-", .), everything()) %>%
    dplyr::rename(Assay = `Data type`) %>%
    filter(`Gene symbol` %in% de_genes_list[[each]]) %>%
    # clean up naming for plotting
    mutate(Assay = case_when(Assay == "proteo" ~ "Whole Cell Proteomics",
                             Assay == "rna" ~ "RNA-Seq"),
           Assay = fct_relevel(Assay, c("RNA-Seq", "Whole Cell Proteomics")),
           # create new display name for phospho proteins to include phos site
           display_name = case_when(
                                    # add a space after the gene to trick into thinking the rownames are not duplicated
                                    Assay == "Whole Cell Proteomics" ~ paste(`Gene symbol`, " ", sep = " "),
                                    TRUE ~ `Gene symbol`)
    ) %>%
    dplyr::select(display_name, Assay, starts_with("7316")) %>%
    # remove NAs
    select_if(~ !any(is.na(.)))

  # preserve gene names for rownames
  rownames <- cptac_data$display_name

  
  # convert to matrix and then add rownames
  mat <- cptac_data %>%
    dplyr::select(any_of(hist_df$sample_id)) %>%
    #dplyr::select(3:(ncol(cptac_data))) %>%
    as.matrix()
  storage.mode(mat) <- "numeric"
  class(mat)


  # set rownames, convert to matrix
  rownames(mat) <- rownames

  # select row annotations
  row_annot <- cptac_data %>%
    dplyr::select(Assay) %>%
    as.data.frame()

  
  # create anno colors
  anno_col <- list(Assay = c("RNA-Seq" = "#DC3220", "Whole Cell Proteomics" = "#40B0A6"))

  # Heatmap annotation
  row_anno = rowAnnotation(df = row_annot,
                           col = anno_col, 
                           show_legend = TRUE, 
                           show_annotation_name = FALSE)
  
  # get diagnoses
  hist_data <- hist_df %>%
    filter(sample_id %in% colnames(mat)) %>%
    distinct(sample_id, cluster, plot_group) %>%
    rename(Histology = plot_group,
           Cluster = cluster) %>%
    column_to_rownames("sample_id") %>% 
    select(Histology, Cluster) %>%
    .[colnames(mat), , drop = FALSE]
  
  # histology colors
  hist_col <- list("Histology" = 
                c("Ependymoma" =                       "#2200ff",       
                  "Atypical Teratoid Rhabdoid Tumor" = "#4d0d85",       
                  "Other high-grade glioma" =          "#ffccf5",       
                  "Low-grade glioma" =                 "#8f8fbf",       
                  "Meningioma" =                       "#2db398",       
                  "DIPG or DMG" =                      "#ff40d9",       
                  "Medulloblastoma" =                  "#a340ff",       
                  "Other tumor" =                      '#b5b5b5',       
                  "Mesenchymal tumor" =                "#7fbf00",       
                  "Craniopharyngioma" =                "#b2502d",       
                  "Mixed neuronal-glial tumor" =       "#685815",       
                  "Non-neoplastic tumor" =             "#FFF5EB",       
                  "Choroid plexus tumor" =             "#00441B",       
                  "Schwannoma" =                       "#ab7200",       
                  "Neurofibroma plexiform" =           "#e6ac39",       
                  "Other CNS embryonal tumor" =        "#b08ccf",       
                  "Germ cell tumor" =                  "#0074d9"),
                "Cluster" = c("1" = "#B2DF8A",
                              "2" = "#E31A1C",
                              "3" = "#33A02C",
                              "4" = "#A6CEE3",
                              "5" = "#FB9A99",
                              "6" = "#FDBF6F",
                              "7" = "#CAB2D6",
                              "8" = "#FFFF99",
                              "9" = "#1F78B4",
                              "10" = "#B15928",
                              "11" = "#6A3D9A"))
  
  
  column_anno = columnAnnotation(df = hist_data,
                                 col = hist_col,
                                 show_legend = TRUE, 
                                 show_annotation_name = FALSE)

  # Make heatmap without legends
  heat_plot <- Heatmap(mat,
                       name = "Z-score",
                       col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                       cluster_rows = FALSE,
                       row_split = row_annot$Assay, 
                       column_gap = 0.5,
                       show_row_names = TRUE,
                       show_column_names = FALSE,
                       show_heatmap_legend = TRUE,
                       cluster_columns = TRUE, 
                       top_annotation = column_anno,
                       right_annotation = row_anno,
                       row_title = NULL, 
                       column_title = NULL, 
                       row_names_gp = grid::gpar(fontsize = 9),
                       column_title_side = "top",
                       heatmap_legend_param = list(legend_direction = "horizontal", 
                                                   legend_position = "top"))

  heatmap_output_file <- file.path(plots_dir, paste0(each,"-SF_RNA_vs_protein_levels_heatmap.pdf"))
  pdf(heatmap_output_file, width = 10, height = 7)
  draw(heat_plot, heatmap_legend_side = "top")
  dev.off()
}

