################################################################################
# 02-SF-mutations.R
# written by Jo Lynne Rokita
#
# usage: Rscript 02-SF-mutations.R
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(maftools)
  library(vroom)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(gtools)
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "oncoprint")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
results_dir   <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")


## check and create plots/results dir
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

source(file.path(input_dir, "mutation-colors.R"))

## input files
cons_maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
tumor_only_maf_file <- file.path(data_dir,"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz")
clin_file <- file.path(data_dir, "histologies-plot-group.tsv")
sf_file <- file.path(root_dir, "analyses","splicing-factor_dysregulation/input","splicing_factors.txt")
hugo_file <- file.path(input_dir, "hgnc-symbol-check.csv")
cluster_file <- file.path(root_dir, "analyses", 
                          "sample-psi-clustering", "results", 
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

# read in files
histologies_df <- read_tsv(clin_file, guess_max = 100000)

## get splicing factor list + CLKs and SRPKs
hugo_genes <- read_csv(hugo_file, skip = 1) %>%
  pull(`Approved symbol`)
sf_genes <- readLines(sf_file) %>%
  unique()

gene_list <- list("hugo" = hugo_genes, "sf" = sf_genes)
gene_list_names <- names(gene_list)

cluster_df <- read_tsv(cluster_file) %>%
  rename(Kids_First_Biospecimen_ID = sample_id)

cluster7_rna <- cluster_df %>%
  filter(cluster == 7) %>%
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id")])

cluster7_dna <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         match_id %in% cluster7_rna$match_id) %>%
  pull(Kids_First_Biospecimen_ID)

matched_dna_samples <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel")) %>%
  pull(Kids_First_Biospecimen_ID)

sample_list <- list("cluster7" = cluster7_dna, "all" = matched_dna_samples)
sample_list_names <- names(sample_list)


# maf cols to select
maf_cols <- c("Hugo_Symbol", 
              "Chromosome", 
              "Start_Position", 
              "End_Position", 
              "HGVSp_Short",
              "Reference_Allele", 
              "Tumor_Seq_Allele2", 
              "Variant_Classification", 
              "Variant_Type",
              "Tumor_Sample_Barcode",
              "t_ref_count",
              "t_alt_count",
              "Transcript_ID",
              "EXON",
              "PolyPhen",
              "SIFT",
              "gnomad_3_1_1_splice_ai_consequence")


# variant classifications to exclude
exclude_class <- c("Silent", "Intron", "RNA", "IGR", "3'Flank", "5'Flank")

# read in and combine MAFs
cons_maf <- data.table::fread(cons_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols)) %>%
  filter(!Variant_Classification %in% exclude_class)

tumor_only_maf <- data.table::fread(tumor_only_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols))  %>%
  filter(!Variant_Classification %in% exclude_class) 

maf <- cons_maf %>%
  bind_rows(tumor_only_maf) %>% 
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

## filter maf for samples with RNA splicing
for (each in gene_list_names) {
  for (sample in sample_list_names) {

  maf_filtered <- maf %>%
    dplyr::filter(Hugo_Symbol %in% gene_list[[each]],
                  Tumor_Sample_Barcode  %in% sample_list[[sample]],
                  Variant_Classification %in% names(colors)) %>%
    dplyr::mutate(pred_deleterious = case_when((grepl("dam", PolyPhen) | grepl("deleterious\\(", SIFT)) ~ "yes",
                                   PolyPhen == "" & SIFT == "" ~ "yes",
                                   TRUE ~ "no")) #%>%
    #dplyr::filter(pred_deleterious == "yes")
  
    print(nrow(maf_filtered))
    write_tsv(maf_filtered, file.path(paste0(results_dir, "/", each, "-", sample, "-samples.maf")))
    
  }
}
