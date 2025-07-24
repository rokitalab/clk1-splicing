################################################################################
# 02-SF-mutations.R
# written by Jo Lynne Rokita
#
# usage: Rscript 02-SF-mutations.R
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "oncoprint")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
results_dir   <- file.path(analysis_dir, "results")

## check and create plots/results dir
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## input files
cons_maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
tumor_only_maf_file <- file.path(data_dir,"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz")
clin_file <- file.path(data_dir, "histologies-plot-group.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
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

cluster6_samples <- cluster_df %>%
  filter(cluster == 6) %>%
  pull(Kids_First_Biospecimen_ID)

indep_rna_df <- read_tsv(indep_rna_file) %>% 
  dplyr::filter(cohort == 'PBTA') %>%
  # get match id for DNA samples
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id", "plot_group", "RNA_library")]) %>%
  filter(RNA_library == "stranded") %>%
  left_join(cluster_df %>% dplyr::select(Kids_First_Biospecimen_ID, cluster))

matched_dna_samples <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         is.na(RNA_library),
         !is.na(pathology_diagnosis),
         match_id %in% indep_rna_df$match_id) %>%
  pull(Kids_First_Biospecimen_ID)

sample_list <- list("cluster6" = cluster6_samples, "all" = matched_dna_samples)
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

# read in and combine MAFs
cons_maf <- data.table::fread(cons_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols)) 

tumor_only_maf <- data.table::fread(tumor_only_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols)) 

maf <- cons_maf %>%
  bind_rows(tumor_only_maf) %>% 
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

## filter maf for samples with RNA splicing + HGGs
for (each in gene_list_names) {
  for (sample in sample_list_names) {

  maf_filtered <- maf %>%
    dplyr::filter(Hugo_Symbol %in% gene_list[[each]],
                  Tumor_Sample_Barcode  %in% sample_list[[sample]]) %>%
    dplyr::mutate(keep = case_when(Variant_Classification == "Missense_Mutation" & (grepl("dam", PolyPhen) | grepl("deleterious\\(", SIFT)) ~ "yes",
                                   Variant_Classification == "Missense_Mutation" & PolyPhen == "" & SIFT == "" ~ "yes",
                                   Variant_Classification != "Missense_Mutation" ~ "yes",
                                   TRUE ~ "no")) %>%
    dplyr::filter(keep == "yes")
    print(nrow(maf_filtered))
    write_tsv(maf_filtered, file.path(paste0(results_dir, "/", each, "-", sample, "-samples.maf")))
  }
}

