################################################################################
# 06-get-pfam-domains.R
# Script that annotates alternatively spliced exons with Pfam protein domains using annoFuseData package
# written by Ryan Corbett
#
# usage: Rscript 06-get-pfam-domains.R
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("annoFuseData")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# set file paths
func_site_file <- file.path(results_dir, "kinases-functional_sites.tsv")

# wrangle data
func_site_df <- read_tsv(func_site_file) %>%
  dplyr::rename(hgnc_symbol = gene)

# Load Pfam data
bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))

# merge splicing events with protein domains by gene symbol

splice_pfam_df <- func_site_df %>%
  left_join(bioMartDataPfam, by = "hgnc_symbol") %>%
  # create separate columns for single exon start and end coordinates
  dplyr::mutate(exon_start = as.double(unlist(lapply(strsplit(`Exon Coordinates`, "-"), function(x) x[[1]]))),
                exon_end = as.double(unlist(lapply(strsplit(`Exon Coordinates`, "-"), function(x) x[[2]])))) %>%
  dplyr::select(SpliceID, Differential_freq, mean_dPSI, Uniprot, hgnc_symbol,
                Baseline_preference, exon_start, exon_end, Baseline_freq,
                Average_TPM, pfam_id, chromosome_name,
                strand, NAME, DESC, domain_chr, domain_start, domain_end)

# filter for single exons overlapping functional domains
# Overlap is "Yes" if 1) exon is completely within domain boundaries, 2) domain is completely within exon boundaries, or 3) exon start or end lies within domain boundaries
splice_pfam_df <- splice_pfam_df %>%
  dplyr::mutate(overlaps_domain = case_when(
    (exon_start > domain_start & exon_start < domain_end) | (exon_end > domain_start & exon_end < domain_end) | (domain_start > exon_start & domain_end < exon_end) ~ "Yes",
    TRUE ~ "No"
  )) %>%
  # only retain exons overlapping domains
  dplyr::filter(overlaps_domain == "Yes") %>%
  # calculate overlap length
  dplyr::mutate(overlap_len = case_when(
    # if domain completely within exon, overlap is length of domain
    domain_start > exon_start & domain_end < exon_end ~ domain_end - domain_start,
    # if exon is completely within domain, overlap is length of exon
    exon_start > domain_start & exon_end < domain_end ~ exon_end - exon_start,
    # if domain overlaps downstream portion of exon
    exon_start < domain_start & exon_end > domain_start ~ exon_end - domain_start,
    # if domain overlaps upstream portion of exon
    domain_end > exon_start & domain_end < exon_end ~ domain_end - exon_start
  )) %>%
  # calculate exon length and percent overlap 
  dplyr::mutate(exon_len = exon_end - exon_start) %>%
  dplyr::mutate(overlap_perc = overlap_len/exon_len * 100)

# write to output
write_tsv(splice_pfam_df,
          file.path(results_dir, "splice-events-pfam-annotated.tsv"))

