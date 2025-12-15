################################################################################
# 04-functional-sites-kinases.R 
# script that filters the splicing events overlapping uniprot sites to kinases
# written by Ammar Naqvi, Jo Lynne Rokita, Patricia Sullivan
#
# usage: Rscript 04-functional-sites-kinases.R 
################################################################################

suppressPackageStartupMessages({
  library("vroom")
  library("ggplot2")
  library("tidyverse")
  library("clusterProfiler")
  library("msigdbr")
  library("org.Hs.eg.db")
  library("ggrepel")
  library("cowplot")
  library("ggpubr")
  library("annoFuseData")
  library("stringr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")
input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## output files for final plots
ora_dotplot_path <- file.path(plots_dir,"kinases-ora-plot.pdf")
kinases_functional_sites = file.path(results_dir,"kinases-functional_sites.tsv")
sf_kinases_functional_sites = file.path(results_dir,"splicing-factor-kinases-functional_sites.tsv")

## retrieve psi values from tables
file_psi_pos_func <- file.path(results_dir,"splicing_events.total.pos.intersectunip.ggplot.txt")
file_psi_neg_func <- file.path(results_dir,"splicing_events.total.neg.intersectunip.ggplot.txt")

## histologies
hist_indep_rna_df <- vroom(file.path(data_dir,"histologies-plot-group.tsv"))

## clusters
cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")

## splicing factors 
splicing_factors <- readLines(file.path(root_dir, "analyses", "splicing-factor_dysregulation", "input", "splicing_factors.txt"))

## read table of recurrent functional splicing (skipping)
dpsi_unip_pos <- vroom(file_psi_pos_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  mutate(Preference='Inclusion')

## read table of recurrent functional splicing (inclusion) 
dpsi_unip_neg <- vroom(file_psi_neg_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>% 
  mutate(Preference='Skipping')

psi_comb <- rbind(dpsi_unip_neg,dpsi_unip_pos) %>% 
  mutate(Uniprot = case_when(Uniprot == 'DisulfBond' ~ "Disulfide Bond",
                             Uniprot == 'LocSignal' ~ "Localization Signal",
                             .default = Uniprot),
         Uniprot_wrapped = stringr::str_wrap(Uniprot, width = 10)
         )

# get and filter for kinase genes
known_kinase_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) %>%
  dplyr::rename(gene=Gene_Symbol) %>% 
  dplyr::filter(type=='Kinase')

psi_unip_kinase <- dplyr::inner_join(psi_comb, known_kinase_df, by='gene') 

## over-representation analysis for kinases
# get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")

## filter for kegg pathways that are included in the curated gene sets
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG", "CP:BIOCARTA", "TFT:GTRD"))

## skipping vs incl 
kinase_skip_pref <- psi_unip_kinase %>% 
  dplyr::filter(Preference=='Skipping') 
kinase_incl_pref <- psi_unip_kinase %>% 
  dplyr::filter(Preference=='Inclusion') 

kinase_genes <- c(kinase_incl_pref$gene, kinase_skip_pref$gene)

##kegg gene set
ora_results <- enricher(
  gene = unique(kinase_genes), # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

options(enrichplot.colours = c("darkorange","blue"))
enrich_skip_plot <- enrichplot::dotplot(ora_results,
                                        x = "geneRatio",
                                        size = "Count",
                                        color = "p.adjust",
                                        label_format = 30,
                                        showCategory = 15) +   
  labs(y = "Pathway",
       x = "Gene Ratio") +
  theme_Publication() +
  scale_size(name = "Gene Count") +  
  scale_fill_gradient(low = "darkorange", high = "blue", name = "B-H p-value") +
  guides(
    fill = guide_colorbar(title = "B-H p-value", label.position = "right", barwidth = 1, barheight = 4)
  )

ggplot2::ggsave(filename = ora_dotplot_path,
                plot = enrich_skip_plot, 
                width=8,
                height=7,
                device="pdf")


kinase_skip_pref <- kinase_skip_pref %>%  
  dplyr::mutate('Exon Coordinates' = str_match(SpliceID, "(\\w+)\\:(\\d+\\-\\d+)\\_")[, 3]) %>%
  unique() 

kinase_incl_pref <- kinase_incl_pref %>% 
  dplyr::mutate('Exon Coordinates' = str_match(SpliceID, "(\\w+)\\:(\\d+\\-\\d+)\\_")[, 3]) %>%
  unique()

kinase_pref <- rbind(kinase_skip_pref, kinase_incl_pref) %>% 
  dplyr::select(splicing_case,SpliceID,dPSI,Uniprot, gene, Preference,`Exon Coordinates`)


# Load in differential splicing events entries
total_diff_events_se <- vroom(file.path(results_dir,"splice_events.diff.SE.txt"))
total_diff_events_ri <- vroom(file.path(results_dir,"splice_events.diff.RI.txt"))
total_diff_events_a3ss <- vroom(file.path(results_dir,"splice_events.diff.A3SS.txt"))
total_diff_events_a5ss <- vroom(file.path(results_dir,"splice_events.diff.A5SS.txt"))

total_diff_events <- bind_rows(total_diff_events_se, total_diff_events_ri, total_diff_events_a3ss, total_diff_events_a5ss) %>%
  dplyr::rename(SpliceID="Splice ID") %>%
  dplyr::count(SpliceID)  %>% 
  inner_join(kinase_pref,by="SpliceID") %>%
  dplyr::rename('Frequency'=n) %>%
  filter(Frequency > 1) # Add recurrence filter

# read in exp file
expr_file <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")

exp <- readRDS(expr_file) %>%
  dplyr::select(any_of(hist_indep_rna_df$Kids_First_Biospecimen_ID)) %>%
  mutate(gene = rownames(.)) %>%
  dplyr::filter(gene %in% total_diff_events$gene) %>%
  dplyr::rowwise() %>%  # Ensure you use parentheses here
  dplyr::filter(sum(c_across(where(is.numeric))) >= 1) %>%
  dplyr::ungroup() %>%
  pivot_longer(
    cols = -gene,          # Select all columns to pivot
    names_to = "Sample",          # Column name for former column names
    values_to = "TPM"      # Column name for values
  ) %>%
  group_by(gene) %>%                     # Group by gene
  summarize(Average_TPM = mean(TPM, na.rm = TRUE))  # Compute mean TPM per gene

# Add gene names as a column to count_data
total_diff_events_gene <- inner_join(total_diff_events,exp, by='gene') %>%
  dplyr::rename(Differential_freq=Frequency,
                mean_dPSI=dPSI,
                Baseline_preference=Preference) %>%
  mutate(Baseline_preference = case_when(
    Baseline_preference == "Inclusion" ~ "Skipping",  # Replace "Inclusion" with "Skipping"
    Baseline_preference == "Skipping" ~ "Inclusion",  # Replace "Skipping" with "Inclusion"
    TRUE ~ Baseline_preference                      # Keep other values unchanged
  ))


## write kinase results to table
write_tsv(total_diff_events_gene, kinases_functional_sites)

## subset for splicing factors
total_diff_events_gene_sf <- total_diff_events_gene %>%
  filter(gene %in% splicing_factors)

## write splicing factor kinase results to table
write_tsv(total_diff_events_gene_sf, sf_kinases_functional_sites)