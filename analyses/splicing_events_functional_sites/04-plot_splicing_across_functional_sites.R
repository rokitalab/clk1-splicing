################################################################################
# 04-plot_splicing_across_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 04-plot_splicing_across_functional_sites.R 
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
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_dpsi_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites.pdf")
file_dpsi_kinase_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_kinase.pdf")
ora_dotplot_path <- file.path(plots_dir,"kinases-ora-plot.pdf")
kinases_functional_sites = file.path(results_dir,"kinases-functional_sites.tsv")
sf_kinases_functional_sites = file.path(results_dir,"splicing-factor-kinases-functional_sites.tsv")

## retrieve psi values from tables
file_psi_pos_func <- file.path(results_dir,"splicing_events.SE.total.pos.intersectunip.ggplot.txt")
file_psi_neg_func <- file.path(results_dir,"splicing_events.SE.total.neg.intersectunip.ggplot.txt")

## histologies
hist_indep_rna_df <- vroom(file.path(data_dir,"histologies.tsv"))

## clusters
cluster_file <- file.path(root_dir, "analyses",
                          "sample-psi-clustering", "results",
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")
cluster6_bs_id <- read_tsv(cluster_file) %>%
  filter(cluster == 6) %>%
  pull(sample_id)

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


## ggstatplot across functional sites
set.seed(123)
counts_psi_comb <- psi_comb %>% 
  dplyr::count(Preference, Uniprot_wrapped)
plot_dsp <-  ggplot(psi_comb, aes(Uniprot_wrapped, dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 5, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Disulfide\nBond", "Localization\nSignal"),
                                        c("Disulfide\nBond", "Modifications"),
                                        c("Disulfide\nBond", "Other"),
                                        c("Localization\nSignal", "Modifications"),
                                        c("Localization\nSignal", "Other"),
                                        c("Modifications", "Other"))) + 
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A"))  + 
  theme_Publication() + 
  
  labs(y="Percent Spliced In (PSI)", x= "Uniprot-defined Functional Site") + 
  geom_text(data = counts_psi_comb, aes(label = paste("n =",n), x = Uniprot_wrapped, y = 0), vjust = 3, size = 4, hjust=.5) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Angles x-axis text
  ylim(c(-20,170))

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 8, height = 5)
print (plot_dsp)
dev.off()

# get and filter for kinase genes
known_kinase_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) %>%
  dplyr::rename(gene=Gene_Symbol) %>% 
  dplyr::filter(type=='Kinase')

psi_unip_kinase <- dplyr::inner_join(psi_comb, known_kinase_df, by='gene') 
counts_psi_unip_kinase <- psi_unip_kinase %>% dplyr::count(Preference )

## make sina plot
set.seed(45)
kinase_dpsi_plot <- ggplot(psi_unip_kinase,aes(Preference,dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  geom_label_repel(box.padding = 0.5, min.segment.length = 0.5,max.overlaps =Inf, aes(label = gene), data=psi_unip_kinase %>% 
                     subset(gene %in% c("CLK1")), size=2) +
  scale_color_manual(name = "Preference", values = c(Skipping = "black", Inclusion = "black")) + 
  theme_Publication() +
  labs(y="Percent Spliced In (PSI)") + 
  geom_text(data = counts_psi_unip_kinase, aes(label = paste("n =",n), x = Preference, y = 0), vjust = 3, size = 4, hjust=.5) +
  theme(legend.position="none")

pdf(file_dpsi_kinase_plot, 
    width = 4, height = 4)
kinase_dpsi_plot
dev.off()

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

##kegg gene set
ora_incl_results <- enricher(
  gene = unique(kinase_incl_pref$gene), # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_skip_results <- enricher(
  gene = unique(kinase_skip_pref$gene), # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

options(enrichplot.colours = c("darkorange","blue"))
enrich_skip_plot <- enrichplot::dotplot(ora_skip_results,
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


options(enrichplot.colours = c("darkorange","blue"))
enrich_incl_plot <- enrichplot::dotplot(ora_incl_results,
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


## since kegg pathways have and show more enrichment, we will save those
plot_pathways <- plot_grid(enrich_skip_plot,enrich_incl_plot,align="hv",
                                labels = c('Exon Skipping', 'Exon Inclusion'))
## save ORA dotplot
ggplot2::ggsave(filename = ora_dotplot_path,
                plot = plot_pathways, 
                width=17,
                height=7,
                device="pdf")

kinase_skip_pref <- kinase_skip_pref %>%  
  dplyr::mutate('Exon Coordinates' = str_match(SpliceID, "(\\w+)\\:(\\d+\\-\\d+)\\_")[, 3]) %>%
  unique() 

kinase_incl_pref <- kinase_incl_pref %>% 
  dplyr::mutate('Exon Coordinates' = str_match(SpliceID, "(\\w+)\\:(\\d+\\-\\d+)\\_")[, 3]) %>%
  unique()

kinase_pref <- rbind(kinase_skip_pref, kinase_incl_pref) %>% 
  dplyr::select(SpliceID,dPSI,Uniprot, gene, Preference,`Exon Coordinates`)


total_diff_events <- vroom(file.path(results_dir,"splice_events.diff.SE.txt")) %>%
  dplyr::rename(SpliceID="Splice ID") %>%
  dplyr::count(SpliceID)  %>% 
  #inner_join(kinase_pref,by="SpliceID") %>%
  dplyr::rename('Frequency'=n) %>%
  filter(Frequency > 1) # Add recurrence filter
  #dplyr::mutate(Baseline_freq=69-Frequency) # 69 samples in cluster - I don't think its used 

# read in exp file
#clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")
expr_file <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")


## load histologies info for HGG subtypes
#hgg_bs_id <- vroom(clin_file) %>%
#  # Select only "RNA-Seq" samples
#  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
#  pull(Kids_First_Biospecimen_ID)


#histologies_df  <-  read_tsv(clin_file) %>%
#  filter(cohort == "PBTA",
#         experimental_strategy == "RNA-Seq",
#         RNA_library=='stranded',
#         plot_group %in% c("DIPG or DMG", "Other high-grade glioma") ) %>%
#  pull(Kids_First_Biospecimen_ID)


exp <- readRDS(expr_file) %>%
  dplyr::select(any_of(cluster6_bs_id)) %>%
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