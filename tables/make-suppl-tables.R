library(tidyverse)
library(openxlsx)
library(vroom)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses")
table_dir <- file.path(root_dir, "tables")
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(table_dir, "input")


# output directory for supplementary tables
supp_tables_dir <- file.path(table_dir, "output")
if(!dir.exists(supp_tables_dir)){
  dir.create(supp_tables_dir, recursive=TRUE)
}

# define input files
histology_file <- file.path(data_dir, "histologies.tsv")
histology_se_events <- file.path(analysis_dir, "histology-specific-splicing", "results", "unique_events.SE.tsv")
cluster_membership <- file.path(analysis_dir, "sample-psi-clustering", "results", "sample-cluster-metadata-top-5000-events-stranded.tsv")
gsva <- file.path(analysis_dir, "sample-psi-clustering", "results", "all_gsva_de_results_stranded.tsv")
CNS_match_json <- file.path(table_dir, "input", "CNS_primary_site_match.json")

kegg_file <- file.path(root_dir, "analyses", "sample-psi-clustering", "results", "hallmark_kegg_splice_geneset_mrna.rds")
sf_list_file <- file.path(root_dir, "analyses","splicing-factor_dysregulation", "input", "splicing_factors.txt")
hugo_file <- file.path(root_dir, "analyses", "oncoprint", "input", "hgnc-symbol-check.csv")
hugo_maf_file <- file.path(root_dir, "analyses", "oncoprint", "results", "hugo-all-samples.maf")
sf_maf_file <- file.path(root_dir, "analyses", "oncoprint", "results", "sf-all-samples.maf")
deseq2_sf_file <- file.path(analysis_dir, "splicing-factor_dysregulation", "results", "cluster6-diffSFs_sig_genes.txt")

func_sites_es_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.SE.total.neg.intersectunip.ggplot.txt") 
func_sites_ei_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.SE.total.pos.intersectunip.ggplot.txt") 
kinase_func_sites_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing-factor-kinases-functional_sites.tsv")
clk1_ex4_prop <- file.path(analysis_dir, "CLK1-splicing_correlations", "results", "clk1-exon4-proportion.tsv")

deseq2_morph_file <- file.path(root_dir,"analyses/CLK1-splicing-impact-morpholino","results","ctrl_vs_treated.de.tsv")
rmats_tsv_file <- file.path(data_dir,"ctrl-vs-morpholino-merged-rmats.tsv")
clk1_consens_targets_file<-file.path(root_dir,"analyses/KNS42-cell-line/results/clk1_consensus_targets.tsv")
func_sites_SE_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.SE.intersectUnip.ggplot.txt")
func_sites_A5SS_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.A5SS.intersectUnip.ggplot.txt")
func_sites_A3SS_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.A3SS.intersectUnip.ggplot.txt")
func_sites_RI_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.RI.intersectUnip.ggplot.txt")
de_ds_genes_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","de_ds_genes.txt")

cluster_cor_file <- file.path(root_dir, "analyses/splicing-factor_dysregulation/results/se-sbi-sf-expr-correlations.tsv")

func_sites_goi_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results", "differential_splice_by_goi_category.tsv")
primers_file <-  file.path(input_dir,"primers.tsv")
ds_de_crispr_events_file <-  file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results", "ds-de-crispr-events.tsv")

# define suppl output files and sheet names, when appropriate
table_s1_file <- file.path(supp_tables_dir, "TableS1-histologies.xlsx")
table_s2_file <- file.path(supp_tables_dir, "TableS2-histology-specific-splice-events.xlsx")
table_s3_file <- file.path(supp_tables_dir, "TableS3-SF-dysreg.xlsx")
table_s4_file <- file.path(supp_tables_dir, "TableS4-functional-sites.xlsx")
table_s5_file <- file.path(supp_tables_dir, "TableS5-cluster-expr-correlations.xlsx")
table_s6_file <- file.path(supp_tables_dir, "TableS6-CLK1-ex4-splicing-impact-morpholino.xlsx")

## write table for histologies
# Sheet 1: README tab
histology_df <- read_tsv(histology_file, guess_max = 100000)
histology_df_tumor <- histology_df %>%
  filter(!is.na(pathology_diagnosis))

readme <- tribble(
  ~`Histology column`,~Definition,~`Possible values`,
  "age_at_chemo_start","Patient age at chemotherapy start in days","numeric",
  "age_at_diagnosis_days","Patient age at diagnosis in days","numeric",
  "age_at_event_days","Patient age at sample collection event in days","numeric",
  "age_at_radiation_start","Patient age at radiation start in days","numeric",
  "age_last_update_days","Patient age at the last clinical event/update in days","numeric",
  "aliquot_id","External aliquot identifier","alphanumeric",
  "broad_histology","Broad WHO classification of cancer type",paste(unique(histology_df$broad_histology), collapse = "; "),
  "cancer_group","Harmonized cancer groupings for plots",paste(unique(histology_df$cancer_group), collapse = "; "),
  "cancer_predispositions","Reported cancer predisposition syndromes",paste(unique(histology_df$cancer_predispositions), collapse = "; "),
  "cell_line_composition","Cell line media",paste(unique(histology_df$cell_line_composition), collapse = "; "),
  "cell_line_passage","Cell line passage at collection","numeric",
  "clinical_status_at_event","Patient status at the time of sample collection", paste(unique(histology_df$clinical_status_at_event), collapse = "; "),
  "CNS_region","Harmonized brain region based on `primary_site`",paste(unique(histology_df$CNS_region), collapse = "; "),
  "cohort","Scientific cohort",paste(unique(histology_df$cohort), collapse = "; "),
  "cohort_participant_id","Scientific cohort participant ID","C#####-C######",
  "composition","Sample composition",paste(unique(histology_df$composition), collapse = "; "),
  "dkfz_v11_methylation_subclass","v11b6 DKFZ methylation-based CNS tumor subclass","text",
  "dkfz_v11_methylation_subclass_score","v11b6 DKFZ methylation-based CNS tumor subclass score","numeric",
  "dkfz_v12_methylation_subclass","v12b6 DKFZ methylation-based CNS tumor subclass score","text",
  "dkfz_v12_methylation_subclass_score","v12b6 DKFZ methylation-based CNS tumor subclass","numeric",
  "dkfz_v12_methylation_mgmt_status","v12b6 DKFZ MGMT promoter methylation status",paste(unique(histology_df$dkfz_v11_methylation_subclass), collapse = "; "),
  "dkfz_v12_methylation_mgmt_estimated","v12b6 DKFZ MGMT promoter methylation fraction","numeric",
  "EFS_days","Event-free survival in days","numeric",
  "EFS_event_type", "Event considered when calculating EFS", paste(unique(histology_df$EFS_event_type), collapse = "; "),
  "ethnicity","Patient reported ethnicity",paste(unique(histology_df$ethnicity), collapse = "; "),
  "experimental_strategy","Sequencing strategy",paste(unique(histology_df$experimental_strategy), collapse = "; "),
  # leaving this non-programmatic because of the duplicates that would come up (eg two selections in one patient, needing data cleanup)
  "extent_of_tumor_resection","Amount of tumor resected at time of surgical event","Biopsy only;Partial resection;Gross/Near total resection;Not Reported;Unavailable",
  "germline_sex_estimate","Predicted sex of patient based on germline X and Y ratio calculation (described in methods)",paste(unique(histology_df$germline_sex_estimate), collapse = "; "),
  "gtex_group","Tissue Type",paste(unique(histology_df$gtex_group), collapse = "; "),
  "gtex_subgroup","Tissue Subtype",paste(unique(histology_df$gtex_subgroup), collapse = "; "),
  "harmonized_diagnosis","`integrated_diagnosis` if exists or updated and harmonized diagnosis using pathology_free_text_diagnosis information","text",
  "integrated_diagnosis","WHO 2021 diagnosis integrated from pathology diagnosis and molecular subtyping","text",
  "Kids_First_Biospecimen_ID","Biospecimen identifier, Kids First or other cohort","BS_########",
  "Kids_First_Participant_ID","Patient identifier, Kids First or other cohort","PT_########",
  "match_id", "ID used to match experimental strategies within an event per sample composition", "Concatenation of sample_id, tumor descriptor, composition, and cell line composition and passage if applicable",
  "molecular_subtype","Molecular subtype defined by WHO 2021 guidelines","text",
  "molecular_subtype_methyl","DKFZ v12b6 or NIH v2 methylation class aligned to WHO 2021 subtypes","text",
  "normal_fraction","Theta2 normal DNA fraction estimate","numeric",
  "Notes","Free text field describing changes from `pathology_diagnosis` to `integrated_diagnosis` or manner in which molecular_subtype was determined","text",
  "OS_days","Overall survival in days","numeric",
  "OS_status","Overall survival status",paste(unique(histology_df$OS_status), collapse = "; "),
  "pathology_diagnosis","Reported and/or harmonized patient diagnosis from pathology reports","text",
  "pathology_free_text_diagnosis","Free text patient diagnosis from pathology reports","text",
  "primary_site","Bodily site(s) from which specimen was derived","text",
  "race","Patient reported race",paste(unique(histology_df$race), collapse = "; "),
  "reported_gender","Patient reported gender",paste(unique(histology_df$reported_gender), collapse = "; "),
  "RNA_library","Type of RNA-Sequencing library preparation",paste(unique(histology_df$RNA_library), collapse = "; "),
  "sample_id","Event id","alphanumeric",
  "sample_type","Broad sample type",paste(unique(histology_df$sample_type), collapse = "; "),
  "seq_center","Sequencing center",paste(unique(histology_df$seq_center), collapse = "; "),
  "short_histology","Abbreviated `cancer_group` or `broad_histology` for plotting purposes",paste(unique(histology_df$short_histology), collapse = "; "),
  "sub_cohort", "sub-cohort", paste(unique(histology_df$sub_cohort), collapse = "; "),
  "tumor_descriptor","Phase of therapy from which tumor was derived",paste(unique(histology_df$tumor_descriptor), collapse = "; "),
  "tumor_fraction","Theta2 tumor DNA fraction estimate","numeric",
  "tumor_fraction_LUMP","LUMP tumor DNA fraction estimate from methylation","numeric",
  "tumor_fraction_RFpurify_ABSOLUTE","RFpurify ABSOLUTE tumor DNA fraction estimate from methylation","numeric",
  "tumor_fraction_RFpurify_ESTIMATE","RFpurify ESTIMATE tumor DNA fraction estimate from methylation","numeric",
  "tumor_ploidy","Control-FREEC ploidy estimate","numeric"
)
# Sheet 2: Histologies file (histology_df)

# Sheet 3: CNS region definition based on definitions from [Cassie Kline](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/tables/input/CNS_primary_site_match.json) with additional manual review of [HGG primary_site](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1025).
# This is integrated upstream of molecular subtyping

# CNS_region ~ primary_site matches
cns_regions_df <- purrr::imap_dfr(jsonlite::fromJSON(CNS_match_json),
                                  function(x, name) { tibble::tibble(CNS_region = name, primary_site = paste0(x, collapse = ';')) })


# Combine and output
list_s1_table <- list(A_README = readme,
                      B_histologies_file = histology_df,
                      C_CNS_region_definitions = cns_regions_df
)
write.xlsx(list_s1_table,
           table_s1_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 2 histology specific splicing events
## sheet 1, se events splicing
se_events_df <- vroom(histology_se_events)

## sheet 2, cluster membership
cluster_membership_df <- read_tsv(cluster_membership) %>%
  dplyr::select(sample_id, plot_group,
                RNA_library, molecular_subtype,
                cluster) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id,
                Histology = plot_group) %>%
  # clean up NAs for excel output
  dplyr::mutate(molecular_subtype = case_when(
    is.na(molecular_subtype) ~ "N/A",
    TRUE ~ molecular_subtype
  ))

## sheet 3, gsva scores
gsva_df <- read_tsv(gsva) %>%
  dplyr::mutate(Filename = as.numeric(stringr::str_extract(Filename, "\\d+"))) %>%
  dplyr::rename(Cluster = Filename) %>%
  arrange(Cluster)
          
# Combine and output
list_s2_table <- list(A_se_events = se_events_df,
                      B_clust_memb = cluster_membership_df,
                      C_gsva_scores = gsva_df)

write.xlsx(list_s2_table,
           table_s2_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 3 Splicing factor dysregulation
#tab 1, kegg list
kegg_rds <- readRDS(kegg_file)
kegg_list <- kegg_rds$KEGG_SPLICEOSOME
#tab 2, hugo list
hugo_list <- read_csv(hugo_file, skip = 1) %>%
  pull(`Approved symbol`)
#tab 3, sf list
sf_list <- readLines(sf_list_file) 

# tab 4
hugo_maf <- read_tsv(hugo_maf_file)

# tab 5
sf_maf <- read_tsv(sf_maf_file)

## tab 6, exon inclusion splicing
deseq_df <- vroom(deseq2_sf_file) %>%
  dplyr::select(-Significant)

# Combine and output
list_s3_table <- list(A_kegg_spliceosome = kegg_list, 
                      B_hugo_spliceosome = hugo_list,
                      G_sf_genes = sf_list,
                      E_hugo_mutations = hugo_maf,
                      F_sf_mutations = sf_maf,
                      G_high_v_low_sbi_deseq2 = deseq_df)

write.xlsx(list_s3_table,
           table_s3_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 4 Differential HGG splicing events impacting functional sites
## sheet 1, exon inclusion events
ds_events_es_df <- vroom(func_sites_es_file)
ds_events_ei_df <- vroom(func_sites_ei_file)
ds_events_A3SS_df <- vroom(func_sites_A3SS_morpho_tsv_file)
ds_events_RI_df <- vroom(func_sites_RI_morpho_tsv_file)


## sheet 2, exon skipping events
es_events_df <- vroom(func_sites_es_file)

##sheet 3 kinases
kinase_events_df <- vroom(kinase_func_sites_file)

## sheet 4 exon 4 proportions
clk1_ex4_prop_df <- vroom(clk1_ex4_prop)

# Combine and output
list_s4_table <- list(A_ds_skipping = ds_events_es_df,
                      B_ds_inclusion = ds_events_ei_df,
                      C_prioritized_sf_kinases = kinase_events_df,
                      D_clk1_ex4_prop = clk1_ex4_prop_df)

write.xlsx(list_s4_table,
           table_s4_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 5 Cluster expression correlations
cluster_cor_df <- read_tsv(cluster_cor_file)

list_s5_table <- list(A_cluster_exp_cor = cluster_cor_df)

write.xlsx(list_s5_table,
           table_s5_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 6 morpholino vs ctrl DESeq2 and rMATs results
deseq2_morpholino_df <- vroom(deseq2_morph_file) %>%
  filter(padj < 0.05) %>%
  dplyr::select(Gene_Symbol,	
                baseMean,	
                log2FoldChange,	
                lfcSE,	
                stat,	
                pvalue,	
                padj)

rmats_df <-  vroom(rmats_tsv_file)
ds_events_SE_df <- vroom(func_sites_SE_morpho_tsv_file)
ds_events_A5SS_df <- vroom(func_sites_A5SS_morpho_tsv_file)
ds_events_A3SS_df <- vroom(func_sites_A3SS_morpho_tsv_file)
ds_events_RI_df <- vroom(func_sites_RI_morpho_tsv_file)
de_ds_genes_file_df <- read_lines(de_ds_genes_file)

primers_df <- vroom(primers_file, delim = "\t")

ds_de_crispr_df <-  vroom(ds_de_crispr_events_file) %>%
  mutate(across(everything(), ~ replace_na(as.character(.), "-")))

consensus_targets_df <- vroom(clk1_consens_targets_file)

list_s6_table <- list(A_deseq2_morp = deseq2_morpholino_df,
                      B_morph_rmats = rmats_df,
                      C_ds_SE = ds_events_SE_df,
                      D_ds_A5SS = ds_events_A5SS_df,
                      E_ds_A3SS = ds_events_A3SS_df,
                      F_ds_RI = ds_events_RI_df,
                      G_de_ds_genes = de_ds_genes_file_df,
                      H_primers = primers_df,
                      I_intersect_de_ds_crispr = ds_de_crispr_df,
                      J_clk1_targets = consensus_targets_df)

write.xlsx(list_s6_table,
           table_s6_file,
           overwrite=TRUE,
           keepNA=TRUE)
