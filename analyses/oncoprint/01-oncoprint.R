################################################################################
# 01-oncoprint.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 01-oncoprint.R 
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
plots_dir   <- file.path(analysis_dir, "plots")
results_dir   <- file.path(analysis_dir, "results")

## check and create plots/results dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## input files
cons_maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
tumor_only_maf_file <- file.path(data_dir,"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz")
clin_file <- file.path(data_dir, "histologies-plot-group.tsv")
goi_file <- file.path(input_dir,"oncoprint-goi-lists-OpenPedCan-gencode-v39.csv")
tmb_file <- file.path(input_dir, "snv-mutation-tmb-coding.tsv")
cnv_file <- file.path(data_dir, "consensus_wgs_plus_freec_wxs_plus_freec_tumor_only.tsv.gz")
fus_file <- file.path(data_dir, "fusion-putative-oncogenic.tsv")
cluster_file <- file.path(root_dir, "analyses", 
                          "sample-psi-clustering", "results", 
                          "sample-cluster-metadata-top-5000-events-stranded.tsv")
sf_file <- file.path(root_dir, "analyses","splicing-factor_dysregulation/input","splicing_factors.txt")
hugo_file <- file.path(input_dir, "hgnc-symbol-check.csv")

# read in files
histologies_df <- read_tsv(clin_file, guess_max = 100000)

## color for barplot
source(file.path(input_dir, "mutation-colors.R"))

# read in files
histologies_df <- read_tsv(clin_file, guess_max = 100000) %>%
  mutate(cancer_predisposition = case_when(cancer_predispositions == "Neurofibromatosis, Type 1 (NF-1)" ~ "NF-1",
                                           cancer_predispositions == "Li-Fraumeni syndrome (TP53)" ~ "LFS",
                                           cancer_predispositions == "Other inherited conditions NOS" ~ "Other",
                                           Kids_First_Participant_ID == "PT_3CHB9PK5" ~ "CMMRD",
                                           Kids_First_Participant_ID == "PT_JNEV57VK" ~ "LS",
                                           Kids_First_Participant_ID == "PT_ZH3SBJPZ" ~ NA_character_,
                                           TRUE ~ NA_character_))

## get splicing factor list + CLKs and SRPKs
hugo_genes <- read_csv(hugo_file, skip = 1) %>%
  pull(`Approved symbol`)
sf_genes <- readLines(sf_file) %>%
  unique()

goi <- read_csv(goi_file)%>%
  select(LGAT, `Embryonal tumor`, HGAT, Other) %>%
  pivot_longer(everything(), values_to = "gene") %>%
  pull(gene) %>%
  unique()

# cat genes
all_goi <- unique(c(hugo_genes, sf_genes, goi))

cluster_df <- read_tsv(cluster_file) %>%
  rename(Kids_First_Biospecimen_ID = sample_id)

cohort_hist <- histologies_df %>% 
  filter(experimental_strategy == "RNA-Seq") %>%
  left_join(cluster_df %>% dplyr::select(Kids_First_Biospecimen_ID, cluster))

matched_dna_samples <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         match_id %in% cohort_hist$match_id)

# n DNA samples - 657
print(length(unique(matched_dna_samples$Kids_First_Participant_ID)))

tmb_df <- read_tsv(tmb_file) %>%
  filter(Tumor_Sample_Barcode %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  left_join(histologies_df[c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  mutate(tmb_status = case_when(tmb <10 ~ "Normal",
                                tmb >=10 & tmb < 100 ~ "Hypermutant",
                                tmb >= 100 ~ "Ultra-hypermutant"))

# read in cnv file and reformat to add to maf, collapse to autosomes only bc of false dels in X:
cnv_df <- read_tsv(cnv_file) %>%
  # select all goi, DNA samples of interest
  filter(!grepl("Xp|Xq|Yp|Yq", cytoband),
         gene_symbol %in% all_goi,
         biospecimen_id %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  mutate(Variant_Classification = case_when(status %in% c("amplification", "Amplification") ~ "Amp",
                                            copy_number > 2*ploidy ~ "Amp",
                                            status %in% c("deep deletion", "Deep deletion") ~ "Del",
                                            copy_number == 0 ~ "Del",
                                            TRUE ~ NA_character_)) %>%
  filter(!is.na(Variant_Classification)) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = biospecimen_id,
                Hugo_Symbol = gene_symbol) 

# read in fusion file and reformat to add to maf
fus_df <- read_tsv(fus_file) %>%
  # select all goi, RNA samples of interest
  filter(Sample %in% cohort_hist$Kids_First_Biospecimen_ID) %>%
  mutate(Variant_Classification = "Fusion") %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  select(Kids_First_Biospecimen_ID, Gene1A, Gene1B, Gene2A, Gene2B, Variant_Classification) %>%
  unique() %>%
  pivot_longer(cols = starts_with("Gene"), 
               names_to = "colname", 
               values_to = "Hugo_Symbol") %>%
  filter(!is.na(Hugo_Symbol),
         Hugo_Symbol %in% all_goi) %>%
  select(Kids_First_Biospecimen_ID, Hugo_Symbol, Variant_Classification) %>%
  unique()


# maf cols to select
maf_cols <- c("Hugo_Symbol", 
              "Chromosome", 
              "Start_Position", 
              "End_Position",
              "HGVSg",
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

## filter maf for samples with RNA splicing
maf_filtered <- maf %>%
  dplyr::filter(Hugo_Symbol %in% all_goi,
                Tumor_Sample_Barcode %in% matched_dna_samples$Kids_First_Biospecimen_ID,
                Variant_Classification %in% names(colors)) %>%
  dplyr::mutate(keep = case_when(Variant_Classification == "Missense_Mutation" & (grepl("dam", PolyPhen) & grepl("deleterious\\(", SIFT)) ~ "yes",
                                 Variant_Classification == "Missense_Mutation" & PolyPhen == "" & SIFT == "" ~ "yes",
                                 Variant_Classification != "Missense_Mutation" ~ "yes",
                                 TRUE ~ "no")) %>%
  dplyr::filter(keep == "yes")

collapse_snv_dat <- maf_filtered %>%
  select(Tumor_Sample_Barcode,Hugo_Symbol,Variant_Classification) %>%
  dplyr::group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  bind_rows(cnv_df, fus_df) %>%
  dplyr::summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
                   Variant_Classification=str_c(unique(Variant_Classification),collapse = ",")) %>%
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(-Kids_First_Biospecimen_ID)

# get genes in order of most to least mutations and are in enrichment results
gene_row_order <- collapse_snv_dat %>%
  count(Hugo_Symbol) %>%
  arrange(-n)

# complex heatmap
gene_matrix <- reshape2::acast(collapse_snv_dat,
                               Hugo_Symbol ~ match_id,
                               value.var = "Variant_Classification",
                               fun.aggregate = function(x) paste(unique(x), collapse = ", ")) %>%
  as.data.frame() %>%
  dplyr::mutate_if(is.character, ~replace_na(.,"")) %>%
  # add multi-hits
  mutate(across(everything(), ~if_else(str_detect(., ","), "Multi_Hit", .))) %>%
  rownames_to_column(var = "Hugo_Symbol") %>%
  mutate(Sort_Order = match(Hugo_Symbol, gene_row_order$Hugo_Symbol)) %>%
  arrange(Sort_Order)  %>%
  write_tsv(file.path(results_dir, "onco_matrix.tsv"))

rownames(gene_matrix) <- gene_matrix$Hugo_Symbol 

gene_matrix <- gene_matrix %>%
  select(-c(Sort_Order, Hugo_Symbol)) 

# how many of top 30 genes are sfs? 
sort(unique(intersect(rownames(gene_matrix[1:30,]), unique(c(sf_genes, hugo_genes)))))
#  "CHD2"   "CHD3"   "DDX1"   "FIP1L1" "MACF1"  "PRKDC"  "QKI"

# mutate the dataframe for plotting
histologies_df_sorted <- cohort_hist %>%
  select(match_id, plot_group, cancer_predisposition, reported_gender, molecular_subtype, CNS_region, cluster) %>%
  group_by(match_id, plot_group, cancer_predisposition, reported_gender, molecular_subtype, cluster) %>%
  summarise(CNS_region = str_c(unique(na.omit(CNS_region)), collapse = ","),
            .groups = "drop") %>%
  left_join(unique(tmb_df[,c("tmb_status", "match_id")])) %>%
  filter(match_id %in% names(gene_matrix)) %>%
  # unset rownames
  column_to_rownames("match_id") %>%
  dplyr::mutate(molecular_subtype = gsub(", TP53", "", molecular_subtype),
                molecular_subtype = case_when(grepl("To be classified", molecular_subtype) ~ "To be classified",
                                              TRUE ~ molecular_subtype),
                CNS_region = case_when(CNS_region == "" ~ NA_character_,
                                       TRUE ~ CNS_region),
                tmb_status = case_when(is.na(tmb_status) ~ "Unknown",
                                       TRUE ~ tmb_status)) 

histologies_df_sorted2 <- histologies_df_sorted %>%
  select(reported_gender, cancer_predisposition, plot_group, CNS_region, tmb_status, cluster) %>%
  dplyr::rename("Gender"=reported_gender,
                "Histology" = plot_group,
                "Predisposition" = cancer_predisposition,
                "CNS Region"=CNS_region, 
                "Mutation Status"=tmb_status,
                "Cluster" = cluster) %>%
  dplyr::mutate(Cluster = factor(Cluster)) %>%
  dplyr::mutate(Cluster = fct_relevel(Cluster, mixedsort(as.character(unique(Cluster))))) %>%
  arrange(Cluster)

# write out metadata
histologies_df_sorted2 %>%
  rownames_to_column(var = "match_id") %>%
  left_join(unique(cohort_hist[,c("match_id", "Kids_First_Participant_ID")])) %>%
write_tsv(file.path(results_dir, "oncoprint_sample_metadata.tsv"))

loc_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#882255", "#6699CC")

names(loc_cols) <- c("Hemispheric", "Midline", "Mixed", "Optic pathway", "Other", "Posterior fossa", "Spine", "Suprasellar", "Ventricles")



ha = HeatmapAnnotation(name = "annotation", 
                       df = histologies_df_sorted2,
                       col=list(
                         "Gender" = c("Male" = "#56B4E9",
                                      "Female" = "pink",
                                      "Unknown" = "whitesmoke"),
                         "Histology" = c("Ependymoma" =                       "#2200ff",       
                                         "Atypical teratoid rhabdoid tumor" = "#4d0d85",       
                                         "Other high-grade glioma" =          "#ffccf5",       
                                         "Low-grade glioma" =                 "#8f8fbf",       
                                         "Meningioma" =                       "#2db398",       
                                         "Diffuse midline glioma" =           "#ff40d9",       
                                         "Medulloblastoma" =                  "#a340ff",       
                                         "Rare CNS tumor" =                   '#b5b5b5',       
                                         "Mesenchymal tumor" =                "#7fbf00",       
                                         "Craniopharyngioma" =                "#b2502d",       
                                         "Mixed neuronal-glial tumor" =       "#685815",       
                                         "Non-neoplastic tumor" =             "#FFF5EB",       
                                         "Choroid plexus tumor" =             "#00441B",       
                                         "Schwannoma" =                       "#ab7200",       
                                         "Neurofibroma plexiform" =           "#e6ac39",       
                                         "Other CNS embryonal tumor" =        "#b08ccf",       
                                         "Germ cell tumor" =                  "#0074d9"),
                         "Predisposition" = c("LFS" = "red",
                                              "NF-1" = "black",
                                              "Other" = "grey"),
                         "CNS Region" = loc_cols,
                         "Mutation Status" = c("Normal" = "grey80",
                                               "Hypermutant" = "orange",
                                               "Ultra-hypermutant" = "red",
                                               "Unknown" = "whitesmoke"),
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
                                       "11" = "#6A3D9A")
                       ),
                       annotation_name_side = "right", 
                       annotation_name_gp = gpar(fontsize = 9),
                       na_col = "whitesmoke")

col = colors

gene_matrix_sorted <- gene_matrix %>%
  select(all_of(rownames(histologies_df_sorted2)))

# global option to increase space between heatmap and annotations
ht_opt$ROW_ANNO_PADDING = unit(1.25, "cm")

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "whitesmoke",col="whitesmoke")),
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Missense_Mutation"]),col = NA)),
  Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonsense_Mutation"]),col = NA)),
  Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Del"]), col = NA)),
  Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Ins"]), col = NA)),
  Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Splice_Site"]), col = NA)),
  Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Translation_Start_Site"]), col = NA)),
  Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonstop_Mutation"]),col = NA)),
  In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Del"]),col = NA)),
  In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Ins"]), col = NA)),
  Stop_Codon_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Stop_Codon_Ins"]), col = NA)),
  Start_Codon_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Start_Codon_Del"]), col = NA)),
  Fusion = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Fusion"]),col = NA)),
  Multi_Hit_Fusion  = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit_Fusion"]),col = NA)),
  Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit"]), col = NA)),
  Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Del"]), col = NA)),
  Amp = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Amp"]), col = NA)),
  Loss = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Loss"]), col = NA)))


# sort samples by histology and cluster and regenerate oncoprint
hist_cluster_order <- histologies_df_sorted2 %>%
  arrange(Histology, Cluster) %>%
  rownames_to_column("match_id") %>%
  pull(match_id)

plot_oncoprint_hist_cluster <- oncoPrint(gene_matrix_sorted[1:30,], get_type = function(x) strsplit(x, ",")[[1]],
                            column_names_gp = gpar(fontsize = 9), show_column_names = F,
                            row_names_gp=gpar(fontsize = 11, fontface="italic"),
                            alter_fun = alter_fun,
                            col = col,
                            top_annotation = ha,
                            alter_fun_is_vectorized = TRUE,
                            column_order =  hist_cluster_order)

# Save plot as PDF
pdf(file.path(plots_dir,"oncoprint-hist-cluster.pdf"),
    width = 15, height = 7)
plot_oncoprint_hist_cluster
dev.off()

for (group_n in unique(histologies_df_sorted2$Cluster)) {
  ha_cluster = HeatmapAnnotation(name = "annotation", 
                         df = histologies_df_sorted2 %>% dplyr::filter(Cluster == group_n) %>% dplyr::select(-Cluster),
                         col=list(
                           "Gender" = c("Male" = "#56B4E9",
                                        "Female" = "pink",
                                        "Unknown" = "whitesmoke"),
                           "Histology" = c("Ependymoma" =                       "#2200ff",       
                                           "Atypical teratoid rhabdoid tumor" = "#4d0d85",       
                                           "Other high-grade glioma" =          "#ffccf5",       
                                           "Low-grade glioma" =                 "#8f8fbf",       
                                           "Meningioma" =                       "#2db398",       
                                           "Diffuse midline glioma" =           "#ff40d9",       
                                           "Medulloblastoma" =                  "#a340ff",       
                                           "Rare CNS tumor" =                   '#b5b5b5',       
                                           "Mesenchymal tumor" =                "#7fbf00",       
                                           "Craniopharyngioma" =                "#b2502d",       
                                           "Mixed neuronal-glial tumor" =       "#685815",       
                                           "Non-neoplastic tumor" =             "#FFF5EB",       
                                           "Choroid plexus tumor" =             "#00441B",       
                                           "Schwannoma" =                       "#ab7200",       
                                           "Neurofibroma plexiform" =           "#e6ac39",       
                                           "Other CNS embryonal tumor" =        "#b08ccf",       
                                           "Germ cell tumor" =                  "#0074d9"),
                           "Predisposition" = c("LFS" = "red",
                                                "LS" = "#7fbf00",
                                                "CMMRD" = "#0074d9",
                                                "NF-1" = "black",
                                                "Other" = "grey"),
                           "CNS Region" = loc_cols,
                           "Mutation Status" = c("Normal" = "grey80",
                                                 "Hypermutant" = "orange",
                                                 "Ultra-hypermutant" = "red",
                                                 "Unknown" = "whitesmoke")
                         ),
                         annotation_name_side = "right", 
                         annotation_name_gp = gpar(fontsize = 9),
                         na_col = "whitesmoke")
  
  cluster_ids <- histologies_df_sorted2 %>% 
    dplyr::filter(Cluster == group_n) %>%
    rownames_to_column("match_id") %>%
    arrange(Histology) %>%
    dplyr::pull(match_id) 
  
  gene_matrix_cluster <- gene_matrix_sorted[,colnames(gene_matrix_sorted) %in% cluster_ids]
  
  plot_oncoprint_cluster <- oncoPrint(gene_matrix_cluster[1:25,], get_type = function(x) strsplit(x, ",")[[1]],
                                      row_names_gp=gpar(fontsize = 11, fontface="italic"),
                                      column_names_gp = gpar(fontsize = 9), show_column_names = F,
                                      alter_fun = alter_fun,
                                      col = col,
                                      top_annotation = ha_cluster,
                                      alter_fun_is_vectorized = TRUE,
                                      column_order = cluster_ids)
  
  pdf(NULL)
  pdf(file.path(plots_dir,paste0("oncoprint-hist-cluster", group_n, ".pdf")),
      width = 12, height = 7)
  print(plot_oncoprint_cluster)
  dev.off()
}
