################################################################################
# 01-oncoprint.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 02-oncoprint-SFs.R 
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(maftools)
  library(vroom)
  library('data.table')
  library(ComplexHeatmap)
  library(circlize)
  
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
clin_file <- file.path(root_dir, "analyses", "cohort_summary", "results", "histologies-plot-group.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
tmb_file <- file.path(input_dir, "snv-mutation-tmb-coding.tsv")
cnv_file <- file.path(data_dir, "consensus_wgs_plus_freec_wxs_plus_freec_tumor_only.tsv.gz")
psi_exp_file <- file.path(root_dir, "analyses", "CLK1-splicing_correlations", "results", "clk1-nf1-psi-exp-phos-df.rds")
fus_file <- file.path(data_dir, "fusion-putative-oncogenic.tsv")
goi_file <- file.path(root_dir, "analyses","splicing-factor_dysregulation/input","splicing_factors.txt")
hugo_file <- file.path(input_dir, "hgnc-symbol-check.csv")
sbi_file <- file.path(root_dir, "analyses", "splicing_index", "results", "splicing_index.SE.txt")

## color for barplot
source(file.path(input_dir, "mutation-colors.R"))

## output file
plot_out <- file.path(plots_dir,"oncoprint-splicesosome-SFs.pdf")

# read in files
histologies_df <- read_tsv(clin_file, guess_max = 100000) 

# read in sbi file
sbi_df <- read_tsv(sbi_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(match_id, SI)

## get splicing factor list + CLKs and SRPKs
hugo_genes <- read_csv(hugo_file, skip = 1) %>%
  pull(`Approved symbol`)
goi <- readLines(goi_file) %>%
  c(hugo_genes) %>%
  unique()

indep_rna_df <- vroom(indep_rna_file) %>% 
  dplyr::filter(cohort == 'PBTA') %>%
  # get match id for DNA samples
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id", "plot_group", "RNA_library")]) %>%
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma"),
         RNA_library %in% c("stranded", "poly-A stranded"))

matched_dna_samples <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         is.na(RNA_library),
         !is.na(pathology_diagnosis),
         match_id %in% indep_rna_df$match_id)

tmb_df <- read_tsv(tmb_file) %>%
  filter(Tumor_Sample_Barcode %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  left_join(histologies_df[c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  mutate(tmb_status = case_when(tmb <10 ~ "Normal",
                                tmb >=10 & tmb < 100 ~ "Hypermutant",
                                tmb >= 100 ~ "Ultra-hypermutant")) %>%
  # remove WXS for a sample we have WGS for
  filter(Kids_First_Biospecimen_ID != "BS_QF7M4SHH")

splice_df <-  readRDS(psi_exp_file) %>%
  rownames_to_column("match_id") %>%
  filter(match_id %in% indep_rna_df$match_id) %>%
  left_join(sbi_df) %>%
  # z-score
  dplyr::mutate(`CLK1-201 (Exon4) PSI` = as.numeric(scale(`CLK1-201 (Exon4) PSI`)),
                `NF1-215 PSI` = as.numeric(scale(`NF1-215 PSI`)),
                `CLK1-201` = as.numeric(scale(`CLK1-201`)),
                `Total CLK1` = as.numeric(scale(`Total CLK1`)),
                `Total NF1` = as.numeric(scale(`Total NF1`)),
                SI = as.numeric(scale(SI)))

# read in cnv file and reformat to add to maf
cnv_df <- read_tsv(cnv_file) %>%
  # select only goi, DNA samples of interest
  filter(gene_symbol %in% goi,
         biospecimen_id %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  mutate(Variant_Classification = case_when(status %in% c("amplification", "Amplification") ~ "Amp",
                                            copy_number > 2*ploidy ~ "Amp",
                                            status %in% c("deep deletion", "Deep deletion") ~ "Del",
                                            copy_number == 0 ~ "Del",
                                            TRUE ~ NA_character_)) %>%
  filter(!is.na(Variant_Classification)) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = biospecimen_id,
                Hugo_Symbol = gene_symbol) %>%
  select(Kids_First_Biospecimen_ID, Hugo_Symbol, Variant_Classification)

# read in fusion file and reformat to add to maf
fus_df <- read_tsv(fus_file) %>%
  # select only goi, DNA samples of interest
  filter(Sample %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
  mutate(Variant_Classification = "Fusion") %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  select(Kids_First_Biospecimen_ID, Gene1A, Gene1B, Gene2A, Gene2B, Variant_Classification) %>%
  unique() %>%
  pivot_longer(cols = starts_with("Gene"), 
               names_to = "colname", 
               values_to = "Hugo_Symbol") %>%
  filter(!is.na(Hugo_Symbol),
         Hugo_Symbol %in% goi) %>%
  select(Kids_First_Biospecimen_ID, Hugo_Symbol, Variant_Classification) %>%
  unique()

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
maf_filtered <- maf %>%
  dplyr::filter(Hugo_Symbol %in% goi,
                Tumor_Sample_Barcode  %in% matched_dna_samples$Kids_First_Biospecimen_ID,
                Variant_Classification %in% names(colors)) %>%
  dplyr::mutate(keep = case_when(Variant_Classification == "Missense_Mutation" & (grepl("dam", PolyPhen) | grepl("deleterious\\(", SIFT)) ~ "yes",
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

# get genes in order of most to least mutations
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
  write_tsv(file.path(results_dir, "onco_matrix_SFs.tsv"))

rownames(gene_matrix) <- gene_matrix$Hugo_Symbol 

gene_matrix <- gene_matrix %>%
  select(-c(Sort_Order, Hugo_Symbol))

# how many samples/what percent have alterations?
print(paste(length(unique(collapse_snv_dat$match_id)), "tumors altered in", 
            length(unique(matched_dna_samples$match_id)), "patients"))
print(paste(round(length(unique(collapse_snv_dat$match_id))/length(unique(matched_dna_samples$match_id))*100, 2),"%"))

collapse_snv_dat_mut <- collapse_snv_dat %>%
  filter(!Variant_Classification %in% c("Amp", "Del"))
print(paste(length(unique(collapse_snv_dat_mut$match_id)), "tumors mutated in", 
            length(unique(matched_dna_samples$match_id)), "patients"))
print(paste(round(length(unique(collapse_snv_dat_mut$match_id))/length(unique(matched_dna_samples$match_id))*100, 2),"%"))

# mutate the hgg dataframe for plotting
histologies_df_sorted <- splice_df %>%
  left_join(histologies_df, by = "match_id", relationship = "many-to-many") %>%
  select(match_id, plot_group, reported_gender, molecular_subtype, CNS_region, `CLK1-201 (Exon4) PSI`, SI,
         `NF1-215 PSI`, `CLK1-201`, `Total CLK1`, `Total NF1`, `NF1 pS864`, `NF1 pS2796`, `Total NF1 Protein`) %>%
  unique() %>%
  group_by(match_id, plot_group, reported_gender, molecular_subtype, `CLK1-201 (Exon4) PSI`, SI,
           `NF1-215 PSI`, `CLK1-201`, `Total CLK1`, `Total NF1`, `NF1 pS864`, `NF1 pS2796`, `Total NF1 Protein`) %>%
  summarise(CNS_region = str_c(unique(na.omit(CNS_region)), collapse = ","),
            CLK1_PSI = mean(`CLK1-201 (Exon4) PSI`),
            .groups = "drop") %>%
  left_join(unique(tmb_df[,c("tmb_status", "match_id")])) %>%
  filter(match_id %in% names(gene_matrix)) %>%
  # unset rownames
  column_to_rownames("match_id") %>%
  arrange(SI) %>%
  dplyr::mutate(molecular_subtype = gsub(", TP53", "", molecular_subtype),
                molecular_subtype = case_when(grepl("To be classified", molecular_subtype) ~ "To be classified",
                                              TRUE ~ molecular_subtype),
                CNS_region = case_when(CNS_region == "" ~ NA_character_,
                                       TRUE ~ CNS_region),
                tmb_status = case_when(is.na(tmb_status) ~ "Unknown",
                                       TRUE ~ tmb_status)) 

# get SI high/low
quantiles_si <- quantile(histologies_df_sorted$SI, probs=c(.25, .75), na.rm = TRUE)
lower_si <- quantiles_si[1]
upper_si <- quantiles_si[2]

histologies_df_sorted <- histologies_df_sorted %>%
  mutate(sbi_status = case_when(SI > upper_si ~ "High",
                                SI < lower_si ~ "Low",
                                TRUE ~ "Middle")) %>%
  filter(tmb_status == "Normal")

histologies_df_sorted2 <- histologies_df_sorted %>%
  select(reported_gender, plot_group, molecular_subtype, CNS_region, SI) %>%
  dplyr::rename("Gender"=reported_gender,
                "Histology" = plot_group,
                "Molecular Subtype"=molecular_subtype,
                "CNS Region"=CNS_region, 
                "SBI" = SI) 

# write out metadata
histologies_df_sorted2 %>%
  rownames_to_column(var = "match_id") %>%
  left_join(unique(histologies_df[,c("match_id", "Kids_First_Participant_ID")])) %>%
  write_tsv(file.path(results_dir, "oncoprint_sample_metadata_SFs.tsv"))

loc_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#882255", "#6699CC")

names(loc_cols) <- c("Hemispheric", "Midline", "Mixed", "Optic pathway", "Other", "Posterior fossa", "Spine", "Suprasellar", "Ventricles")


ha = HeatmapAnnotation(name = "annotation", 
                       df = histologies_df_sorted2,
                       col=list(
                         "Gender" = c("Male" = "#56B4E9",
                                      "Female" = "pink",
                                      "Unknown" = "whitesmoke"),
                         "Histology" = c("DIPG or DMG" = "#ff40d9",
                                         "Other high-grade glioma" = "#ffccf5"),
                         "Molecular Subtype" = c("DHG, H3 G35" = "springgreen4",
                                                 "DMG, H3 K28" = "#ff40d9",
                                                 "DIPG, H3 wildtype" = "orchid",
                                                 "HGG, H3 wildtype" = "lightpink",
                                                 "HGG, IDH" = "indianred",
                                                 "IHG, NTRK-altered" = "cornflowerblue",
                                                 "IHG, ALK-altered" = "skyblue4",
                                                 "IHG, ROS1-altered" = "lightblue1",
                                                 "HGG, PXA" = "navy",
                                                 "To be classified" = "whitesmoke"),
                         "CNS Region" = loc_cols,
                         "SBI" = colorRamp2(c(-4, 0, 4), c("darkblue","white", "red"))
                       ),
                       annotation_name_side = "right", 
                       annotation_name_gp = gpar(fontsize = 9),
                       na_col = "whitesmoke")


col = colors
df = histologies_df_sorted2

gene_matrix_sorted <- gene_matrix %>%
  select(all_of(rownames(df)))

# global option to increase space between heatmap and annotations
ht_opt$ROW_ANNO_PADDING = unit(1.25, "cm")

plot_oncoprint <- oncoPrint(gene_matrix_sorted[1:20,], get_type = function(x) strsplit(x, ",")[[1]],
                            column_names_gp = gpar(fontsize = 9), show_column_names = F,
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
                              Loss = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Loss"]), col = NA))),
                            col = col,
                            top_annotation = ha,
                            alter_fun_is_vectorized = TRUE,
                            #bottom_annotation = ha1,
                            column_order =  colnames(gene_matrix_sorted))

# Save plot as PDF
pdf(plot_out, width = 9, height = 3.5)
plot_oncoprint
dev.off()



# create df for enrichment
ids <- histologies_df_sorted %>%
  rownames_to_column(var = "match_id") %>%
  select(match_id, sbi_status)

total_high <- nrow(ids)/2
total_low <- nrow(ids)/2


alteration_counts <- collapse_snv_dat %>%
  full_join(ids) %>%
  filter(sbi_status != "Middle") %>%
  ## group by junction and calculate means
  select(Hugo_Symbol, sbi_status) %>%
  group_by(Hugo_Symbol, sbi_status) %>%
  count() %>%
  ungroup() %>%
  # Spread to wide format to get separate columns for "High" and "Low"
  pivot_wider(names_from = sbi_status, values_from = n, values_fill = list(n = 0)) %>%
  rowwise() %>%
  mutate(
    Fisher_Test = list(
      fisher.test(
        matrix(
          c(High, total_high - High,  # Counts of High and the absence of High
            Low, total_low - Low),    # Counts of Low and the absence of Low
          nrow = 2
        )
      )
    ),
    P_Value = Fisher_Test$p.value
  ) %>%
  select(Hugo_Symbol, High, Low, P_Value) %>%
  ungroup() %>%
  arrange(P_Value) %>%
  write_tsv(file.path(results_dir, "sbi_high_low_mutation_counts_SFs.tsv"))



