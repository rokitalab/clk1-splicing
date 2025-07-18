################################################################################
# 01-generate-cohort-summary-circos-plot.R
# Generate circos summary plot of cohort used
# Author: Ammar Naqvi, Shehbeel Arif, Jo Lynne Rokita
# Usage: Rscript 01-generate-cohort-summary-circos-plot.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "cohort_summary")
input_dir   <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(plots_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## output files for final plots
file_circos_plot <- file.path(analysis_dir, "plots", "cohort_circos.pdf")


# Load datasets and pre-process
hist_df <- read_tsv(file.path(data_dir,"histologies.tsv"), guess_max = 100000) %>%
  # mutate ages at dx (7316-851 is 17)
  mutate(age_at_diagnosis_days = case_when(Kids_First_Participant_ID == "PT_AEDWCP8Z" ~
                                             as.integer(365.25*17),
                                           Kids_First_Participant_ID == "PT_6PYNBA9C" ~ as.integer(365.25*5/12),
                                           TRUE ~ age_at_diagnosis_days)
  ) %>%
  # filter
  filter(cohort == "PBTA",
         !broad_histology %in% c("Eye tumor", "Hematologic malignancy", "Mixed tumor"),
         !is.na(pathology_diagnosis),
         !composition %in% c("Derived Cell Line", "Patient Derived Xenograft")) %>%
  # collapse reported gender to 3 groups
  mutate(reported_gender = case_when(reported_gender == "Not Reported" ~ "Unknown",
                                     TRUE ~ reported_gender),
         germline_sex_estimate = case_when(is.na(germline_sex_estimate) ~ "Unknown",
                                        TRUE ~ germline_sex_estimate),
         # update 7316-3066
         broad_histology = case_when(sample_id == "7316-3066" ~ "Tumor of cranial and paraspinal nerves", 
                                     broad_histology == "Other" ~ "Other tumor",
                                     cancer_group == "Glial-neuronal tumor" ~ "Neuronal and mixed neuronal-glial tumor",
                                     cancer_group %in% c("Cavernoma", "Malignant peripheral nerve sheath tumor") ~ "Benign tumor",
                                     TRUE ~ broad_histology),
         cancer_group = case_when(sample_id == "7316-3066" ~ "Neurofibroma/Plexiform",
                                  grepl("xanthogranuloma", pathology_free_text_diagnosis) & broad_histology == "Histiocytic tumor" & pathology_diagnosis == "Other" ~ "Juvenile xanthogranuloma",
                                  broad_histology == "Choroid plexus tumor" ~ "Choroid plexus tumor",
                                  cancer_group == "Glial-neuronal tumor" ~ "Glial-neuronal tumor NOS",
                                  cancer_group == "Low-grade glioma" ~ "Low-grade glioma/astrocytoma",
                                  cancer_group %in% c("High-grade glioma", "Astrocytoma", "Astroblastoma", "Glioblastoma", "Diffuse hemispheric glioma",
                                                      "Infant-type hemispheric glioma") ~ "High-grade glioma/astrocytoma",
                                  broad_histology == "Ependymal tumor" ~ "Ependymoma",
                                  broad_histology == "Other tumor" ~ "Other tumor",
                                  broad_histology == "Diffuse astrocytic and oligodendroglial tumor" & (is.na(cancer_group) | cancer_group == "Oligodendroglioma") ~ "High-grade glioma/astrocytoma",
                                  cancer_group %in% c("Non-germinomatous germ cell tumor", "Diffuse leptomeningeal glioneuronal tumor", "Malignant peripheral nerve sheath tumor") ~ NA_character_,
                                  broad_histology == "Meningioma" ~ "Meningioma",
                                  cancer_group == "Perineuroma" ~ "Neurofibroma/Plexiform",
                                  is.na(cancer_group) & broad_histology == "Tumor of cranial and paraspinal nerves" ~ "Neurofibroma/Plexiform",
                                  TRUE ~ cancer_group))
hist_update <- hist_df %>%
  # update PNOC and mioncoseq DIPGs
  mutate(cancer_group = case_when(pathology_diagnosis == "Brainstem glioma- Diffuse intrinsic pontine glioma" & 
                                    sub_cohort %in% c("PNOC", "Mioncoseq") &
                                    !grepl("IDH", molecular_subtype) &
                                    cancer_group == "High-grade glioma/astrocytoma" ~ "Diffuse intrinsic pontine glioma", 
         TRUE ~ cancer_group)) %>%
  mutate(molecular_subtype = case_when(cancer_group == "Diffuse intrinsic pontine glioma" ~ 
                                    gsub("HGG", "DIPG", molecular_subtype),
                                  TRUE ~ molecular_subtype))

# find and store all NAs in "age_at_diagnosis_days" that are not PNOC
hist_NAs_df <- hist_update %>%
  filter(is.na(age_at_diagnosis_days) & sub_cohort != "PNOC")

# Apply the filter only when sub_cohort is not "PNOC" and remove under 40
under40_df <- hist_update %>%
  dplyr::filter(sub_cohort == "PNOC" | 
           (sub_cohort != "PNOC" & age_at_diagnosis_days < (365.25*40)) |
              (sub_cohort != "PNOC" & is.na(age_at_diagnosis_days) & age_at_event_days < (365.25*40))
  )

# add cancer/plot group mapping file 
map_file <- read_tsv(file.path(input_dir, "plot-mapping.tsv"))

# add plot mapping file and old plot groups, export this.
combined_plot_map <- under40_df %>%
  full_join(map_file, by = c("broad_histology", "cancer_group")) %>%
  select(names(map_file)) %>%
  unique() %>%
  arrange(broad_histology) %>%
  write_tsv(file.path(results_dir, "plot_mapping.tsv"))

# add plot mapping to histology df
combined_hist_map <- under40_df %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) %>%
  write_tsv(file.path(results_dir, "histologies-plot-group.tsv"))

combined_hist_map <- combined_hist_map %>%
  filter(experimental_strategy == "RNA-Seq",
         RNA_library == "stranded",
         !is.na(pathology_diagnosis))

## filter using independent specimens file 
independent_specimens_df <- read_tsv(file.path(data_dir,"independent-specimens.rnaseqpanel.primary.tsv")) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% combined_hist_map$Kids_First_Biospecimen_ID)

intersect(hist_NAs_df$Kids_First_Biospecimen_ID, independent_specimens_df$Kids_First_Biospecimen_ID)

# Merge both meta datasets
hist_indep_df <- combined_hist_map %>%
  right_join(independent_specimens_df, by="Kids_First_Biospecimen_ID") %>% 
  unique()

uniq_plot_cols <- combined_plot_map %>%
  select(plot_group, plot_group_hex) %>%
  unique()

hist_indep_df <- hist_indep_df %>%
  # add labeling that shows wrapped groups with `(n=X)`
  dplyr::count(plot_group) %>%
  dplyr::distinct() %>%
  dplyr::left_join(uniq_plot_cols) %>%
  # Create wrapped with (n=X) factor column for cancer groups
  dplyr::mutate(plot_group_n = glue::glue("{plot_group} (N={n})")) %>%
  dplyr::inner_join(hist_indep_df) %>%
  unique() %>%
  # Columns of interest: Kids_First_Biospecimen_ID, plot group, CNS_region, OS status, reported_gender
  dplyr::select(Kids_First_Biospecimen_ID, plot_group_n, plot_group, plot_group_hex, CNS_region, reported_gender) %>%
  remove_rownames() %>% 
  column_to_rownames("Kids_First_Biospecimen_ID") %>%
  arrange(plot_group_n, plot_group_hex, CNS_region, reported_gender) 

uniq_plot_cols <- uniq_plot_cols %>%
  left_join(hist_indep_df[,c("plot_group_n", "plot_group_hex")]) %>%
  unique() %>%
  arrange(plot_group)

# color palette for histology
cols <- uniq_plot_cols$plot_group_hex
names(cols) <- uniq_plot_cols$plot_group_n
#sorted_cols <- cols[order(names(cols))]

# remove hex from df
hist_indep_df <- hist_indep_df %>%
  select(-c(plot_group_hex, plot_group))
  
# Create histology factors for splitting the Circos plot
split <- factor(hist_indep_df$plot_group)

# colors for others
gender_cols <- c( "pink", "#56B4E9", "lightgrey")
names(gender_cols) <- c("Female", "Male", "Unknown")

loc_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#882255", "#6699CC")
names(loc_cols) <- c(sort(unique(hist_indep_df$CNS_region)))

cols_all <- c(cols, loc_cols, gender_cols)

pdf(file_circos_plot, width = 7, height = 7)
# Reset Circos plot
circos.clear()

# Make base Circos plot
circos.par(points.overflow.warning = FALSE)
circos.heatmap(hist_indep_df,
               split=split, # Specify how to split plot into sectors
               col = unlist(cols_all))#, # Add in color list
      
# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = "black")
}

## adding legends
lgd_plot_group = Legend(title = "Histology", 
                        at =  names(cols), 
                        legend_gp = gpar(fill = cols))

lgd_tum_loc = Legend(title = "Tumor location", 
                       at =  names(loc_cols),
                       legend_gp = gpar(fill = loc_cols))

lgd_gender = Legend(title = "Gender", 
                    at =  names(gender_cols), 
                    legend_gp = gpar(fill = gender_cols))

## Placement of legend based on device size
h = dev.size()[2]
circle_size = unit(1, "snpc")

# Create legends
lgd_list1 = packLegend(lgd_plot_group, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_tum_loc, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")

# Add legends on plot
draw(lgd_list1, x = unit(80, "mm"), y = unit(90, "mm"))
draw(lgd_list2, x = unit(129, "mm"), y = unit(95, "mm"))
draw(lgd_list3, x = unit(125, "mm"), y = unit(63, "mm"))
circos.clear()
dev.off()



