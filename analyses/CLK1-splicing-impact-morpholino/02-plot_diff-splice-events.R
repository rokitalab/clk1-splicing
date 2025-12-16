################################################################################
# 02-plot_diff-splice-events.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 02-plot_diff-splice-events.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggrepel")
  library("vroom")
  library("ggpubr")
  library("annoFuseData")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
input_dir <- file.path(analysis_dir, "input")

results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


##theme for all plots
# source function for theme for plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
de_events_plot_path <- file.path(plots_dir, "clk1_diff_splice_events_by_type.pdf")
  
## get and setup input
## rmats file
rmats_merged_file  <- file.path(data_dir,"ctrl-vs-morpholino-merged-rmats.tsv")

## extract strong splicing changes 
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
                 filter(FDR < 0.05 & PValue < 0.05) 

## extract strong differential splicing cases (dPSI >= |.10|) and use CLK1 high exon 4 as a reference (eg. -dPSI means there more inclusion in CLK1 exon 4)
splicing_df_skip <- splicing_df %>% filter(IncLevelDifference  >= .10) %>% mutate(Preference="Skipping")
splicing_df_incl <- splicing_df %>% filter(IncLevelDifference <= -.10) %>% mutate(Preference="Inclusion",
                                                                                IncLevelDifference = abs(IncLevelDifference) )

psi_comb <- rbind(splicing_df_incl,splicing_df_skip) %>%
  dplyr::mutate(
    splice_id = case_when(splicing_case == "RI" ~ glue::glue("{GeneID}:{riExonStart_0base}-{riExonEnd}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}"),
                          splicing_case == "SE" ~ glue::glue("{GeneID}:{exonStart_0base}-{exonEnd}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}"),
                          splicing_case == "MXE" ~ glue::glue("{GeneID}:{`1stExonStart_0base`}-{`1stExonEnd`}_{`2ndExonStart_0base`}-{`2ndExonEnd`}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}"),
                          TRUE ~ glue::glue("{GeneID}:{longExonStart_0base}-{longExonEnd}_{shortES}-{shortEE}_{flankingES}-{flankingEE}"))
  )

## splicing cases
splicing_case <- psi_comb %>% 
  count(Preference, splicing_case)

counts_psi_comb <- psi_comb %>% 
  dplyr::count(splicing_case, Preference)  %>%
  mutate(splicing_case = factor(splicing_case, levels = c("RI", "A5SS", "A3SS", "MXE", "SE")))

table(splicing_case$Preference, splicing_case$splicing_case)

# Create the side-by-side bar plot with custom colors
pdf(de_events_plot_path, height = 3, width = 5.5)
ggplot(counts_psi_comb, aes(x = Preference, y = n, fill = splicing_case)) +
  geom_bar(stat = "identity", position = position_stack()) +
  scale_fill_manual(name = "Splicing Case", values = c(SE = "#0C7BDC", MXE = "orchid3", RI = "#00A5A5", A3SS = "#FFC20A", A5SS = "#E04C2F"))+ 
  labs(x = "Preference",
       y = "Differential Splice Events") +
  theme_Publication() + 
  coord_flip()+
  ylim(c(0,5000)) +
  theme(legend.position = c(0.4, 1.04), # shift legend to the left to fit
        legend.direction = "horizontal") +
  guides(fill = guide_legend(reverse = TRUE))
dev.off()

# annotate significant splice events for TSG/Oncogene
annots <- read_tsv(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) %>%
  dplyr::rename(geneSymbol = Gene_Symbol) %>%
  mutate(annotation = case_when(grepl("Oncogene|TumorSuppressorGene", type) ~ "Onco_TSG",
                              type %in% c("Kinase", "Kinase, CosmicCensus", "Kinase, TranscriptionFactor", "CosmicCensus, Kinase CosmicCensus, TranscriptionFactor") ~ "Kinase",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(annotation)) %>% 
  select(geneSymbol, annotation) %>%
  unique()

# export significant splice events
psi_comb_select <- psi_comb %>% 
  dplyr::select(splicing_case, geneSymbol, PValue, FDR, IncLevelDifference, exonStart_0base, exonEnd, 
                           "1stExonStart_0base",'1stExonEnd',
                           '2ndExonStart_0base','2ndExonEnd','riExonStart_0base', 'riExonEnd',"upstreamES", 
                           "upstreamEE","downstreamES","downstreamEE",
                           "longExonStart_0base","longExonEnd",
                           "shortES","shortEE",
                           "flankingES","flankingEE","Preference") %>%
  left_join(annots)

write_tsv(psi_comb_select, 
          file=file.path(results_dir,"splice-events-significant.tsv"))

