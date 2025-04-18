################################################################################
# 01-plot_cell-proliferation-assay.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_cell-proliferation-assay.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("vroom")
  library("tidyverse")
  library("ggpubr")
  library("rstatix")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "CLK1-inhibition-assays")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output for plot
file_line_plot = file.path(plots_dir,"cell-prolif.pdf")

## input file
cell_prolif_res_file <- file.path(input_dir,"2024.05.30_KNS-42_6D_cirtuvivint_nogroups.txt")

incucyte_data <- read_tsv(cell_prolif_res_file, skip=10, col_names = TRUE) 

# Convert Date Time column to datetime
incucyte_data <- incucyte_data %>%
  mutate(`Date Time` = as.POSIXct(`Date Time`, format = "%m/%d/%Y %I:%M:%S %p"))

# Reshape data into long format for plotting
incucyte_data_long <- incucyte_data %>%
  pivot_longer(cols = -c(`Date Time`, Elapsed), names_to = "Treatment", values_to = "Measurement") %>%
  mutate(`Date Time` = format(`Date Time`, "%H:%M:%S")) %>%
  dplyr::rename("Time"=`Date Time`) %>%
  filter(!grepl("Std", Treatment)) %>%
  mutate(Treatment = sub(".*(?=CIrtuvivint)", "", Treatment, perl = TRUE),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment),
         Treatment = str_replace(Treatment, "KNS-42\\s*\\d+K / well ", ""),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment))

# Subset the data to only the treatments of interest for stats
filtered_df <- incucyte_data_long %>%
  mutate(Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, 
                                 "Fresh DMEM 100%",
                                 "DMSO vehicle 2%", 
                                 "Cirtuvivint 0.01 µM",
                                 "Cirtuvivint 0.05 µM", 
                                 "Cirtuvivint 0.5 µM", 
                                 "Cirtuvivint 0.1 µM",
                                 "Cirtuvivint 1 µM",
                                 "Cirtuvivint 5 µM",
                                 "Cirtuvivint 10 µM"))

# Compute mean and standard error for each Treatment and Time combination
mean_se_df <- filtered_df %>%
  group_by(Elapsed, Treatment) %>%
  summarise(
    Mean = mean(Measurement, na.rm = TRUE),
    SE = sd(Measurement, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Find the maximum mean value for each elapsed time
max_mean_df <- mean_se_df %>%
  group_by(Elapsed) %>%
  summarise(MaxMean = max(Mean, na.rm = TRUE))

# Plot the means with error bars as a line graph
plot_prolif <- ggplot(mean_se_df, aes(x = as.numeric(Elapsed), y = Mean, group = Treatment, color = Treatment)) +
  geom_line() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.5) +
  geom_point(size = 0.3) +
  labs(x = "Elapsed Time (hours)", y = "KNS-42 Confluence (%)") +
  theme_Publication() 
  

pdf(file_line_plot, width = 7, height = 4)
print(plot_prolif)
dev.off()



