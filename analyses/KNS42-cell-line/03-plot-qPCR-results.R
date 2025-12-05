################################################################################
# 03-plot-qPCR-res.R
# script that generates barplot of fold change for qRT-PCR
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 03-plot-qPCR-res.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("tidyverse")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "KNS42-cell-line")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir,"theme_for_plots.R"))

## output for plot
file_plot = file.path(plots_dir,"qPCR-morp.pdf")

qpcr_res_file = file.path(input_dir,"qpcr-results-raw-ct.csv")
qpcr_fc_df <- read_csv(qpcr_res_file) 

# -------------------------------
# Compute per-replicate FC
# -------------------------------

# Get HPRT per replicate for dCT
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------
# Compute per-replicate FC
# -------------------------------
hprt_ct <- qpcr_fc_df %>%
  filter(Primers == "HPRT") %>%
  select(Treatment, CT)

qpcr_dct_rep <- qpcr_fc_df %>%
  filter(Primers != "HPRT") %>%
  left_join(hprt_ct, by = "Treatment", suffix = c("", "_HPRT")) %>%
  mutate(dCT = CT - CT_HPRT)

control_dct <- qpcr_dct_rep %>%
  filter(Treatment == "Control") %>%
  select(Primers, dCT) %>%
  rename(control_dCT = dCT)

qpcr_ddct_rep <- qpcr_dct_rep %>%
  left_join(control_dct, by = "Primers") %>%
  mutate(
    ddCT = dCT - control_dCT,
    FC = 2^-ddCT
  ) %>%
  mutate(
    Primers = case_when(
      Primers == "Exons 3-5 set 1" ~ "Exons 3-5 (A)",
      Primers == "Exons 3-5 set 2" ~ "Exons 3-5 (B)",
      TRUE ~ Primers
    )
  )

qpcr_ddct_rep$Treatment <- factor(qpcr_ddct_rep$Treatment, levels = c("Control","1uM","5uM","10uM"))

# -------------------------------
# Summary for plotting (mean + SE)
# -------------------------------
qpcr_summary <- qpcr_ddct_rep %>%
  group_by(Primers, Treatment) %>%
  summarise(
    FC_mean = mean(FC),
    SE_FC   = sd(FC)/sqrt(n()),
    .groups = "drop"
  )

# -------------------------------
# Paired Wilcoxon test: 10uM vs 1uM per primer
# -------------------------------
sign_test_df <- qpcr_ddct_rep %>%
  filter(Treatment %in% c("1uM","10uM")) %>%
  group_by(Primers) %>%
  summarise(
    p_val = wilcox.test(FC[Treatment=="10uM"], FC[Treatment=="1uM"], paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )


# -------------------------------
# Plot
# -------------------------------
qpcr_plot <- ggplot(qpcr_summary, aes(x = Primers, y = FC_mean, fill = Treatment)) + 
  geom_bar(position = position_dodge(width = 0.9), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = FC_mean - SE_FC, ymax = FC_mean + SE_FC),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(data = sign_test_df,
            aes(x = Primers, y = max(qpcr_summary$FC_mean) + 0.15, label = p_label),
            inherit.aes = FALSE, vjust = 0, fontface = "bold") +
  ylab("Fold Change") + 
  xlab(expression(bold(bolditalic("CLK1")~"exon-exon junction"))) +
  theme_Publication() +
  scale_fill_manual(values = c("lightgrey", "lightblue", "#0C7BDC", "blue3")) +
  guides(fill = guide_legend(title = "Morpholino\nTreatment")) +
  ylim(c(0, max(qpcr_summary$FC_mean) + 0.3))

# Save plot
pdf(file_plot, width = 6.5, height = 4)
print(qpcr_plot)
dev.off()
