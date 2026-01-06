################################################################################
# 03-plot-qPCR-results.R
# script that generates barplot of fold change for qRT-PCR
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 03-plot-qPCR-results.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("tidyverse")
  library("rstatix")
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
hprt_ct <- qpcr_fc_df %>%
  dplyr::filter(Primers == "HPRT") %>%
  group_by(Treatment) %>%
  summarise(mean_hprt = mean(CT), .groups = "drop")

qpcr_dct_rep <- qpcr_fc_df %>%
  dplyr::filter(Primers != "HPRT") %>%
  left_join(hprt_ct, by = "Treatment") %>%
  dplyr::mutate(dCT = CT - mean_hprt) %>%
  dplyr::mutate(
    Primers = case_when(
      Primers == "Exons 3-5 set 1" ~ "Exons 3-5 (A)",
      Primers == "Exons 3-5 set 2" ~ "Exons 3-5 (B)",
      TRUE ~ Primers
    )
  )

# -------------------------------
# Compute FC
# -------------------------------
fc <- qpcr_dct_rep %>%
  group_by(Primers, Treatment) %>%
  summarise(mean_dCT = mean(dCT, na.rm = TRUE),
            se_dCT   = sd(dCT, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  group_by(Primers) %>%
  mutate(
    dCT_ref     = mean_dCT[Treatment == "Control"][1],
    se_ref      = se_dCT[Treatment == "Control"][1],
    
    ddCT        = mean_dCT - dCT_ref,
    se_ddCT     = sqrt(se_dCT^2 + se_ref^2),
    
    fold_change = 2^(-ddCT),
    SE_FC       = fold_change * log(2) * se_ddCT
  ) %>%
  ungroup()

fc$Treatment <- factor(fc$Treatment, levels = c("Control","1uM","5uM","10uM"))

# -------------------------------
# T test
# -------------------------------
pvals_vs_ctrl <- qpcr_dct_rep %>%
  group_by(Primers) %>%
  t_test(dCT ~ Treatment, ref.group = "Control") %>%
  adjust_pvalue(method = "BH") %>%
  dplyr::rename(Treatment = group2)

sign_test_df <- fc %>%
  left_join(
    pvals_vs_ctrl %>% select(Primers, Treatment, p.adj.signif),
    by = c("Primers", "Treatment")
  ) %>%
  group_by(Primers, Treatment) %>%
  mutate(
    p_label = if_else(Treatment == "Control", " ", p.adj.signif),
    y_pos = fold_change + SE_FC + 0.08
  ) %>%
  ungroup()

sign_test_df$Treatment <- factor(sign_test_df$Treatment, levels = c("Control","1uM","5uM","10uM"))


# -------------------------------
# Plot
# -------------------------------
qpcr_plot <- ggplot(fc, aes(x = Primers, y = fold_change, fill = Treatment)) + 
  geom_bar(position = position_dodge(width = 0.9), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = fold_change - SE_FC, ymax = fold_change + SE_FC),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(
    data = sign_test_df,
    aes(x = Primers, y = y_pos, label = p_label, group = Treatment),
    position = position_dodge(width = 0.9),
    vjust = 0,
    fontface = "bold",
    size = 4,
    inherit.aes = FALSE
  ) +
  ylab("Fold Change") + 
  xlab(expression(bold(bolditalic("CLK1")~"exon-exon junction"))) +
  theme_Publication() +
  scale_fill_manual(values = c("lightgrey", "lightblue", "#0C7BDC", "blue3")) +
  guides(fill = guide_legend(title = "Morpholino\nTreatment"))

# Save plot
pdf(file_plot, width = 6.5, height = 4)
print(qpcr_plot)
dev.off()
