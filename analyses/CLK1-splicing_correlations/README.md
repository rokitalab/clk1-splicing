# CLK1 Splicing Correlations

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza), Patricia Sullivan (@pj-sullivan)

The purpose of this module is to correlate CLK1 Exon 4 splicing with splicing burden, normals, RNA expression, and proteomics

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `00-subset-CLK1-from-rMATs.R` subset rMATs output for CLK1
* `01-plot_CLK1_highExon4_vs_lowExon4_and_SBI.R` correlate high vs low levels of overall splicing burden with CLK1 exon 4 inclusion levels
* `02-plot_splicing_vs_expr.R` correlate CLK1 exon 4 inclusion levels with RNA expression, including CLK1, SRSF1, SRSF2, SRSF10
* `03-CLK-SRSF-expr-correlations.R` generates scatter plots of CLK and SRSF transcript abundance vs. CLK1 exon 4 PSI and transcript abundance
* `04-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R` generates summary heatmap of spearman correlation coefficients between CLK1 PSI and transcript abundance against CLK and SRSF RNA, total protein, and phosphoprotein expression
* `05-plot-Ex4-PSI-all-tumors.R` plot CLK1 Ex 4 PSIs across tumors
* `06-plot-kinase-PSI-variances-tumors.R` plot CLK1 Ex 4 and all splicing factor kinase PSIs that are functional and differential
* `07-plot-clk1-ex4-transcripts-normals.R` plot CLK1 Ex 4 transcripts proportion in tumors and control samples
* `08-plot-clk1-ex4-transcripts-by-age.R` plots CLK1 exon 4 PSI variations across tumors by age
* `09-volcano-plot-high-vs-low-Ex4-PSI.R` compare high vs low Exon 4 PSI samples

## Directory structure
```
.
├── 00-subset-CLK1-from-rMATs.R
├── 01-plot_CLK1_highExon4_vs_lowExon4_and_SBI.R
├── 02-plot_splicing_vs_expr.R
├── 03-CLK-SRSF-expr-correlations.R
├── 04-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R
├── 05-plot-Ex4-PSI-all-tumors.R
├── 06-plot-kinase-PSI-variances-tumors.R
├── 07-plot-clk1-ex4-transcripts-normals.R
├── 08-plot-clk1-ex4-transcripts-by-age.R
├── 09-volcano-plot-high-vs-low-Ex4-PSI.R
├── README.md
├── archive
│   ├── 08-CLK1-impact-NF1-splicing.Rmd
│   ├── 08-CLK1-impact-NF1-splicing.html
│   ├── 09-clk1-nf1-protein-correlations.R
│   └── 10-clk1-nf1-single-sample-heatmap.R
├── input
│   ├── BA_KFWTGZPC_20250506.filtered.SE.MATS.JC.txt
│   ├── BA_KFWTGZPC_20250506.rsem.isoforms.results.gz
│   ├── CPTAC3-pbt.xls
│   └── gtex-samples-by-age.tsv
├── plots
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK1-Ex4-PSI-clusters.pdf
│   ├── CLK1-Ex4-PSI-cohort.pdf
│   ├── CLK1-Ex4-range-across-ctrls.pdf
│   ├── CLK1-Ex4-range-across.pdf
│   ├── CLK1-Ex4-sdev-across.pdf
│   ├── CLK1-psi-expr-correlation-heatmap.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_201_exp_DMG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_201_exp_HGG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_exp_DMG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_exp_HGG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_201_exp_DMG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_201_exp_HGG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_exp_DMG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_exp_HGG.pdf
│   ├── CLK1_exp_v_SRPK1_exp_all_hgg.pdf
│   ├── CLK1_exp_v_SRPK1_exp_midline_hgg.pdf
│   ├── CLK1_exp_v_SRPK1_exp_other_hgg.pdf
│   ├── CLK1_exp_vs_CLK1_psi_full_cohort.pdf
│   ├── CLK1_exp_vs_CLK1_psi_hgg_all.pdf
│   ├── CLK1_exp_vs_CLK1_psi_hgg_midline.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── DMGs_volcano_plot.pdf
│   ├── Other HGGs_volcano_plot.pdf
│   ├── PBTA_volcano_plot.pdf
│   ├── PSI-range-kinses.pdf
│   ├── SRPK_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_full_cohort.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_hgg_all.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_hgg_midline.pdf
│   ├── SRSF11_exp_vs_CLK1_psi_full_cohort.pdf
│   ├── SRSF11_exp_vs_CLK1_psi_hgg_all.pdf
│   ├── SRSF11_exp_vs_CLK1_psi_hgg_midline.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_full_cohort.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_hgg_all.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_hgg_midline.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_full_cohort.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_hgg_all.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_hgg_midline.pdf
│   ├── SRSF_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── all HGGs_volcano_plot.pdf
│   ├── all_SBI_high_vs_low_CLK1_stranded.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_stranded.pdf
│   ├── clk1_ex4-tpm-ctrls-summary.pdf
│   ├── clk1ex4-tpm-tumor-age-bin-perm-test.pdf
│   ├── clk1ex4-tpm-tumor-age-bin-wc-test.pdf
│   ├── dmg_SBI_high_vs_low_CLK1_stranded.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-DMGs.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-Other HGGs.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-PBTA.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-all HGGs.pdf
│   └── other_hgg_SBI_high_vs_low_CLK1_stranded.pdf
├── results
│   ├── DMGs_gene_sign_list.tsv
│   ├── Other HGGs_gene_sign_list.tsv
│   ├── PBTA_gene_sign_list.tsv
│   ├── all HGGs_gene_sign_list.tsv
│   ├── clk1-exon4-psi-hgg.tsv
│   ├── clk1-exon4-psi-normals-stats.tsv
│   ├── clk1-exon4-psi.tsv
│   ├── clk1-splice-events-rmats.tsv
│   └── hgg-dmg-clk-srsf-expression-phosphorylation.tsv
├── run_module.sh
└── util
    └── function-create-scatter-plot.R
```
