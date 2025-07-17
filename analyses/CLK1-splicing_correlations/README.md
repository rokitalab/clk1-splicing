# CLK1 Splicing Correlations

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to correlate CLK1 Exon 4 splicing with splicing
burden, RNA expression and proteomics

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
./data/clk1-splice-events-rmats.tsv.gz
./data/gene-counts-rsem-expected_count-collapsed.rds
./analyses/CLK1-splicing_correlations/input/CPTAC3-pbt.xls
../splicing_index/results/splicing_index.SE.txt
../splicing-factor_dysregulation/input/splicing_factors.txt
./data/cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
./data/cptac-protein-imputed-prot-expression-abundance.tsv.gz
./data/hope-protein-imputed-phospho-expression-abundance.tsv.gz
./data/hope-protein-imputed-prot-expression-abundance.tsv.gz
```

## Folder content
* `run_module.sh` shell script to run analysis
* `00-subset-CLK1-from-rMATs.R` subset rMATs output for CLK1
* `01-plot_highExon4_vs_lowExon4_and_SBI.R` correlate high vs low levels of overall splicing burden with CLK1 exon 4 inclusion levels
* `02-plot_splicing_vs_expr.R` correlate CLK1 exon 4 inclusion levels with RNA expression, including CLK1, SRSF1, SRSF2, SRSF10
* `03-plot_CLK1-Ex4-splicing_vs_SRSF1-expr.R` correlate CLK1 RNA expression with RNA SRSFs, phospho, and total proteomics across brain tumor types
* `04-CLK1_PSI_plots.R` generates stacked barplot of relative inclusion/skipping of exon 4 across midline HGG tumors with differential splicing in CLK1.
* `05-CLK-SRSF-expr-correlations.R` generates scatter plots of CLK and SRSF transcript abundance vs. CLK1 exon 4 PSI and transcript abundance
* `06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R` generates summary heatmap of spearman correlation coefficients between CLK1 PSI and transcript abundance against CLK and SRSF RNA, total protein, and phosphoprotein expression
* `07-plot-clk1ex4-hgg-normals.R` plots CLK1 exon 4 transcript expression for HGGs, GTEX normals, fetal brain, and TGEN pediatric normals.
* `08-plot-Ex4-PSI-all-tumors.R` plot CLK1 Ex 4 PSIs across tumors
* `09-plot-kinase-PSI-variances-tumors.R` plot CLK1 Ex 4 and all splicing factor kinase PSIs that are functional and differential
* `10-plot-clk1-ex4-transcripts-normals.R` plot CLK1 Ex 4 transcripts proportion in tumors and control samples
* `11-plot-clk1-ex4-transcripts-by-age-gtex.R` plots CLK1 exon 4 PSI variations across gtex by age
* `12-plot-clk1-ex4-transcripts-by-age.R` plots CLK1 exon 4 PSI variations across tumors by age
* `13-volcano-plot-high-vs-low-Ex4-PSI.R` compare high vs low Exon 4 PSI samples

## Directory structure
```
.
├── 00-subset-CLK1-from-rMATs.R
├── 01-plot_highExon4_vs_lowExon4_and_SBI.R
├── 02-plot_splicing_vs_expr.R
├── 03-plot_SR-phosp_vs_CLK1-RNA.R
├── 04-CLK1_PSI_plots.R
├── 05-CLK-SRSF-expr-correlations.R
├── 06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R
├── 07-plot-clk1ex4-hgg-normals.R
├── 08-plot-Ex4-PSI-all-tumors.R
├── 09-plot-kinase-PSI-variances-tumors.R
├── 10-plot-clk1-ex4-transcripts-normals.R
├── 11-plot-clk1-ex4-transcripts-by-age-gtex.R
├── 12-plot-clk1-ex4-transcripts-by-age.R
├── 13-volcano-plot-high-vs-low-Ex4-PSI.R
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
│   ├── All_volcano_plot.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
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
│   ├── CLK1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── CLK1_exp_vs_CLK1_psi_midline_hgg.pdf
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
│   ├── PSI-range-kinses.pdf
│   ├── SRPK_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── SR_phos_CLK1_exp_heatmap.pdf
│   ├── all_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── clk1_ex4-tpm-ctrls-summary.pdf
│   ├── clk1_ex4-tpm-phgg-ctrls.pdf
│   ├── clk1ex4-tpm-tumor-age-bin-perm-test.pdf
│   ├── clk1ex4-tpm-tumor-age-bin-wc-test.pdf
│   ├── clk4-tpm-gtex-age-bin-perm-test.pdf
│   ├── clk4-tpm-gtex-age-bin-wc-test.pdf
│   ├── dmg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-All.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-DMGs.pdf
│   ├── high_low_ex4_diff-genes-ora-dotplot-Other HGGs.pdf
│   └── other_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
├── results
│   ├── All_gene_sign_list.tsv
│   ├── DMGs_gene_sign_list.tsv
│   ├── Other HGGs_gene_sign_list.tsv
│   ├── all_hgg-mean_clk1_psi.txt
│   ├── clk1-exon4-proportion.tsv
│   ├── clk1-exon4-psi-hgg.tsv
│   ├── clk1-exon4-psi.tsv
│   ├── clk1-splice-events-rmats.tsv
│   ├── dmg-mean_clk1_psi.txt
│   ├── hgg-dmg-clk-srsf-expression-phosphorylation.tsv
│   ├── nf1-splice-events-rmats.tsv
│   └── other_hgg-mean_clk1_psi.txt
├── run_module.sh
└── util
    └── function-create-scatter-plot.R
```
