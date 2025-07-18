# Splicing Factor Dysregulation

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to identify known splicing factors [(ensemble-annotated)](https://genome.cshlp.org/content/26/6/732.long) that are dysregulated in high splicing burden tumors, potentially explaining the global levels of splicing dysregulation.

The input gene list is taken from supplementary table 2 of the publication linked above. The HGNC symbol was taken for all genes encoding a known or putative RNA binding protein (RBP).

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-plot-diffExp_highlowSBI.R` performs differential gene expression on high vs low SBI tumors
* `02-plot_SFs_rna_vs_prot.R` generates heatmap of dysregulated splicing factors showing RNA and protein levels across brain tumors  
* `03-plot_CLK_fam_expression.R` generates boxplot of CLK1,CLK2,CLK3,CLK4, and SRPK1 expression (TPM)
* `04-plot_sbi_vs_splicing_factor_expr.R` plot splicing burden vs. splicing factor gene expression by cluster

## Directory structure
```
.
├── 01-plot-diffExp_highlowSBI.R
├── 02-plot_SFs_rna_vs_prot.R
├── 03-plot_CLK_fam_expression.R
├── 04-plot_sbi_vs_splicing_factor_expr.R
├── input
│   ├── CPTAC3-pbt_SF_family.xls
│   ├── CPTAC3-pbt_SFs.xls
│   └── splicing_factors.txt
├── plots
│   ├── all_hgg-barplot-sbi-SFs.pdf
│   ├── all_hgg-CLK1-tpms-CLK-SPRK1-kinases.pdf
│   ├── all_hgg-enhancedVolcano-sbi.pdf
│   ├── all_hgg-SF_RNA_vs_protein_levels_heatmap.pdf
│   ├── cluster6-barplot-sbi-SFs.pdf
│   ├── cluster6-CLK1-tpms-CLK-SPRK1-kinases.pdf
│   ├── cluster6-enhancedVolcano-sbi.pdf
│   ├── cluster6-SF_RNA_vs_protein_levels_heatmap.pdf
│   ├── dmg-barplot-sbi-SFs.pdf
│   ├── dmg-enhancedVolcano-sbi.pdf
│   ├── other_hgg-barplot-sbi-SFs.pdf
│   ├── other_hgg-enhancedVolcano-sbi.pdf
│   ├── sbi-sf-correlation-heatmap-byCluster.pdf
│   ├── sbi-vs-clk1-ex4-psi-cluster6.pdf
│   ├── sbi-vs-clk1-ex4-psi-other-clusters.pdf
│   ├── sbi-vs-clk1-tpm-byHist.pdf
│   ├── sbi-vs-clk1-tpm-cluster6.pdf
│   └── sbi-vs-clk1-tpm-other-clusters.pdf
├── README.md
├── results
│   ├── all_hgg-diffSFs_sig_genes.txt
│   ├── cluster6-diffSFs_sig_genes.txt
│   ├── dmg-diffSFs_sig_genes.txt
│   └── other_hgg-diffSFs_sig_genes.txt
└── run_module.sh
```
