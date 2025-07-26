# Oncoprint

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess CLK1 exon 4 PSI levels in the context of mutations

## Usage
### How to Run:
```
bash run-oncoprint.sh
```

#### Input:
* MAF - T/N: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz``` <br>
* MAF - T only: ```./data/snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz``` <br>
* clinical info: ```./data/histologies.tsv```  <br>
* genes of interest: ```input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv``` <br>
* independent specimen file (RNA): ```./data/independent-specimens.rnaseqpanel.primary.tsv```  <br>
* mutation colors: ```input/mutation-colors.R``` <br>
* splicing events: ```./data/splice-events-rmats.tsv.gz```  <br>
* tumor mutation burden: ```./input/snv-mutation-tmb-coding.tsv``` <br>

#### Output:
```plots/oncoprint.pdf``` <br>

## Scripts
* `01-oncoprint.R` generates oncoprint with mutation frequencies with CLK1 exon 4 PSI, gender, molecular subtype, CNS region and mutation status information across pediatric HGGs, as well as enrichment of CLK1 high/low tumors by gene alteration 
* `02-SF-mutations.R` extracts mutations in HUGO spliceosome or all SFs for cluster 6 or all samples.

## Directory Structure
```
.
├── 01-oncoprint.R
├── 02-SF-mutations.R
├── README.md
├── input
│   ├── hgnc-symbol-check.csv
│   ├── mutation-colors.R
│   ├── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
│   └── snv-mutation-tmb-coding.tsv
├── plots
│   ├── oncoprint-clk1-psi.pdf
│   ├── oncoprint-cluster.pdf
│   ├── oncoprint-hist-cluster.pdf
│   ├── oncoprint-hist-cluster1.pdf
│   ├── oncoprint-hist-cluster10.pdf
│   ├── oncoprint-hist-cluster11.pdf
│   ├── oncoprint-hist-cluster2.pdf
│   ├── oncoprint-hist-cluster3.pdf
│   ├── oncoprint-hist-cluster4.pdf
│   ├── oncoprint-hist-cluster5.pdf
│   ├── oncoprint-hist-cluster6.pdf
│   ├── oncoprint-hist-cluster7.pdf
│   ├── oncoprint-hist-cluster8.pdf
│   └── oncoprint-hist-cluster9.pdf
├── results
│   ├── clk1_high_low_mutation_counts.tsv
│   ├── hgg_cluster_mutation_count.tsv
│   ├── hgg_cluster_mutation_freq.tsv
│   ├── hugo-all-samples.maf
│   ├── hugo-cluster6-samples.maf
│   ├── onco_matrix.tsv
│   ├── onco_matrix_SFs.tsv
│   ├── oncoprint_sample_metadata.tsv
│   ├── oncoprint_sample_metadata_SFs.tsv
│   ├── sbi_high_low_mutation_counts_SFs.tsv
│   ├── sf-all-samples.maf
│   └── sf-cluster6-samples.maf
└── run-oncoprint.sh
```
