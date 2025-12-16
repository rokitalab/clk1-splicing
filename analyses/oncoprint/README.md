# Oncoprint

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess the cohort's mutation status

## Usage
### How to Run:
```
bash run-oncoprint.sh
```

#### Input:
* MAF - T/N: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz``` <br>
* MAF - T only: ```./data/snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz``` <br>
* clinical info: ```./data/histologies-plot-group.tsv```  <br>
* genes of interest: ```input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv``` <br>
* mutation colors: ```input/mutation-colors.R``` <br>
* splicing events: ```./data/splice-events-rmats.tsv.gz```  <br>
* tumor mutation burden: ```./input/snv-mutation-tmb-coding.tsv``` <br>

## Scripts
* `01-oncoprint.R` generates oncoprint with gender, molecular subtype, CNS region and mutation status information across the cohort
* `02-SF-mutations.R` extracts mutations in HUGO spliceosome or all SFs for cluster 7 and all samples.

## Directory Structure
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
│   ├── hugo-all-samples.maf
│   ├── hugo-cluster7-samples.maf
│   ├── onco_matrix.tsv
│   ├── oncoprint_sample_metadata.tsv
│   ├── sf-all-samples.maf
│   └── sf-cluster7-samples.maf
└── run-oncoprint.sh
```
