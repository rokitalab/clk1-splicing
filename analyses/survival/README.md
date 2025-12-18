# Survival by splicing burden

Module authors: Ryan Corbett (@rjcorb), Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess survival in CNS tumor pediatric patients by cluster, overall splicing burden index, spliceosome GVSA score, and CLK1 TPM and exon 4 PSI.

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run-survival-module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
../splicing_index/results/splicing_index.total.txt
```

## Folder content
* `run-survival-module.sh` shell script to run analysis
* `01-prepare-survival.Rmd` Merge splicing data with relevant histology data and define survival features
* `02-survival_by_cluster.Rmd` Assess survival by splicing cluster assignment in all histologies, LGG, and HGG
* `03-survival_by_cluster_sbi_spliceosome_gsva.Rmd` Plot survival by cluster and SI or GSVA status
* `04-survival_by_cluster_clk1.Rmd` Plot survival by cluster and CLK1 TPM and PSI
