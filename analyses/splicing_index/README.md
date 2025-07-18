# Splicing Index

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to compute the splicing index of each tumor (proportion of mis-spliced events)

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files (`data` folder, or TMB in `input` folder):
```
histologies.tsv
splice-events-rmats.tsv.gz
independent-specimens.rnaseqpanel.primary-plus.tsv
independent-specimens.rnaseqpanel.primary-plus.tsv
snv-mutation-tmb-coding.tsv
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-generate_splicing-index_and_diff-events_table.pl` processes rMATs output and computes splicing burden index per sample stratified by splicing case
* `02-plot_splicing_burden_index.R` takes splicing index burden values and generates CDF plot categorized by histologies and denotes median values
* `03-identify_and_plot_histologies_by_SBI.R` identifies high vs low SBI tumors categorized by histologies
* `04-plot_total-splicing-cases.R` plots total splicing cases across samples  
* `05-plot-tmb-vs-sbi.R` plots splicing burden vs TMB by mutation status and cancer group
* `06-plot-gsva-score-vs-sbi.R` plots splicing burden vs GSVA score for the splicosome by histology or cluster
* `07-plot-sbi-high-low.R` plots splicing burden high vs low samples by histology

## Directory structure
```
.
├── 01-generate_splicing-index_and_diff-events_table.pl
├── 02-plot_splicing_burden_index.R
├── 03-identify_and_plot_histologies_by_SBI.R
├── 04-plot_total-splicing-cases.R
├── 05-plot-tmb-vs-sbi.R
├── 06-plot-gsva-score-vs-sbi.R
├── 07-plot-sbi-high-low.R
├── README.md
├── input
│   └── snv-mutation-tmb-coding.tsv
├── plots
│   ├── boxplot_sbi-tmb-by-cg.pdf
│   ├── boxplot_sbi-tmb-by-mutation-status.pdf
│   ├── corplot-sbi-vs-gsva-spliceosome-by-cluster.pdf
│   ├── corplot-sbi-vs-gsva-spliceosome-by-hist.pdf
│   ├── corplot-sbi-vs-gsva-spliceosome.pdf
│   ├── corplot_sbi-tmb-by-cg.pdf
│   ├── corplot_sbi-tmb.pdf
│   ├── hist_by_sbi_level_barplot.pdf
│   ├── sbi-plot-A3SS-boxplot.pdf
│   ├── sbi-plot-A3SS.pdf
│   ├── sbi-plot-A5SS-boxplot.pdf
│   ├── sbi-plot-A5SS.pdf
│   ├── sbi-plot-RI-boxplot.pdf
│   ├── sbi-plot-RI.pdf
│   ├── sbi-plot-SE-boxplot.pdf
│   ├── sbi-plot-SE.pdf
│   ├── scale-sbi-plot-SE-boxplot.pdf
│   └── splice-types.pdf
├── results
│   ├── splice_events.diff.A3SS.txt.gz
│   ├── splice_events.diff.A5SS.txt.gz
│   ├── splice_events.diff.RI.txt.gz
│   ├── splice_events.diff.SE.txt.gz
│   ├── splicing_index.A3SS.txt
│   ├── splicing_index.A5SS.txt
│   ├── splicing_index.RI.txt
│   └── splicing_index.SE.txt
└── run_module.sh
```
