# KNS42 Cell Line
This module investigates the dependency of CLK1 in pediatric HGG cell line KNS42 relative to all other brain tumor cell lines using DepMap portal data, and generates visualizations of KNS42 morpholino treated cells and functional validations (experimental) of CLK1 exon 4 splicing on cell proliferation.  

## Usage
To run analysis module from the command line:
```
bash run_module.sh
```

### Input
```
input/CLK1_CRISPR_depmap_score.csv
input/OmicsDefaultModelProfiles.csv
input/cell_prolif_res.tsv
input/qpcr-raw-triplicates.tsv
```

## Folder Content
```01-plot-and-process-depmap-data.R``` is script plotting CLK1 dependency scores and CLK1 transcript expression for brain tumor cell lines available in the DepMap Portal<br>
```02-plot_cell-proliferation-assay-res.R``` plots the results from the cell proliferation assay of treated vs ctrl/untreated cells<br>
```03-plot-qPCR-results.R``` visualizes qRT-PCR results
```04-prioritization-depmap-morph.R``` investigates the prioritization of functional splice events, clk1 targets, and depmap dependencies

Current DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions
DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. Dataset. https://doi.org/10.25452/figshare.plus.27993248.v1


## Directory Structure
```
.
├── 01-plot-and-process-depmap-data.R
├── 02-plot_cell-viability-assay-res.R
├── 03-plot-qPCR-results.R
├── 04-prioritization-depmap-morph.R
├── README.md
├── input
│   ├── 2023-03-22 162604_JLRmod_20240218.xls
│   ├── CLK1_CRISPR_depmap_score.csv
│   ├── CRISPRGeneEffect.csv
│   ├── Model.csv
│   ├── OmicsDefaultModelProfiles.csv
│   ├── cell_prolif_res.tsv
│   ├── qpcr-results-raw-ct.csv
│   └── tam_etal_clk1_targets.txt
├── plots
│   ├── cell_viability-barplot.pdf
│   ├── depmap_score_CLK1_vs_score_KNS42.pdf
│   ├── depmap_score_all_cell_lines.pdf
│   ├── depmap_score_cns_cell_lines.pdf
│   └── qPCR-morp.pdf
├── results
│   ├── clk1_consensus_targets.tsv
│   ├── mean_ped_glioma_crispr_scores_func_kinase_splice_events.csv
│   └── splice_genes_functional.tsv
└── run_module.sh
```