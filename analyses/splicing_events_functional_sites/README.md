# Splicing Overlapping Functional Sites

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing events that result in loss/gain of functional sites (as defined by Uniprot) amongst identified differential splicing events

## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run shell script to make tables and subsequent plots below**
```
./run_module.sh
```

*Input files:*
```
data/histologies.tsv
data/splice-events-rmats.tsv.gz
```

*Output files:*
```
results/splice_events.diff.SE.HGG.txt
results/splicing_events.SE.total.HGG.pos.bed
results/splicing_events.SE.total.HGG.neg.bed
results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt

```

## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `01-get-uniprot.sh` downloads relevant uniprot files
* `01-extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `02-run_bedtools_intersect.sh` runs bedtools intersect to find exon coordinates corresponding to Uniprot sites
* `03-format_for_ggplot.sh` formats and appends file into table for plotting
* `04-plot_splicing_across_functional_sites.R` generates ggplot violin plots of average dPSI per event identidied overlapping a functional site, outputting to `plots/*png`
* `05-plot-splice-patterns` generates plots for visualizing splicing event types into `plots` folder

## Directory structure
```
.
├── 00-get-uniprot.sh
├── 01-extract_recurrent_splicing_events_cluster.pl
├── 01-extract_recurrent_splicing_events_hgg.pl
├── 02-run_bedtools_intersect.sh
├── 03-format_for_ggplot.pl
├── 04-plot_splicing_across_functional_sites.R
├── 05-plot-splice-patterns.R
├── README.md
├── input
│   ├── UP000005640_9606_disulfid.bed
│   ├── UP000005640_9606_domain.bed
│   ├── UP000005640_9606_mod_res.bed
│   ├── UP000005640_9606_signal.bed
│   └── gene_lists.tsv
├── plots
│   ├── dPSI_across_functional_sites.pdf
│   ├── dPSI_across_functional_sites_kinase.pdf
│   ├── kinases-ora-plot.pdf
│   └── splicing_pattern_plot.pdf
├── results
│   ├── Naqvi-CLK1-workflow.pdf
│   ├── kinases-functional_sites.tsv
│   ├── splice_events.diff.A3SS.txt
│   ├── splice_events.diff.A5SS.txt
│   ├── splice_events.diff.RI.txt
│   ├── splice_events.diff.SE.txt
│   ├── splicing-factor-kinases-functional_sites.tsv
│   ├── splicing_events.SE.total.neg.intersectunip.ggplot.txt
│   └── splicing_events.SE.total.pos.intersectunip.ggplot.txt
└── run_module.sh
```
