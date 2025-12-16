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
data/histologies-plot-group.tsv
data/splice-events-rmats.tsv.gz
```

*Output files:*
```
results/splicing_events.total.neg.intersectunip.ggplot.txt # All overlapping events < -2 z-scores, "skipping events"
results/splicing_events.total.pos.intersectunip.ggplot.txt # All overlapping events > 2 z-scores, "inclusion events"
results/kinases-functional_sites.tsv # Functional sites in kinases
results/splicing-factor-kinases-functional_sites.tsv # Functional sites in kinase splicing factors
```

## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `01-get-uniprot.sh` downloads relevant uniprot files
* `01-extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `02-run_bedtools_intersect.sh` runs bedtools intersect to find exon coordinates corresponding to Uniprot sites
* `03-format_for_ggplot.sh` formats and appends file into table for plotting
* `04-functional-sites-kinases.R` subsets the functional site results down to kinases
* `05-plot-splice-patterns` generates plots for visualizing splicing event types into `plots` folder

## Directory structure
```
.
├── 00-get-uniprot.sh
├── 01-extract_recurrent_splicing_events_cluster.pl
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
│   ├── kinases-ora-plot.pdf
│   └── splicing_pattern_plot.pdf
├── results
│   ├── kinases-functional_sites.tsv
│   ├── line_count.txt
│   ├── splice_events.diff.A3SS.txt
│   ├── splice_events.diff.A5SS.txt
│   ├── splice_events.diff.RI.txt
│   ├── splice_events.diff.SE.txt
│   ├── splicing-factor-kinases-functional_sites.tsv
│   ├── splicing_events.A3SS.total.neg.bed
│   ├── splicing_events.A3SS.total.pos.bed
│   ├── splicing_events.A5SS.total.neg.bed
│   ├── splicing_events.A5SS.total.pos.bed
│   ├── splicing_events.RI.total.neg.bed
│   ├── splicing_events.RI.total.pos.bed
│   ├── splicing_events.SE.total.neg.bed
│   ├── splicing_events.SE.total.pos.bed
│   ├── splicing_events.total.neg.intersectunip.ggplot.txt
│   └── splicing_events.total.pos.intersectunip.ggplot.txt
└── run_module.sh

```
