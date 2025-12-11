#!/bin/sh

## histology and rmats file
hist_file="../cohort_summary/results/histologies-plot-group.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"

## N≥2 histology specific splicing based on splicing index computations for all events
perl 01-generate_hist_spec_events_tab.pl $hist_file $rmats_file

# upset plot from the total recurrent events file
Rscript --vanilla 02-plot_histology-specific_splicing_events.R

# plot normalized unique events by histology from the total recurrent events file
Rscript --vanilla 03-plot-histology-specific-norm-events.R

rm Rplots.pdf
