#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

## histology input file (column orders important)
#input_file="../cohort_summary/results/histologies-plot-group.tsv"
#primary_specimens="../../data/independent-specimens.rnaseqpanel.primary.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
cluster_file="../sample-psi-clustering/results/sample-cluster-metadata-top-5000-events-stranded.tsv"

echo "input files:" $cluster_file ;
echo $rmats_file ;

## Download uniprot files
bash 00-get-uniprot.sh

## Process rMATS files given clusters.
perl 01-extract_recurrent_splicing_events_cluster.pl $cluster_file $rmats_file RI
perl 01-extract_recurrent_splicing_events_cluster.pl $cluster_file $rmats_file A3SS
perl 01-extract_recurrent_splicing_events_cluster.pl $cluster_file $rmats_file A5SS
perl 01-extract_recurrent_splicing_events_cluster.pl $cluster_file $rmats_file SE

echo "bedtools intersect...";
bash 02-run_bedtools_intersect.sh

echo "make tab for ggplot ...";
perl 03-format_for_ggplot.pl

## make plots
echo "make plots ...";
Rscript 04-plot_splicing_across_functional_sites.R

## make plots
echo "plot splice patterns";
Rscript --vanilla 05-plot-splice-patterns.R

##rm intermediatery files
rm results/splicing_events*.wo.txt
rm results/*bed
