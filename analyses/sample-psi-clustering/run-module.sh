#!/bin/bash

set -e
set -o pipefail

if [ -f "results/pbta-splice-event-psis.RDS" ]; then
    echo "Found psi matrix. Proceeding..."
else
    echo "psi matrix file does not exist. Running 00-create-psi-matrix.R..."
    Rscript --vanilla 00-create-psi-matrix.R
fi

# perform clustering
Rscript --vanilla 01-sample_clustering.R

# create geneset rds
Rscript --vanilla 02-create-geneset-rds.R

# assess hist and subtype distribution across clusters
Rscript --vanilla 03-plot-histology-distr-across-clusters.R

# assess pathway enrichment across clusters
Rscript --vanilla 04-diff_pathways_per_cluster.R

# collate spliceosome GSVA data
bash 05-generate-gsva-summary-file.sh