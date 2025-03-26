#!/bin/bash

set -e
set -o pipefail

if [ -f "results/pbta-splice-event-psis.RDS" ]; then
    echo "Found psi matrix. Proceeding..."
else
    echo "psi matrix file does not exist. Running 00-prepare_rmats.sh..."
    Rscript --vanilla 00-create-psi-matrix.R
fi

# perform clustering
Rscript --vanilla 01-sample_clustering.R

# assess hist and subtype distribution across clusters
Rscript --vanilla 02-plot-histology-distr-across-clusters.R

# assess SBI group distribution across clusters
Rscript --vanilla 03-plot-sbi-with-cluster-mem.R

# assess pathway enrichment across clusters
Rscript --vanilla 04-diff_pathways_per_cluster.R

# collate spliceosome GSVA data
bash 05-generate-spliceosome-summary-file.sh