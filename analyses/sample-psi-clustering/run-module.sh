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