#!/bin/sh
#
# Script to run all analysis modules
# Author: Ammar S Naqvi
# Date: 2024/05/10
#


set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

## cohort summary
echo "----------------------------------"
echo "cohort summary"
echo "----------------------------------"
cd ${analyses_dir}/cohort_summary
Rscript 01-generate-cohort-summary-circos-plot.R

## stranded-polyA-assessment
echo "----------------------------------"
echo stranded-polyA-assessment
echo "----------------------------------"
cd ${analyses_dir}/stranded-polyA-assessment
bash run_module.sh

## histology-specific splice events
echo "----------------------------------"
echo "histology-specific splice events"
echo "----------------------------------"
cd ${analyses_dir}/histology-specific-splicing
bash run_module.sh

## splicing index
echo "----------------------------------"
echo "splicing index"
cd ${analyses_dir}/splicing_index
bash run_module.sh

## psi clustering
echo "----------------------------------"
echo "psi clustering"
cd ${analyses_dir}/sample-psi-clustering
bash run_module.sh

## survival
echo "----------------------------------"
echo "survival"
echo "----------------------------------"
cd ${analyses_dir}/survival
bash run-survival-module.sh

## splicing-factor_dysregulation
echo "----------------------------------"
echo "splicing factor dysregulation"
echo "----------------------------------"
cd ${analyses_dir}/splicing-factor_dysregulation
bash run_module.sh

## splicing events functional sites
echo "----------------------------------"
echo "splicing events functional sites"
echo "----------------------------------"
cd ${analyses_dir}/splicing_events_functional_sites
bash run_module.sh

## CLK1 splicing correlations
echo "----------------------------------"
echo "CLK1 splicing correlations"
echo "----------------------------------"
cd ${analyses_dir}/CLK1-splicing_correlations
bash run_module.sh

## oncoprint
echo "----------------------------------"
echo "oncoprint"
echo "----------------------------------"
cd ${analyses_dir}/oncoprint
bash run-oncoprint.sh

## long-read-CLK1-validation
echo "----------------------------------"
echo "long-read-CLK1-validation"
cd ${analyses_dir}/long-read-CLK1-validation
bash run_module.sh

## KNS42-cell-line
echo "----------------------------------"
echo "KNS42-cell-line"
echo "----------------------------------"
cd ${analyses_dir}/KNS42-cell-line
bash run_module.sh

## CLK1-splicing-impact-morpholino
echo "----------------------------------"
echo "CLK1-splicing-impact-morpholino"
echo "----------------------------------"
cd ${analyses_dir}/CLK1-splicing-impact-morpholino
bash run_module.sh
