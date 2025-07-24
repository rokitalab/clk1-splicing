#!/bin/bash

set -e
set -o pipefail

# run prepare-survival script
Rscript -e "rmarkdown::render('01-prepare-survival.Rmd')"

# run run-survival scripts
Rscript -e "rmarkdown::render('02-run-survival-SIgroup.Rmd')"

# run plot-survival script
Rscript --vanilla 03-plot-survival.R

# run HGG survival by CLK1 status
Rscript -e "rmarkdown::render('04-survival-hgg-clk1-status.Rmd')"

# survival analyses by splicing cluster
Rscript -e "rmarkdown::render('05-survival_by_cluster.Rmd')"

# run survival by CLK1 status
# Rscript -e "rmarkdown::render('06-run-survival-clk1-status-all.Rmd')"

# plot survival by CLK1 status
# Rscript --vanilla 07-plot-survival-clk1-status-all.R

# run survival by cluster assignment and SBI or KEGG spliceosome GSVA
Rscript -e "rmarkdown::render('08-survival_by_cluster_sbi_spliceosome_gsva.Rmd')"