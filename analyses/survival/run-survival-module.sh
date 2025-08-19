#!/bin/bash

set -e
set -o pipefail

# run prepare-survival script
Rscript -e "rmarkdown::render('01-prepare-survival.Rmd')"

# survival analyses by splicing cluster
Rscript -e "rmarkdown::render('02-survival_by_cluster.Rmd')"

# run survival by cluster assignment and SBI or KEGG spliceosome GSVA
Rscript -e "rmarkdown::render('03-survival_by_cluster_sbi_spliceosome_gsva.Rmd')"

# Run splicing cluster assignment, splicing burden, gsva, CLK1
Rscript -e "rmarkdown::render('04-survival_by_cluster_clk1.Rmd')"