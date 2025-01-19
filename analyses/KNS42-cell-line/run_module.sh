#!/bin/sh

## DepMap portal cell-lines analyses
Rscript --vanilla 01-plot-and-process-depmap-data.R

## cell prolif assay
Rscript --vanilla 02-plot_cell-viability-assay-res.R

## PCR assay
Rscript --vanilla 03-plot-qPCR-results.R

## prioritization, depmap and morph targets
Rscript --vanilla04-prioritization-depmap-morph.R
