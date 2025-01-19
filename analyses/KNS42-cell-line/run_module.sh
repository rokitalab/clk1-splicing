#!/bin/sh

## DepMap portal cell-lines analyses
Rscript --vanilla 01-plot-and-process-depmap-data.R

## cell prolif assay
Rscript --vanilla 02-plot_cell-viability-assay-res.R

## PCR assay
Rscript --vanilla 03-plot-qPCR-results.R

## prioritization, depmap and morph targets
# first download the CRISPR file into data
FILE_URL="https://plus.figshare.com/ndownloader/files/51064667"
OUTPUT_FILE="../../data/CRISPRGeneEffect.csv"

# Download the file using wget or curl
echo "Downloading file from $FILE_URL..."

# Using wget
wget -O "$OUTPUT_FILE" "$FILE_URL"

if [ $? -eq 0 ]; then
    echo "Download complete: $OUTPUT_FILE"
else
    echo "Error: Failed to download the file."
    exit 1
fi

Rscript --vanilla 04-prioritization-depmap-morph.R
