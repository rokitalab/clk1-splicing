#!/bin/sh

## plot oncoprint
echo "plotting oncoprint"
Rscript --vanilla 01-oncoprint.R

#Commenting out as it uses a file not available anymore
#Rscript --vanilla 02-oncoprint-SFs.R
