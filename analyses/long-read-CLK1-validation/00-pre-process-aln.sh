#!/bin/sh

# make results folder
mkdir -p results

cat input/7316_1763.CLK1.aln |  grep -P "chr2\t" | grep -P "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > results/7316_1763.CLK1.processed.txt
cat input/7316_1769.CLK1.aln | grep -P "chr2\t" | grep -P "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > results/7316_1769.CLK1.processed.txt
cat input/KNS42.CLK1.aln | grep -P "chr2\t" | grep -P "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > results/KNS42.CLK1.processed.txt
