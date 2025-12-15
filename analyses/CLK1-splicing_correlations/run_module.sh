## subset rMATs for CLK1
Rscript --vanilla 00-subset-CLK1-from-rMATs.R

## plot HGGs with high vs low exon 4
Rscript --vanilla 01-plot_CLK1_highExon4_vs_lowExon4_and_SBI.R

## plot correlations of HGG splicing vs expr
Rscript --vanilla 02-plot_splicing_vs_expr.R

# run CLK1-SRSF RNA expression correlation script
Rscript --vanilla 03-CLK-SRSF-expr-correlations.R

# run CLK1-SRSF protein/phosphoprotein expression correlation script
Rscript --vanilla 04-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R

# plot CLK1 Ex 4 PSI
Rscript --vanilla 05-plot-Ex4-PSI-all-tumors.R

# plot CLK1 Ex 4 and all functionally hit splicing factor kinase PSIs
Rscript --vanilla 06-plot-kinase-PSI-variances-tumors.R

# plot CLK1 Ex 4 transcripts proportion in tumors and control samples
Rscript --vanilla 07-plot-clk1-ex4-transcripts-normals.R

# Plot CLK1 exon 4 PSI variations across tumors by age
Rscript --vanilla 08-plot-clk1-ex4-transcripts-by-age.R

# Compare high vs low Exon 4 PSI samples
Rscript --vanilla 09-volcano-plot-high-vs-low-Ex4-PSI.R
