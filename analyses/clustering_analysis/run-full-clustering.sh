# 0) scripts to create matrix and input
# pre-process and create neccessary input for clustering
perl code/00-create_matrix_of_PSI_SE_gene.pl ../../data/histologies-plot-group.tsv ../../data/splice-events-rmats.tsv.gz ../../data/independent-specimens.rnaseqpanel.primary.tsv input/pan_cancer_splicing_SE.gene.tsv
Rscript --vanilla code/01-convert_to_rds.R
Rscript --vanilla code/02-create-inputs.R

# 1) run script for optimal clustering
Rscript code/03-optimal-clustering.R \
--input_mat "input/pan_cancer_splicing_SE.gene.rds" \
--cluster_algorithm "hc, km, pam" \
--cluster_distance "pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski" \
--filter_expr FALSE \
--protein_coding_only FALSE \
--feature_selection "dip.test" \
--transformation_type "none" \
--max_k 15

# 1) PBTA splicing data
# get ccp clustering output for a specific combination of distance + algorithm + % variable genes
# for the splicing dataset, pam + binary k = 12, 0% based on optimal clustering
Rscript code/04-get-clustering-output.R \
--input_mat "input/pan_cancer_splicing_SE.gene.rds" \
--data_type "non_expr" \
--var_genes "0" \
--cluster_algorithm "pam" \
--cluster_distance "binary" \
--prefix "non_expr_pan_cancer_splice_subset"

# get differential genes per cluster and perform pre-ranked gsea using those genes
# this was done for k = 12 with pam + binary + 0% genes
Rscript code/05-diff-genes-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_binary_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_binary_0_ccp.rds" \
--n_cluster "binary" \
--gene_set "input/hallmark_splice_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_pam_binary_0" \
--output_dir "output/diff_genes"

# get differential pathways per cluster using GSVA
# this was done for k = 2 with pam  + binary + 0% genes
Rscript code/06-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_binary_0_matrix.rds" \
--input_clin "../../data/histologies-plot-group.tsv" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_binary_0_ccp.rds" \
--n_cluster "12" \
--gene_set "input/hallmark_splice_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_pam_binary_0" \
--output_dir "output/diff_pathways"

# get heatmap of the CCP matrix
Rscript code/07-plot-clustering-heatmap.R \
--ccp_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_binary_0_ccp.rds" \
--input_clin "../../data/histologies-plot-group.tsv" \
--n_cluster "12" \
--prefix "non_expr_pan_cancer_splice_subset"

##plot cluster members categorized by SBI high vs low
Rscript code/08-plot-sbi_with_cluster-mem.R

## plot histologies across clusters
Rscript code/09-plot-histology-distr-across-clusters.R
