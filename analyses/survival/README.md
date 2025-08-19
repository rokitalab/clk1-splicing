# Survival by splicing burden

Module authors: Ryan Corbett (@rjcorb), Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess survival in CNS tumor pediatric patients by cluster, overall splicing burden index, spliceosome GVSA score, and CLK1 TPM and exon 4 PSI.

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run-survival-module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
../splicing_index/results/splicing_index.SE.txt
```

## Folder content
* `run-survival-module.sh` shell script to run analysis
* `01-prepare-survival.Rmd` Merge splicing data with relevant histology data and define survival features
* `02-survival_by_cluster.Rmd` Assess survival by splicing cluster assignment in all histologies, LGG, and HGG
* `03-survival_by_cluster_sbi_spliceosome_gsva.Rmd` Plot survival by cluster and SI or GSVA status
* `04-survival_by_cluster_clk1.Rmd` Plot survival by cluster and CLK1 TPM and PSI

## Directory structure
```
.
├── 01-prepare-survival.nb.html
├── 01-prepare-survival.Rmd
├── 02-survival_by_cluster.nb.html
├── 02-survival_by_cluster.Rmd
├── 03-survival_by_cluster_sbi_spliceosome_gsva.nb.html
├── 03-survival_by_cluster_sbi_spliceosome_gsva.Rmd
├── 04-survival_by_cluster_clk1.nb.html
├── 04-survival_by_cluster_clk1.Rmd
├── plots
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score_CLK1_Ex4_TPM.pdf
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score_SBI_SE.pdf
│   ├── forest_add_EFS_cluster6_histology_resection_clk1_psi_group.pdf
│   ├── forest_add_EFS_cluster6_histology_resection_clk1_psi.pdf
│   ├── forest_add_EFS_cluster6_histology_resection_clk1_tpm_group.pdf
│   ├── forest_add_EFS_cluster6_histology_resection_clk1_tpm.pdf
│   ├── forest_add_EFS_cluster6_histology_resection_spliceosome_group.pdf
│   ├── forest_add_EFS_HGG_subtype_cluster_assignment_SBI.pdf
│   ├── forest_add_EFS_HGG_subtype_cluster_assignment_spliceosome_score.pdf
│   ├── forest_add_EFS_HGG_subtype_cluster_assignment.pdf
│   ├── forest_add_EFS_LGG_resection_subtype_cluster_assignment_SBI.pdf
│   ├── forest_add_EFS_LGG_resection_subtype_cluster_assignment_spliceosome_score.pdf
│   ├── forest_add_EFS_LGG_resection_subtype_cluster_assignment.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_CLK1_ex4_PSI.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_CLK1_Ex4_TPM.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_SBI_CLK1_Ex4_TPM.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_SBI.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_spliceosome_score_CLK1_Ex4_TPM.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment_spliceosome_score.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment.pdf
│   ├── forest_add_OS_cluster6_histology_resection_clk1_psi_group.pdf
│   ├── forest_add_OS_cluster6_histology_resection_clk1_psi.pdf
│   ├── forest_add_OS_cluster6_histology_resection_clk1_tpm_group.pdf
│   ├── forest_add_OS_cluster6_histology_resection_clk1_tpm.pdf
│   ├── forest_add_OS_cluster6_histology_resection_si.pdf
│   ├── forest_add_OS_HGG_subtype_cluster_assignment_SBI.pdf
│   ├── forest_add_OS_HGG_subtype_cluster_assignment_spliceosome_score.pdf
│   ├── forest_add_OS_HGG_subtype_cluster_assignment.pdf
│   ├── forest_add_OS_LGG_resection_subtype_cluster_assignment.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_CLK1_ex4_PSI.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_CLK1_Ex4_TPM.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_SBI_CLK1_Ex4_TPM.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_SBI.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_spliceosome_score_CLK1_Ex4_TPM.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment_spliceosome_score.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment.pdf
│   ├── forest_int_EFS_resection_lgg_group_cluster_assignment_CLK1_ex4_PSI.pdf
│   ├── forest_int_EFS_resection_lgg_group_cluster_assignment_CLK1_Ex4_TPM.pdf
│   ├── forest_int_EFS_resection_lgg_group_cluster_assignment_SBI_SE.pdf
│   ├── forest_int_EFS_resection_lgg_group_cluster_assignment_spliceosome_gsva_score.pdf
│   ├── forest_int_EFS_resection_lgg_group_cluster_clk1_age.pdf
│   ├── forest_int_OS_resection_lgg_group_cluster_assignment_CLK1_ex4_PSI.pdf
│   ├── forest_int_OS_resection_lgg_group_cluster_assignment_CLK1_Ex4_TPM.pdf
│   ├── forest_int_OS_resection_lgg_group_cluster_assignment_SBI_SE.pdf
│   ├── forest_int_OS_resection_lgg_group_cluster_assignment_spliceosome_gsva_score.pdf
│   ├── forest_int_OS_resection_lgg_group_cluster_clk1_age.pdf
│   ├── km_cluster6_EFS_clk1_psi_group.pdf
│   ├── km_cluster6_EFS_clk1_tpm_group.pdf
│   ├── km_cluster6_EFS_sbi_group.pdf
│   ├── km_cluster6_EFS_splice_group.pdf
│   ├── km_cluster6_OS_clk1_psi_group.pdf
│   ├── km_cluster6_OS_clk1_tpm_group.pdf
│   ├── km_cluster6_OS_sbi_group.pdf
│   ├── km_cluster6_OS_splice_group.pdf
│   ├── km_EFS_cluster_assignment.pdf
│   ├── km_hgg_EFS_cluster_assignment.pdf
│   ├── km_hgg_EFS_spliceosome_score.pdf
│   ├── km_hgg_OS_cluster_assignment.pdf
│   ├── km_hgg_OS_spliceosome_score.pdf
│   ├── km_lgg_EFS_cluster_assignment.pdf
│   ├── km_lgg_OS_cluster_assignment.pdf
│   └── km_OS_cluster_assignment.pdf
├── README.md
├── results
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_CLK1_ex4_PSI.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_CLK1_Ex4_TPM.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_SBI_CLK1_Ex4_TPM.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_SBI.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_spliceosome_score_CLK1_Ex4_TPM.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster_spliceosome_score.RDS
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster.RDS
│   ├── cox_EFS_additive_terms_subtype_cluster_clk1_psi_group.RDS
│   ├── cox_EFS_additive_terms_subtype_cluster_clk1_psi.RDS
│   ├── cox_EFS_additive_terms_subtype_cluster_clk1_tpm_group.RDS
│   ├── cox_EFS_additive_terms_subtype_cluster_clk1_tpm.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_CLK1_ex4_PSI.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_CLK1_Ex4_TPM_age.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_CLK1_Ex4_TPM.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_SBI_SE.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score_CLK1_Ex4_TPM.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score_SBI_SE.RDS
│   ├── cox_EFS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score.RDS
│   ├── cox_hgg_EFS_additive_terms_subtype_cluster_SBI.RDS
│   ├── cox_hgg_EFS_additive_terms_subtype_cluster_spliceosome_score.RDS
│   ├── cox_hgg_EFS_additive_terms_subtype_cluster.RDS
│   ├── cox_hgg_OS_additive_terms_subtype_cluster_SBI.RDS
│   ├── cox_hgg_OS_additive_terms_subtype_cluster_si_group.RDS
│   ├── cox_hgg_OS_additive_terms_subtype_cluster_spliceosome_score.RDS
│   ├── cox_hgg_OS_additive_terms_subtype_cluster.RDS
│   ├── cox_lgg_EFS_additive_terms_resection_subtype_cluster_SBI.RDS
│   ├── cox_lgg_EFS_additive_terms_resection_subtype_cluster_spliceosome_score.RDS
│   ├── cox_lgg_EFS_additive_terms_resection_subtype_cluster.RDS
│   ├── cox_lgg_OS_additive_terms_resection_subtype_cluster.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_CLK1_ex4_PSI.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_CLK1_Ex4_TPM.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_SBI_CLK1_Ex4_TPM.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_SBI.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_spliceosome_score_CLK1_Ex4_TPM.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster_spliceosome_score.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster.RDS
│   ├── cox_OS_additive_terms_subtype_cluster_clk1_psi_group.RDS
│   ├── cox_OS_additive_terms_subtype_cluster_clk1_psi.RDS
│   ├── cox_OS_additive_terms_subtype_cluster_clk1_tpm_group.RDS
│   ├── cox_OS_additive_terms_subtype_cluster_clk1_tpm.RDS
│   ├── cox_OS_interaction_terms_resection_lgg_group_cluster_clk1_age.RDS
│   ├── cox_OS_interaction_terms_resection_lgg_group_cluster_CLK1_ex4_PSI.RDS
│   ├── cox_OS_interaction_terms_resection_lgg_group_cluster_CLK1_Ex4_TPM.RDS
│   ├── cox_OS_interaction_terms_resection_lgg_group_cluster_SBI_SE.RDS
│   ├── cox_OS_interaction_terms_resection_lgg_group_cluster_spliceosome_gsva_score.RDS
│   ├── logrank_cluster6_EFS_clk1_psi_group.RDS
│   ├── logrank_cluster6_EFS_clk1_tpm_group.RDS
│   ├── logrank_cluster6_EFS_SBI.RDS
│   ├── logrank_cluster6_EFS_splice_group.RDS
│   ├── logrank_cluster6_OS_clk1_psi_group.RDS
│   ├── logrank_cluster6_OS_clk1_tpm_group.RDS
│   ├── logrank_cluster6_OS_SBI.RDS
│   ├── logrank_cluster6_OS_splice_group.RDS
│   ├── logrank_hgg_EFS_cluster_assignment.RDS
│   ├── logrank_hgg_EFS_splice_group.RDS
│   ├── logrank_hgg_OS_cluster_assignment.RDS
│   ├── logrank_hgg_OS_splice_group.RDS
│   ├── logrank_lgg_EFS_cluster_assignment.RDS
│   ├── logrank_lgg_OS_cluster_assignment.RDS
│   ├── logrank_OS_cluster_assignment.RDS
│   ├── logrankEFS_cluster_assignment.RDS
│   ├── splicing_indices_with_survival.tsv
│   └── subtypes-for-survival.tsv
├── run-survival-module.sh
└── util
    └── survival_models.R
```
