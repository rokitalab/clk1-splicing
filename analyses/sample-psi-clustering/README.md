# Sample clustering by splice event PSI

Module authors: Ryan Corbett (@rjcorb)

The purpose of this module is cluster primary tumors by splice events with the most variable PSIs

## Usage
### script to run analysis
```
bash run-module.sh
```


## Folder content
* 00-create-psi-matrix.R; create PSI matrix across primary tumors
* 01-sample_clustering.R; cluster samples by most variable PSIs and by RNA-seq library type
* 02-plot-histology-distr-across-clusters.R; assess hist and subtype distribution across clusters
* 03-plot-sbi-with-cluster-mem.R; assess SBI group distribution across clusters
* 04-diff_pathways_per_cluster.R; assess pathway enrichment across clusters
* 05-generate-spliceosome-summary-file.sh; collate spliceosome GSVA data

# Input files
* lgg-braf-fusion-breakpoint-annotation.tsv; LGG BRAF fusion-positive BS IDs annotated with KIAA1549::BRAF fusion breakpoint group. Pulled from [pbta-germline-somatic repo](https://github.com/diskin-lab-chop/pbta-germline-somatic/blob/main/analyses/survival/input/lgg-braf-fusion-breakpoint-annotation.tsv)
* mb_shh_molecular_subtypes.tsv; molecularly-defined SHH-MB subgroups. Pulled from [OpenPedCan repo](https://github.com/rokitalab/OpenPedCan-Project-CNH/blob/dev/analyses/molecular-subtyping-MB/results/MB_molecular_subtype.tsv)

## Directory structure
```
.
├── 00-create-psi-matrix.R
├── 01-sample_clustering.R
├── 02-plot-histology-distr-across-clusters.R
├── 03-plot-sbi-with-cluster-mem.R
├── 04-diff_pathways_per_cluster.R
├── 05-generate-spliceosome-summary-file.sh
├── README.md
├── input
│   ├── lgg-braf-fusion-breakpoint-annotation.tsv
│   ├── mb_shh_molecular_subtypes.tsv
│   └── subtype_hex.tsv
├── plots
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-all_libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-poly-A_stranded.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-all_libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── cluster_membership-subtypes_poly-A-stranded.pdf
│   ├── cluster_membership-subtypes_stranded.pdf
│   ├── cluster_membership_poly-A-stranded.pdf
│   ├── cluster_membership_sbi_group_poly-A-stranded.pdf
│   ├── cluster_membership_sbi_group_stranded.pdf
│   ├── cluster_membership_stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-all_libraries.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-poly-A_stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-all_libraries.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-poly-A_stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── lgg-braf-fusion-breakpoint-cluster-enr-top-5000-events-poly-A-stranded.pdf
│   ├── lgg-braf-fusion-breakpoint-cluster-enr-top-5000-events-stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-all_libraries.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-poly-A_stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-all_libraries.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-poly-A_stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── mb-group34-subgroup-cluster-enr-top-5000-events-poly-A-stranded.pdf
│   ├── mb-group34-subgroup-cluster-enr-top-5000-events-stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-all_libraries.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-poly-A_stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-all_libraries.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-poly-A_stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── mb-shh-subgroup-cluster-enr-top-5000-events-stranded.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-all_libraries.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-poly-A_stranded.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-stranded.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-all_libraries.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-poly-A_stranded.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-stranded.pdf
│   ├── sample-psi-heatmap-top-1000-events-all_libraries.pdf
│   ├── sample-psi-heatmap-top-1000-events-poly-A_stranded.pdf
│   ├── sample-psi-heatmap-top-1000-events-stranded.pdf
│   ├── sample-psi-heatmap-top-5000-events-all_libraries.pdf
│   ├── sample-psi-heatmap-top-5000-events-poly-A_stranded.pdf
│   ├── sample-psi-heatmap-top-5000-events-stranded.pdf
│   ├── top5_pathways_poly-A-stranded.pdf
│   └── top5_pathways_stranded.pdf
├── results
│   ├── cluster_10_pathway_stranded.tsv
│   ├── cluster_11_pathway_stranded.tsv
│   ├── cluster_1_pathway_poly-A-stranded.tsv
│   ├── cluster_1_pathway_stranded.tsv
│   ├── cluster_2_pathway_poly-A-stranded.tsv
│   ├── cluster_2_pathway_stranded.tsv
│   ├── cluster_3_pathway_poly-A-stranded.tsv
│   ├── cluster_3_pathway_stranded.tsv
│   ├── cluster_4_pathway_poly-A-stranded.tsv
│   ├── cluster_4_pathway_stranded.tsv
│   ├── cluster_5_pathway_poly-A-stranded.tsv
│   ├── cluster_5_pathway_stranded.tsv
│   ├── cluster_6_pathway_poly-A-stranded.tsv
│   ├── cluster_6_pathway_stranded.tsv
│   ├── cluster_7_pathway_poly-A-stranded.tsv
│   ├── cluster_7_pathway_stranded.tsv
│   ├── cluster_8_pathway_poly-A-stranded.tsv
│   ├── cluster_8_pathway_stranded.tsv
│   ├── cluster_9_pathway_poly-A-stranded.tsv
│   ├── cluster_9_pathway_stranded.tsv
│   ├── gsva_output_poly-A-stranded.tsv
│   ├── gsva_output_stranded.tsv
│   ├── kegg_spliceosome_summary.tsv
│   ├── lgg-braf-fusion-cluster-membership-poly-A-stranded.tsv
│   ├── lgg-braf-fusion-cluster-membership-stranded.tsv
│   ├── mb-subgroup-cluster-membership-poly-A-stranded.tsv
│   ├── mb-subgroup-cluster-membership-stranded.tsv
│   ├── psi-matrix-top-1000-events-all_libraries.rds
│   ├── psi-matrix-top-1000-events-poly-A_stranded.rds
│   ├── psi-matrix-top-1000-events-stranded.rds
│   ├── psi-matrix-top-5000-events-all_libraries.rds
│   ├── psi-matrix-top-5000-events-poly-A_stranded.rds
│   ├── psi-matrix-top-5000-events-stranded.rds
│   ├── sample-cluster-metadata-top-1000-events-all_libraries.tsv
│   ├── sample-cluster-metadata-top-1000-events-poly-A_stranded.tsv
│   ├── sample-cluster-metadata-top-1000-events-stranded.tsv
│   ├── sample-cluster-metadata-top-5000-events-all_libraries.tsv
│   ├── sample-cluster-metadata-top-5000-events-poly-A_stranded.tsv
│   └── sample-cluster-metadata-top-5000-events-stranded.tsv
├── run-module.sh
└── util
    └── heatmap_function.R
```
