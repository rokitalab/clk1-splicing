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

## Directory structure
```
.
├── 00-create-psi-matrix.R
├── 01-sample_clustering.R
├── 02-plot-histology-distr-across-clusters.R
├── 03-plot-sbi-with-cluster-mem.R
├── 04-diff_pathways_per_cluster.R
├── README.md
├── input
│   └── subtype_hex.tsv
├── plots
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-all libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-all libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── cluster_membership-subtypes_poly-A-stranded.pdf
│   ├── cluster_membership-subtypes_stranded.pdf
│   ├── cluster_membership_poly-A-stranded.pdf
│   ├── cluster_membership_sbi_group_poly-A-stranded.pdf
│   ├── cluster_membership_sbi_group_stranded.pdf
│   ├── cluster_membership_stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-all libraries.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-poly-A stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-all libraries.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-poly-A stranded.pdf
│   ├── hgg-dmg-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-all libraries.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-poly-A stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-all libraries.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-poly-A stranded.pdf
│   ├── lgg-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-all libraries.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-poly-A stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-all libraries.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-poly-A stranded.pdf
│   ├── mb-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-all libraries.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-poly-A stranded.pdf
│   ├── sample-cluster-histology-enr-top-1000-events-stranded.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-all libraries.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-poly-A stranded.pdf
│   ├── sample-cluster-histology-enr-top-5000-events-stranded.pdf
│   ├── sample-psi-heatmap-top-1000-events-all libraries.pdf
│   ├── sample-psi-heatmap-top-1000-events-poly-A stranded.pdf
│   ├── sample-psi-heatmap-top-1000-events-stranded.pdf
│   ├── sample-psi-heatmap-top-5000-events-all libraries.pdf
│   ├── sample-psi-heatmap-top-5000-events-poly-A stranded.pdf
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
│   ├── pbta-splice-event-psis.RDS
│   ├── psi-matrix-top-1000-events-all libraries.rds
│   ├── psi-matrix-top-1000-events-poly-A stranded.rds
│   ├── psi-matrix-top-1000-events-stranded.rds
│   ├── psi-matrix-top-5000-events-all libraries.rds
│   ├── psi-matrix-top-5000-events-poly-A stranded.rds
│   ├── psi-matrix-top-5000-events-stranded.rds
│   ├── sample-cluster-metadata-top-1000-events-all libraries.tsv
│   ├── sample-cluster-metadata-top-1000-events-poly-A stranded.tsv
│   ├── sample-cluster-metadata-top-1000-events-stranded.tsv
│   ├── sample-cluster-metadata-top-5000-events-all libraries.tsv
│   ├── sample-cluster-metadata-top-5000-events-poly-A stranded.tsv
│   └── sample-cluster-metadata-top-5000-events-stranded.tsv
├── run-module.sh
└── util
    └── heatmap_function.R
```
