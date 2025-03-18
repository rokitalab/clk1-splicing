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
├── plots
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-all libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-poly-A stranded.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-1000-events-stranded.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-all libraries.pdf
│   ├── atrt-sample-cluster-subtype-enr-top-5000-events-stranded.pdf
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
│   └── sample-psi-heatmap-top-5000-events-stranded.pdf
├── results
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
