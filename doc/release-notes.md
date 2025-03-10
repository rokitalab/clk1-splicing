# release notes
## current releease (v11)
- Data release date: 2024-09-18
- Data release update (CNV consensus file): 2024-07-07
- OpenPedCan (OPC) data release date: 2024-01-03 (v13)
	- v14 pre-release `histologies.tsv` was created in [PR 542](https://github.com/d3b-center/OpenPedCan-analysis/pull/542) and is used here
- status: available

Modifications:
- `GSE243682-normal-brain-histologies.tsv`: Publicly available pediatric normal brain samples (GSE243682) metadata
- `GSE243682-normal-splice-events-rmats.tsv.gz`: rMATs results from publicly available pediatric normal brain samples (GSE243682)   
- `GSE73721-normal-histologies.tsv` Publicly available pediatric normal astrocyte samples (GSE73721) metadata
- `GSE73721-normal-rna-isoform-expression-rsem-tpm.rds` transcript expression results from publicly available normal astrocyte samples (GSE73721)
- `GSE73721-normal-splice-events-rmats.tsv.gz` rMATs results from publicly available pediatric normal brain samples (GSE73721)   
- `evodevo-histologies.tsv` Evo-Devo (developmental-normal) metadata
- `evodevo_rna-isoform-expression-rsem-tpm.rds` transcript expression results of Evo-Devo (developmental-normal) metadata
- `gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz`: rMATs results of Evo-Devo (developmental-normal)
- `gtex-harmonized-isoform-expression-rsem-tpm.rds`: transcript expression results of GTex samples

## current release (v10)
- Data release date: 2024-09-18
- Data release update (CNV consensus file): 2024-07-07
- OpenPedCan (OPC) data release date: 2024-01-03 (v13)
	- v14 pre-release `histologies.tsv` was created in [PR 542](https://github.com/d3b-center/OpenPedCan-analysis/pull/542) and is used here
- status: available

Modifications:
- `histologies-plot_group.tsv` was re-generated to remove samples corresponding to patients >= 40 years of age and updated PNOC with pathology dx of DIPG to have DIPG cancer_group
v10
.
├── CLK1-CRISPR-DepMap-score.csv
├── clk1-splice-events-rmats.tsv
├── control-rna-isoform-expression-rsem-counts-tpm.rds
├── OmicsDefaultModelProfiles.csv
├── OmicsExpressionTranscriptsTPMLogp1Profile.csv
├── consensus_wgs_plus_freec_wxs_plus_freec_tumor_only.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
├── ctrl-vs-morpholino-merged-rmats.tsv
├── ctrl_vs_morpho.rsem.genes.results.tsv
├── fusion-putative-oncogenic.tsv
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
├── histologies-plot-group.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── morpholno.merged.rmats.tsv
├── nf1-splice-events-rmats.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-expected-counts.rds
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
└── splice-events-rmats.tsv.gz

## current release (v9)
- Data release date: 2024-07-01
- Data release update (CNV consensus file): 2024-07-07
- OpenPedCan (OPC) data release date: 2024-01-03 (v13)
	- v14 pre-release `histologies.tsv` was created in [PR 542](https://github.com/d3b-center/OpenPedCan-analysis/pull/542) and is used here
- status: available

Modifications:
- `snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz` was re-generated (OPC pre-v16 release) and now contains SIFT/PolyPhen scores
- OPC v15 modified CNV consensus file `consensus_wgs_plus_freec_wxs_plus_freec_tumor_only.tsv.gz` was created in [PR 591](https://github.com/d3b-center/OpenPedCan-analysis/pull/591) and is used here

v9
.
├── CLK1-CRISPR-DepMap-score.csv
├── clk1-splice-events-rmats.tsv
├── control-rna-isoform-expression-rsem-counts-tpm.rds
├── OmicsDefaultModelProfiles.csv
├── OmicsExpressionTranscriptsTPMLogp1Profile.csv
├── consensus_wgs_plus_freec_wxs_plus_freec_tumor_only.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
├── ctrl-vs-morpholino-merged-rmats.tsv
├── ctrl_vs_morpho.rsem.genes.results.tsv
├── fusion-putative-oncogenic.tsv
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
├── histologies-plot-group.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── morpholno.merged.rmats.tsv
├── nf1-splice-events-rmats.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-expected-counts.rds
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
└── splice-events-rmats.tsv.gz

## previous release (v8)
- Data release date: 2024-04-05
- OpenPedCan data release date: 2024-01-03 (v13)
	- v14 pre-release `histologies.tsv` was created in [PR 542](https://github.com/d3b-center/OpenPedCan-analysis/pull/542) and is used here
- status: available

Additional files:
- `clk1-splice-events-rmats.tsv` and `nf1-splice-events-rmats.tsv` files from `CLK1-splicing_correlations` module results folder added
- `histologies-plot-group.tsv` from `cohort_summary` module added
- `snv-mutation-tmb-coding.tsv` and `snv-mutation-tmb-all.tsv` added from [PR 560](https://github.com/d3b-center/OpenPedCan-analysis/pull/560/commits/256ef04c3af141a8289a0efa68cf63b8aa4e30a4) to include missed samples
- `rna-isoform-expression-rsem-tpm.rds` and `rna-isoform-expression-rsem-expected-counts.rds` added, but subsetted for PBTA samples
- `gene-counts-rsem-expected_count-collapsed.rds` subsetted for PBTA samples
- `gene-expression-rsem-tpm-collapsed.rds` subsetted for PBTA samples
- `GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz` from GTEX portal
- `control-rna-isoform-expression-rsem-counts-tpm.rds` GENCODE v39 processed controls

v8
.
├── CLK1-CRISPR-DepMap-score.csv
├── clk1-splice-events-rmats.tsv
├── control-rna-isoform-expression-rsem-counts-tpm.rds
├── OmicsDefaultModelProfiles.csv
├── OmicsExpressionTranscriptsTPMLogp1Profile.csv
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
├── ctrl-vs-morpholino-merged-rmats.tsv
├── ctrl_vs_morpho.rsem.genes.results.tsv
├── fusion-putative-oncogenic.tsv
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
├── histologies-plot-group.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── morpholno.merged.rmats.tsv
├── nf1-splice-events-rmats.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-expected-counts.rds
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
└── splice-events-rmats.tsv.gz

## previous release (v7)
- Data release date: 2024-01-08
- OpenPedCan data release date: 2024-01-03 (v13)
	- v14 pre-release `histologies.tsv` was created in [PR 542](https://github.com/d3b-center/OpenPedCan-analysis/pull/542) and is used here
- status: available

Additional files:
- `splice-events-rmats-pbta.tsv.gz`: rmats subsetted by PBTA cohort
- `snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz`: maf for tumor only samples
- `cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz`: Phospho-Proteomics samples from CPTAC
- `cptac-protein-imputed-prot-expression-abundance.tsv.gz` Whole Cell Proteomics samples from CPTAC
- `gbm-protein-imputed-phospho-expression-abundance.tsv.gz` gbm protein-phsop expression tsv
- `gbm-protein-imputed-prot-expression-abundance.tsv.gz` gbm protein expression tsv
- `hope-protein-imputed-phospho-expression-abundance.tsv.gz` Phospho-Proteomics samples from HOPE
- `hope-protein-imputed-prot-expression-abundance.tsv.gz` Whole Cell Proteomics samples from HOPE
```
v7
.
├── CLK1-CRISPR-DepMap-score.csv
├── OmicsDefaultModelProfiles.csv
├── OmicsExpressionTranscriptsTPMLogp1Profile.csv
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
├── ctrl-vs-morpholino-merged-rmats.tsv
├── ctrl_vs_morpho.rsem.genes.results.tsv
├── fusion-putative-oncogenic.tsv
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── morpholno.merged.rmats.tsv
├── release-notes.md
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
└── splice-events-rmats-pbta.tsv.gz
```


## previous release (v6)
- Data release date: 2023-09-08
- OpenPedCan data release date: 2023-04-30 (v12)
- status: available

Additional files:
- `morpholno.merged.rmats.tsv` : merged matrix of rMATs splice events, comparing morphilino treated vs untreated KNS42 cells
- `ctrl_vs_morpho.rsem.genes.results.tsv`: count data matrix of morphilino treated and untreated KNS42 cells
- `CLK1_CRISPR_depmap_score.csv` : DepMap portal database file for dependency
- `OmicsDefaultModelProfiles.csv` : Profile mappings of cell lines and omics data

```
v6
├── cnv-consensus.seg.gz
├── consensus_wgs_plus_cnvkit_wxs.tsv.gz
├── fusion-putative-oncogenic.tsv
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── md5sum.txt
├── rMATS_merged.comparison.tsv.gz
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-coding.tsv
└── splice-events-rmats.tsv.gz
├── morpholno.merged.rmats.tsv
├── ctrl_vs_morpho.rsem.genes.results.tsv
├── CLK1_CRISPR_depmap_score.csv
└── OmicsDefaultModelProfiles.csv
```

## previous release (v5)
- Data release date: 2023-05-26
- OpenPedCan data release date: 2023-04-30 (v12)
- status: available

Additional files:
- `splice-events-rmats.tsv.gz`: merged matrix of rMATs splice events, run as single sample
- `rMATS_merged.comparison.tsv.gz`: Midline HGG sample splice events compared to non-diseased brainstem, run using rMATs
- `release-notes.md`: specific release notes to this repository

```
v5
├── cnv-consensus.seg.gz
├── consensus_wgs_plus_cnvkit_wxs.tsv.gz
├── fusion-putative-oncogenic.tsv
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── md5sum.txt
├── rMATS_merged.comparison.tsv.gz
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-coding.tsv
└── splice-events-rmats.tsv.gz
```

## previous release (v4)
- Data release date: 2023-04-30
- OpenPedCan data release date: 2023-04-30 (v12)
- status: available

Additional files:
- `splice-events-rmats.tsv.gz`: merged matrix of rMATs splice events, run as single sample

```
v4
├── cnv-consensus.seg.gz
├── consensus_wgs_plus_cnvkit_wxs.tsv.gz
├── fusion-putative-oncogenic.tsv
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── md5sum.txt
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-coding.tsv
└── splice-events-rmats.tsv.gz
```
