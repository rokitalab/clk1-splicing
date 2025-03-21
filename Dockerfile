FROM rocker/tidyverse:4.4.0
LABEL maintainer="Ammar S. Naqvi (naqvia@chop.edu)"
WORKDIR /rocker-build/


### Install apt-getable packages to start
#########################################
RUN apt-get -y update && apt-get install -y bedtools \
	build-essential \
	bzip2 \
	cpanminus \
	cmake \
	curl \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgdal-dev \
        libglpk40 \
        libglpk-dev \
        libhdf5-dev \
	libgmp-dev \
	liblzma-dev \
	libmagick++-dev \
	libmpfr-dev \
	libncurses5-dev \
	libproj-dev \
	libreadline-dev \
	libssl-dev \
	libv8-dev \
	libxt-dev \
	libudunits2-dev \
	zlib1g-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6

RUN R -e "remotes::install_github('YuLab-SMU/yulab.utils',dependencies = TRUE)"

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
#RUN R -e "BiocManager::install(version = '3.19')"
RUN R -e "options(repos = c(CRAN = 'https://cran.rstudio.com/')); BiocManager::install(version = '3.19', ask = FALSE)"

# Install packages
RUN R -e 'BiocManager::install(c( \
  "AnnotationDbi", \
  "Biobase", \
  "Biostrings", \
  "broom", \
  "circlize", \
  "COINr", \
  "coin", \
  "ComplexHeatmap", \
  "ConsensusClusterPlus", \
  "corrplot", \
  "cowplot", \
  "DGCA", \
  "DESeq2", \
  "DOSE", \
  "diptest", \
  "edgeR", \
  "EnhancedVolcano", \
  "factoextra", \
  "fgsea", \
  "fpc", \
  "ggpubr", \
  "ggstatsplot", \
  "ggthemes", \
  "ggVennDiagram", \
  "gridExtra", \
  "GO.db", \
  "GSVA", \
  "Hmisc", \
  "hrbrthemes", \
  "impute", \
  "limma", \
  "lspline", \
  "maftools",\
  "msigdbr", \
  "NOISeq", \
  "optparse", \
  "org.Hs.eg.db", \
  "PMCMRplus", \
  "pheatmap", \
  "pwalign", \
  "preprocessCore", \
  "reshape2", \
  "rstatix", \
  "rtracklayer", \
  "R.utils", \
  "sva", \
  "survival", \
  "survminer", \
  "sva", \
  "WGCNA", \
  "VennDiagram",\
  "UpSetR" \
))'


## install GitHub packages
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('d3b-center/annoFuseData', ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66', dependencies = TRUE)"
#RUN R -e "remotes::install_github('rcastelo/GSVA', ref = 'df9001cfd07017001dfba07a3099e6b7dc5ce324', dependencies = TRUE)"
RUN R -e "remotes::install_github('andymckenzie/DGCA', ref = '075fc79a32df3955e75e262b8269c257d8ffac9c', dependencies = TRUE)"
RUN R -e "remotes::install_github('YuLab-SMU/yulab.utils',dependencies = TRUE)"

# install clusterProfiler
RUN R -e 'BiocManager::install("clusterProfiler")'


# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
