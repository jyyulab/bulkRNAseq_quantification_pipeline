---
title: I-Software Installation
layout: default
nav_order: 1
parent: Pipeline Setup
---

### Part I: Software installation

We selcted the tools for this pipeline mainly based on two considerations: 1) they are well-established and widely-used; 2) they can work together with each other. We finally managed to complie all tools in one single conda environment. To install them:

1. Install **conda** if it's not available yet. This can be done by following this [tutorial](https://www.anaconda.com/docs/getting-started/getting-started).

2. Create a conda env for this pipeline:

   ``` shell
   conda create --prefix /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025 python=3.9 r-base=4.4
   ```

3. Install key software and dependencies:

   ``` shell
   ## activate the conda env
   conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025
   
   ## install tools from bioconda channel
   conda install -c bioconda rseqc=5.0.4 fastp=1.0.1 biobambam=2.0.185 samtools=1.22.1 fastqc=0.12.1 bedtools=2.31.1 bowtie2=2.5.4 rsem=1.3.3 star=2.7.11b salmon=1.10.3 cutadapt=5.1 htseq=2.0.9 ucsc-genepredtobed=482 ucsc-gtftogenepred=482
   
   ## install dependencies from conda-forge
   conda install -c conda-forge pandoc=3.7.0.2
   
   ## install R packages from r channel
   conda install -c r r-rmarkdown=2.29 r-ggplot2=3.5.2 r-dplyr=1.1.4 r-envstats=3.1.0 r-kableextra=1.4.0 r-rjson=0.2.23 r-cowplot=1.2.0 r-plotly=4.11.0
   
   ## deactivate conda env
   conda deactivate
   ```
