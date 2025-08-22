# Bulk RNA-Seq Quantification Pipeline 2025

## Overview

![Picture](/Users/qpan/Library/CloudStorage/OneDrive-St.JudeChildren\'sResearchHospital/Typora/bulkRNAseq_pipeline/overview.png)

As shown above, this pipeline contains three parts:

1. Data preprocessing: it generates the standard-in-format, clean-in-sequence FASTQ files that can be directly proceed to quantification analysis.
2. Quantification: it supports three well-established and widely-used quantifiers:
  1) Salmon: an **alignment-free quantifier** with **wicked-fast speed** and **comarable accuracy** (https://salmon.readthedocs.io/en/latest/salmon.html).
  2) RSEM: an **alignment-based quantifier** with **high accuracy**. It has been used as **gold standard** in many benchmarking studies (https://github.com/bli25/RSEM_tutorial).
  3) STAR: supports spliced transcript alignment (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). It is recommended by GDC (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline).
3. Summary: it generates a html report of quantification analyis, including alignment statistics, correlation analysis, gene body coverage visualization, etc.

## Pipeline setup

We recommend to set up the pipeline in Conda environment.

1. Install conda if it's not available yet. https://www.anaconda.com/

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
   conda install -c r r-rmarkdown=2.29 r-ggplot2=3.5.2 r-dplyr=1.1.4 r-envstats=3.1.0 r-kableextra=1.4.0 r-rjson=0.2.23 r-cowplot=1.2.0
   
   ## deactivate conda env
   conda deactivate
   ```



## Tutorial

A detailed tutorial to set up and run this pipeline can be found here: https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/.
