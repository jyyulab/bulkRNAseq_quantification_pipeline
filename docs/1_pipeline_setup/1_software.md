---
title: I-Software Installation
layout: default
nav_order: 1
parent: Pipeline Setup
---

## Part I: Software Installation

As mentioned before, we have manged to install and maintain all required tools and dependencies in one single conda environment. To set up this computational environment:

1. Make sure [**conda**](https://docs.conda.io/projects/conda/en/latest/index.html) is already available in your case. If not, please install it following this [guidance](https://www.anaconda.com/docs/getting-started/getting-started).

   ``` bash
   module load conda/202402 # version 24.1.2
   ```

   NOTE: Though **we always encourage the users to create the conda environment and install the required tools and dependencies manually following #2 and #3 blow**, we do provide the **`environment.yml`** file (available [here](https://github.com/jyyulab/bulkRNAseq_quantification_pipeline/blob/main/scripts/environment.yml)) so that the users can directly check the tools and dependencies and easily create the conda environment using: 

   ``` bash
   conda env create -f environment.yml -n bulkRNAseq_2025 -p /your-own-path/bulkRNAseq_2025
   ```

2. Create a conda environment for this pipeline:

   ``` shell
   conda create --name bulkRNAseq_2025 --prefix /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025 python=3.9 r-base=4.4
   ```

   - Please modify the `--prefix` accordingly. If ignore it, the conda enviroment will be installed in the default directory.
   - Make sure `python=3.9` and `r-base=4.4` are specified. There might be some compatiable issues for the default python and R version.

3. Install key software and dependencies:

   ``` shell
   ## activate the conda env
   conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025 # please modify accordingly
   
   ## install tools from bioconda channel
   conda install -c bioconda rseqc=5.0.4 fastp=1.0.1 biobambam=2.0.185 samtools=1.22.1 fastqc=0.12.1 bedtools=2.31.1 bowtie2=2.5.4 rsem=1.3.3 star=2.7.11b salmon=1.10.3 cutadapt=5.1 htseq=2.0.9 ucsc-genepredtobed=482 ucsc-gtftogenepred=482
   
   ## install dependencies from conda-forge
   conda install -c conda-forge pandoc=3.7.0.2
   
   ## install R packages from r channel
   conda install -c r r-rmarkdown=2.29 r-ggplot2=3.5.2 r-dplyr=1.1.4 r-envstats=3.1.0 r-kableextra=1.4.0 r-rjson=0.2.23 r-cowplot=1.2.0 r-plotly=4.11.0
   
   ## deactivate conda env
   conda deactivate
   ```

4. Deposite the scripts

   Thearetically, the scripts of this pipeline can be saved in any path that you can access. However, we HIGHLY recommend the users to deposite them in this path: `/your-conda-env-path/pipeline/scripts`. In the showcase, it's `/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/scripts`.

   These scripts, by their applications, can be divided into three groups:

   - `/scripts/run/*.pl`: These scripts are used to run the pipeline step-by-step.
   - `/scripts/setup/*.pl`: These scripts are used in pipeline setup.
   - `/scripts/qc/*.Rmd`: These scripts are used to generate the quality control reports.

   ``` bash
   ## All scripts are available here: https://github.com/jyyulab/bulkRNAseq_quantification_pipeline/tree/main/scripts.
   
   ## To download them:
   ## 1) visit: https://download-directory.github.io/
   ## 2) paste this link: https://github.com/jyyulab/bulkRNAseq_quantification_pipeline/tree/main/scripts
   ## 3) press Enter
   ## 4) unzip the downloaded file and save them to your folder
   ```

   
