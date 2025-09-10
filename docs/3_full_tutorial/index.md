---
title: Full Tutorial
layout: default
nav_order: 4
permalink: /docs/3_full_tutorial/index
---

# Full Tutorial for Runing the Pipeline

---

Welcome to the **Full Tutorial** for the bulk RNA-seq quantification pipeline! This tutorial provides **comprehensive documentation** and **step-by-step** guidance for this pipeline. and it is intended for users who are new to this field. From this tutorial, you will find detailed instructions to the **rationale**, **inputs**, **outputs** and **key parameters** for each step. If you want to quickly run this pipeline on your own data, please refer to the [**Quick Tutorial**](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/2_quick_tutorial/quick_tutorial).



### Before you get started

1. **Ensure this pipeline is set up on your system.**

   This pipeline, including its **tools**, **scripts**, and **databases**, can be set up and maintained in a single conda environment. We have set up this pipeline in the conda environment below, with four pre-built genome assemblies: **hg38**, **hg19**, **mm39** and **mm10**. You can activate it using the following commands:

   ``` bash
   module load conda3/202402 # conda version: 24.1.2
   conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025 # change it to your conda environment accordingly
   ```

   - If you can access to this environment but require a different reference genome assembly, please refer to the [Database Preparation](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/1_pipeline_setup/2_database.html) tutorial.
   - To set up your own conda environment for this pipeline, please refer to the [Pipeline Setup](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/pipeline_setup/index) tutorial. It typically takes **~3** hours to complete.

2. **Ensure you have input data available**

   This pipeline accepts both **FASTQ files** and pre-computated alignment (**BAM/SAM**) files. We always encourage you to start walking through this pipeline with your own data. However, if you do not have your own input data to start with, you can use the test data provided with this pipeline (see the table below):

   - For those who can access to the pre-built conda environment, the test data is available at: `/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata`. 
   - For other users, please download them from [Zenodo](https://zenodo.org/records/17088549).
   
   |         | Format                     | Library Type | Phred Encoding | Species | Number of Reads | Preparation                |
   | ------- | -------------------------- | ------------ | -------------- | ------- | --------------- | -------------------------- |
   | sample1 | FASTQ                      | Paired-end   | Phred+33       | human   | 19,410,373 * 2  | Downsampled from real data |
   | sample2 | FASTQ                      | Single-end   | Phred+33       | mouse   | 16,005,450      | Downsampled from real data |
   | sample3 | FASTQ, multiple lanes      | Paired-end   | Phred+33       | human   | 19,410,373 * 2  | Split by lane from sample1 |
   | sample4 | FASTQ, multiple lanes      | Single-end   | Phred+33       | mouse   | 16,005,450      | Split by lane from sample2 |
   | sample5 | BAM/SAM (to genome)        | Paired-end   | Phred+33       | human   | 13,856,075 * 2  | Downsampled from real data |
   | sample6 | BAM/SAM (to transcriptome) | Single-end   | Phred+33       | mouse   | 13,533,162      | Downsampled from real data |
