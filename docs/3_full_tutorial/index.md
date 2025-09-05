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

1. Please make sure you have set up this pipeline in your end.

   This pipeline, including its **tools**, **databases**, and **scripts**, can be set up and maintained in a single conda environment. We have set up this pipeline in the conda environment below, with four pre-built genome assemblies: hg38, hg19, mm39 and mm10. You can activate it using the following commands:

   ``` bash
   module load conda3/202402 # conda version: 24.1.2
   conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025 # change it to your conda environment accordingly
   ```

   - If you can access to it but require a different reference genome assembly, please refer to the [Database Preparation](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/1_pipeline_setup/2_database.html) tutorial.
   - To set up a conda environment for this pipeline in you end, please refer to the [Pipeline Setup](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/pipeline_setup/index) tutorial. It typically takes **~3** hours to complete.

2. Please make sure you have the input data in your hands.

   This pipeline accepts both **FASTQ files** and pre-computated alignment (**BAM/SAM**) files. If you do not have input samples to start with, you can use the test data we provide in this pipeline:

   | Sample ID | Species | Library Type | File Format           | # Raw Reads               | Preparation                |
   | --------- | ------- | ------------ | --------------------- | ------------------------- | -------------------------- |
   | sample1   | Human   | PE-100       | FASTQ                 | 38,820,746 (19,410,373*2) | Downsampled from real data |
   | sample2   | Mouse   | SE-50        | FASTQ                 | 16,005,450                | Downsampled from real data |
   | sample3   | Human   | PE-100       | FASTQ, split by lanes | 38,820,746 (19,410,373*2) | Split from sample1         |
   | sample4   | Mouse   | SE-50        | FASTQ, split by lanes | 16,005,450                | Split from sample2         |
   | sample5   | Human   | PE-100       | BAM                   | 27,712,150 (13,856,075*2) | Downsampled from real data |
   | sample6   | Mouse   | SE-50        | BAM                   | 13,533,152                | Downsampled from real data |
   
   
