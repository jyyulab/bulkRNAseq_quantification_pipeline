# Bulk RNA-Seq Quantification Pipeline 2025

## Overview

![Picture](docs/figures/overview.png)

As shown above, this pipeline contains three stages:

#### 1. Data preprocessing ####

In this stage, the pipeline intakes raw inputs of variable formats (e.g., FASTQ, BAM/SAM or FASTA), and generates the **standard-in-format**, **clean-in-sequence** FASTQ files that can be directly used for quantification analysis.

#### 2. Quantification

In this stage, the pipeline generates the quantification meansurements at both gene- and transcript-levels. It supports three well-established and widely-used quantifiers:

- [**<u>Salmon</u>**](https://salmon.readthedocs.io/en/latest/salmon.html): an **alignment-free quantifier** with **wicked-fast speed** and **comarable accuracy**.

- [**<u>RSEM</u>**](https://github.com/bli25/RSEM_tutorial): an **alignment-based quantifier** with **high accuracy**. It has been used as **gold standard** in many benchmarking studies.

- [**<u>STAR</u>**](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf): an **alignment-based quantifier** featured by **spliced transcripts alignment. This is the tool used by [GDC mRNA quantification analysis pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline).

#### 3. Summary

In this stage, the pipeline generates **an HTML report** of quantification analyis for each sample, including alignment statistics, correlation analysis, gene body coverage visualizations, and more. When multiple samples are provided, it can also produce a universal report summarizing statistics of all samples, as well as  a master gene expression matrix that can be directly used for **NetBID** analysis.

## Tutorial

A detailed tutorial to ***set up*** and ***run*** this pipeline can be found here: https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/.
