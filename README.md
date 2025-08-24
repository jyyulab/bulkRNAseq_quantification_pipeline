# Bulk RNA-Seq Quantification Pipeline 2025

## Overview

![Picture](docs/figures/overview.png)

As shown above, this pipeline contains three parts:

1. Data preprocessing: it generates the standard-in-format, clean-in-sequence FASTQ files that can be directly proceed to quantification analysis.
2. Quantification: it supports three well-established and widely-used quantifiers:
  1) Salmon: an **alignment-free quantifier** with **wicked-fast speed** and **comarable accuracy** (https://salmon.readthedocs.io/en/latest/salmon.html).
  2) RSEM: an **alignment-based quantifier** with **high accuracy**. It has been used as **gold standard** in many benchmarking studies (https://github.com/bli25/RSEM_tutorial).
  3) STAR: supports spliced transcript alignment (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). It is recommended by GDC (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline).
3. Summary: it generates a html report of quantification analyis, including alignment statistics, correlation analysis, gene body coverage visualization, etc.

## Tutorial

A detailed tutorial to set up and run this pipeline can be found here: https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/.
