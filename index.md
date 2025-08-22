---
title: Overview
layout: default
nav_order: 1
permalink: /index
---

Bulk RNA sequencing (RNA-Seq) is a highly sensitive and accurate tool for meansuring expression across the transcriptome. In addition to the transcriptome quantification, RNA-Seq also allows researchers to detect new splicing junctions (e.g. TOPHAP/TOPHAP2-regtools), novel transcripts (e.g. Cufflinks), gene fusion (e.g. STAR-Fusion, Arriba), single nucleotide variants (e.g. STAR-GATK), and other features. **<u>This pipeline is for transcriptome quantification purpose only</u>.**

The current bulk RNA-Seq quantification methods can be grouped into two categories, **alignment-based** and **alignment-free**, as summarized in the table below. 

|            | Alignment-based methods                                      | Alignment-free methods                                       |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Definition | Methods that quanfity from the **alignments to either transcritpome or genome** (i.e. BAM/SAM files). The coordinates of mapped reads are provided. | Methdos that quantify from **read-transcript matches**. The coordinates of mapped reads are **NOT** provided. |
| Principle  | "Seed and extend" for aligners.<br />Quantifiers vary in rules and weighting methods to count reads. | "*k-mer*s" based indexing;<br />Multiple models to identify read-transcript matches,<br /> e.g. SEMEs for  Salmon, T-DBG for Kallisto |
| Examples   | **Aligner**: Bowtie/Bowtie2, STAR, BWA, HISAT2, TopHat2, et. al.<br />**Quantifier**: RSEM, HTSeq, featureCounts, IsoEM, Cufflinks et. al. | Salmon, Kallisto, Sailfish, Fleximer, RNA-Skim, RapMap, et. al. |
| Accuracy   | High                                                         | a little bit lower or equal                                  |
| Speed      | Slow, a few hours for a typical run                          | Super-fast, a few minutes for a type run                     |

<u>**To ensure the accuracy of quantification, we employ one signature method from each of these two categories for cross-validation**:</u> **1)** **RSEM**, the most highly cited alignment-based method which shows the highest accuracy in most benchmarks; **2)** **Salmon**, one wicked-fast and highly-accurate alignment-free method which is recently further enhanced by integrating selective alignment and decoy sequences. We also introduce **3) STAR**, another alignment-based method recomended by GDC, as an optional method.

Below is an overview of the pipelines, which contains three sections: 1) **Preprocessing**: to prepare the standard inputs for quantification analysis; 2) **Quantification**: to **quantify** the gene and transcript expression in both alignment-based and alignment-free methods; 3) **Summarization**: to compile the **expression matrices** at both gene and transcript levels, and generate the **quanlity control report**.

![Picture1](/Users/qpan/Desktop/Picture1.png)

To serve better, we have:

* Complied all tools required into one single conda environment, which can be easily launched by:

  ```bash
  module load conda3/202210
  conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023
  ```

* Updated all the index libraries to the latest version

  | Genome | GENCODE | Path                                                         |
  | ------ | ------- | ------------------------------------------------------------ |
  | hg38   | v43     | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43 |
  | hg19   | v43     | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43 |
  | mm39   | vM32    | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32 |
  | mm10   | vM25    | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25 |

* Prepared the test cases of multi-format inputs

  | Cases   | Library Type | File Format           | Path                                                         |
  | ------- | ------------ | --------------------- | ------------------------------------------------------------ |
  | Sample1 | Paired-end   | FASTQ                 | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample1 |
  | Sample2 | Single-end   | FASTQ                 | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample2 |
  | Sample3 | Paired-end   | FASTQ, multiple lanes | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample3 |
  | Sample4 | Single-end   | FASTQ, multiple lanes | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample4 |
  | Sample5 | Paired-end   | BAM                   | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample5 |
  | Sample6 | Single-end   | BAM                   | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample6 |

For traning purpose, we will go through the pipelines in a step-by-step way with real cases as listed above. 