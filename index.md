---
title: Overview
layout: default
nav_order: 1
permalink: /index
---

### Scope of this pipeline

---

Bulk RNA sequencing (RNA-Seq) is a highly sensitive and accurate tool for meansuring expression across the transcriptome. In addition to the transcriptome quantification, RNA-Seq also allows researchers to detect new splicing junctions (e.g. TOPHAP/TOPHAP2-regtools), novel transcripts (e.g. Cufflinks), gene fusion (e.g. STAR-Fusion, Arriba), single nucleotide variants (e.g. STAR-GATK), and other features. **<u>This pipeline is for transcriptome quantification purpose only</u>.**

### The core question

---

**How can we ensure the accuracy of bulk RNA-seq quantification?** In this pipeline, we address this by **applying multiple quantification methods for cross-validation**, thereby increasing the reliability of the results.

The current bulk RNA-Seq quantification methods can be grouped into two categories, **alignment-based** and **alignment-free**, as summarized in the table below. 

|            | Alignment-based methods                                      | Alignment-free methods                                       |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Definition | Methods that quanfity from the **alignments to either transcritpome or genome** (i.e. BAM/SAM files). The coordinates of mapped reads are provided. | Methdos that quantify from **read-transcript matches**. The coordinates of mapped reads are **NOT** provided. |
| Principle  | "Seed and extend" for aligners.<br />Quantifiers vary in rules and weighting methods to count reads. | "*k-mer*s" based indexing;<br />Multiple models to identify read-transcript matches,<br /> e.g. SEMEs for  Salmon, T-DBG for Kallisto |
| Examples   | **Aligner**: Bowtie/Bowtie2, STAR, BWA, HISAT2, TopHat2, et. al.<br />**Quantifier**: RSEM, HTSeq, featureCounts, IsoEM, Cufflinks et. al. | Salmon, Kallisto, Sailfish, Fleximer, RNA-Skim, RapMap, et. al. |
| Accuracy   | High                                                         | a little bit lower or equal                                  |
| Speed      | Slow, a few hours for a typical run                          | Super-fast, a few minutes for a type run                     |

In this pipeline, we provides three quantification methods covering both categories:

- [**Salmon**](https://salmon.readthedocs.io/en/latest/salmon.html): one **wicked-fast** and **highly-accurate** alignment-free method which is recently further enhanced by integrating **selective alignment** and **decoy sequences**
- [**RSEM**](https://github.com/bli25/RSEM_tutorial): the most highly cited alignment-based method which shows the highest accuracy in most benchmarks
- [**STAR**](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf): another **alignment-based** quantifier featured by **spliced transcripts alignment. This is the tool used by [GDC mRNA quantification analysis pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline)**

### Overview

---

![Picture](docs/figures/overview.png)

#### Three stages

- **Preprocessing**: to prepare the standard inputs for quantification analysis
- **Quantification**: to **quantify** the gene and transcript expression using both alignment-free and alignment-based methods
- **Summarization**: to compile the **expression matrices** at both gene and transcript levels, and generate the **quanlity control report**.

### Two species

In this pipeline, we have pre-generated the index libraries for two species: **human** and **mouse**, as listed below. For other species, you will need to generate the index libraries by yourself following the Pipeline Setup tutorial.

| Genome | GENCODE release | Release date | Path                                                         |
| ------ | --------------- | ------------ | ------------------------------------------------------------ |
| hg38   | v48             | 05.2025      | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48 |
| hg19   | v48lift37       | 05.2025      | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release48 |
| mm39   | vM37            | 05.2025      | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM37 |
| mm10   | vM25            | 04.2020      | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25 |

### One environment 

We have managed to complie all tools required in this pipeline into one single conda environment. You can easily set it up following our tutorial.

For St Jude HPC users, it can be easily launched by:

```bash
module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025
```