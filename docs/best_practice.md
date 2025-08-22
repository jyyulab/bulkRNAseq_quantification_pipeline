---
title: Best practice
layout: default
nav_order: 6
permalink: /docs/best_practice
---

The step-by-step tutorial shown above is set to get you familiar with this pipeline. In practice, it could be tedious to run this pipeline in a step-by-step manor, especially in the cases with multiple samples. For the best practice, we designed **a sample table-centered strategy** to run this pipeline.

### 1. Sample Table

As shown below, sample table is a tab-delimited text file with 6 columns:

* **sampleID**: name of samples, which can only contain letters, numbers and underscores. They should NOT start with numbers.
* **libraryType**: [PE | SE]. For the alignment files, like BAM or SAM, please use samtools view -c -f 1 input.bam to tell. See above for how.
* **phredMethod**: [Phred33 | Phred64]. Those data generated in recent 5 years should be Phred33. FastQC can tell it. See above for how.
* **reference**: [hg38 | hg19 | mm39 | mm10]. Reference genomes supported by this pipeline.
* **input**: files for quantification. This pipeline supports FASTQ, BAM and SAM.
* **output**: path to save the output files.

![image-20230901180916287](/Users/qpan/Library/Application Support/typora-user-images/image-20230901180916287.png)

That's it! This table is the only one you need to prepare. You can generate it in either of these ways:

* any coding language you prefer, e.g. BASH, R, Python, Perl et. al.
* Excel or VIM. And for you convince, we have a templete avaible on HPC: /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/temp_bulkRNAseq.txt. You can simply copy it to your folder and edit it in VIM. (To insert TAB in VIM: on INSERT mode press control + v + TAB)

### 2. To run the pipeline

```bash
## 0. activate the conda env
module load conda3/202210
conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023

## 1. preprocessing
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/preProcessing.pl sampleTable.txt

## 2. adapterTrimming
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/adapterTrimming.pl sampleTable.txt

## 3. quantify by Salmon
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/quantSalmon.pl sampleTable.txt

## 4. genebody coverage
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/geneCoverage.pl sampleTable.txt

## 5. quantification summary
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/quantSummary.pl sampleTable.txt

```

