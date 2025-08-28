---
title: Full Tutorial
layout: default
nav_order: 4
permalink: /docs/3_full_tutorial/index
---

### Step-by-Step Tutorial

---

In this section, we provide a **step-by-step tutorial** to run this pipeline. This is mainly for learning purpose.
For each step, we will introduce:

- why we need it?
- the specific commands to run it and the key arguments
- key outputs

### A quick run of this pipeline

---

In practice, you usually have multiple samples to analyze, and it would be tedious and less efficient to do it in a line-by-line manor. Instead, for the best practice, we desigend a sample table-centered strategy to run this pipeline quickly and easily.

#### I. Prepare the sample table

Below is an example of the sample table for this pipeline:

![image](sampleTable_template.png)

  It is a **tab-delimited text file with 6 columns**:

1) **<u>sampleID</u>**: name of samples. Some rules apply:
   - Should contain **letters**, **numbers** or **underscores** only
   - Should NOT **start with numbers**

2. **<u>libraryType</u>**: type of libraries, `[PE | SE]`. 

   - My input files are in BAM/SAM format, and I'm not sure about the library type of them? The command below can tell that:

     ``` bash
     ## To tell the BAM/SAM files are single- or paired-end
     samtools view -c -f 1 input.bam
     # This command counts the matching records in the bam/sam file.
     # It returns 0 for single-end sequeing. Otherwise, the input bam/sam file is paired-end.
     ```

3. **<u>phredMethod</u>**: Phred quality score encoding method, `[Phred33 | Phred64]`. Not sure about the answer? These two rules can help tell that:

   - Phred64 was retired in late 2011. Data genrated after that should be in Phred33.

   - Use the FastQC to tell that: 

     ``` bash
     ## To tell the Phred quality score encoding method in FASTQ/BAM/SAM files
     fastqc input.fq.gz # for FASTQ files
     fastqc input.bam # for BAM/SAM files
     # This command generates a html report. In the "Basic Statistics" section, there is a measure called "Endcoding":
     # "Sanger / Illumina 1.9" indicates Phred33, while "Illumina 1.5 or lower" indicates Phred64.
     ```

4. **<u>reference</u>**: reference genome assembly, `[hg38 | hg19 | mm39 | mm10]`. 
   - We recommend to use `hg38` for human samples, and `mm39` for mouse samples. The assembly, `hg19` and `mm10`, are mainly used to match some legacy data.
   - For other species or genome assembly, you will need to manually create the required reference files following [this tutorial](https://jyyulab.github.io/bulkRNAseq_quantification_pipeline/docs/pipeline_setup/reference.html).

5. **<u>input</u>**: input files for quantification. This pipeline accepts:

   - **<u>Standard FASTQ files</u>**: both paired-end (sample1) and single-end (sample2).

   - **<u>FASTQ files of multple lanes</u>**: both paired-end (sample3) and single-end (sample4).

   - **<u>BAM/SAM files</u>**: single file only (sample 5 and sample6). For samples with  multiple BAM/SAM files (usually splited ones), please merge them first:

     ``` bash
     samtools merge -o output.bam input_1.bam input_2.bam ...
     ```

6. **<u>output</u>**: path to save the output files. This pipeline will create a folder named by the `sampleID` under this directory.

That's it! This table is the only one you need to prepare. You can generate it in either of these ways:

* Any coding language you prefer, e.g. BASH, R, Python, Perl *et. al.*
* Excel or VIM. And for you convince, we have a templete avaible on HPC: /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/git_repo/testdata/sampleTable.testdata.txt. You can simply copy it to your folder and edit it in VIM. (To insert TAB in VIM: on INSERT mode press control + v + TAB)



### II. Run the script of each step



``` bash
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





It's mainly for learning purpose, and you are welcome to go throught them line by line, mainly for learning purposes. In practice, you usually have multiple samples to run, and it could be tedious to run  

There are two parts to set up this pipeline:

- **Software installation**: to install all tools required by this pipeline.
- **Reference preparation**: to generate reference files used in this pipeline.

For Yu Lab members, we have set this pipeline up on HPC:

- All the tools have been compiled in one conda environment, which can be launched by:

  ``` bash
  module load conda3/202402
  conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025
  ```

- We have generated the files for four reference genomes: hg38, hg19, mm39 and mm10.


  | Genome            | GENCODE release | Release date | Ensembl release | Path                                                         |
  | ----------------- | --------------- | ------------ | --------------- | ------------------------------------------------------------ |
  | hg38 (GRCh38.p14) | v48             | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48 |
  | hg19 (GRCh37.p13) | v48lift37*      | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release48 |
  | mm39 (GRCm39)     | vM37            | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM37 |
  | mm10 (GRCm38.p6)  | vM25            | 04.2020**    | v100            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25 |

  *: The updates for the hg19/GRCh37 genome assembly have stopped in 2013. However, gene annotation continue to be updated by mapping the comprehensive gene annotations originally created for the GRCh38/hg38 reference chromosomes onto GRCh37 primary assembly using [gencode-backmap](https://github.com/diekhans/gencode-backmap) .

  **: The updates for the mm10/GRCm38 genome assembly and gene annotation have stopped in 2019.

You should run this tutorial only when you want to set up this pipeline locally or complie other reference genome / gene annotation.