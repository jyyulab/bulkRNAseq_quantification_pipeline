---
title: Preprocessing
layout: default
nav_order: 1
parent: Full Tutorial
---

# I: Preprocessing

---

## Why is this necessary?

In a pipeline, it's important to **standardize the inputs** for each step. However, in the real cases, the input data can **very in many aspects**, such as **file format** (e.g., FASTQ, BAM/SAM, FASTA), **Phred quality score encoding method** (e.g., Phred+33, Phred+64), and others. And, most of the RNA-seq data quantifiers or aligners require the input data free with adapter and low-quality sequences. So, we need to preprocess the input data to generate **standard-in-format**, **clean-in-sequence** FASTQ files which can be directly proceed to subsequent quantification analysis.

## 1. Data format standardization

Usually, you have two FASTQ files (R1, R2) for each paired-end sequencing sample (e.g. sample1), or one single FASTQ file for each single-end sequencing sample (e.g. sample2). If so, you are good to move to the **Adapter Trimming** step.

> **NOTE**: There is an ***exceedingly rare*** situation that your input data is in **interleaved FASTQ** format. In this format, both **mate1** and **mate2** reads are combined in a single FASTQ file. This mean, though you just have one FASTQ file per sample, the library type is paired-end, not single-end. You can split the interleaved FASTQ file using the command below:
>
> ``` bash
> ## To split an interleaved FASTQ file
> fastp --interleaved_in --in1 interleaved.fq --out1 fqRaw_R1.fq.gz --out2 fqRaw_R2.fq.gz
> ```
>
> 

However, if this is not the case, you will need to generate the raw FASTQ files by yourself:

* If you start with multiple FASTQ files for each mate, usually generated in different lanes, you need to **merge** them:

  ```bash
  ## For paired-end sequencing
  dir_sample3=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample3
  cat $dir_sample3/fqRaw_R1_L001.fq.gz $dir_sample3/fqRaw_R1_L002.fq.gz $dir_sample3/fqRaw_R1_L003.fq.gz $dir_sample3/fqRaw_R1_L004.fq.gz > /path-to-save-outputs/sample3/fqRaw_R1.fq.gz
  cat $dir_sample3/fqRaw_R2_L001.fq.gz $dir_sample3/fqRaw_R2_L002.fq.gz $dir_sample3/fqRaw_R2_L003.fq.gz $dir_sample3/fqRaw_R2_L004.fq.gz > /path-to-save-outputs/sample3/your_path/fqRaw_R2.fq.gz
  
  ## For single-end sequencing
  dir_sample4=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample4
  cat $dir_sample4/fqRaw_L001.fq.gz $dir_sample4/fqRaw_L002.fq.gz $dir_sample4/fqRaw_L003.fq.gz $dir_sample4/fqRaw_L004.fq.gz > /path-to-save-outputs/sample4/fqRaw.fq.gz
  ```

  >  **NOTE:** For the paired-end data, the files of lanes **MUST BE IN THE SAME ORDER** between **mate1** and **mate2**.

* If you start with BAM/SAM files, usually collected from other sources, you need to **convert** them:

  We previously used **`bedtools bamtofastq`** to convert BAM/SAM to FASTQ files, in which the BAM/SAM files MUST BE SORTED BY NAME. As recommended by GDC, we now move to **`bamtofastq`** integrated in the toolset named **`biobambam2`**. This command does not need the input BAM/SAM files sorted. **The only thing you need to figure out before running it is that the input BAM/SAM file was generated from single-end or paired-end sequencing.** This could be done by **`samtools view -c -f 1 input.bam`** which counts the matching records. It returns 0 for single-end sequencing or a non-zero integer for paired-end sequencing.

  ```bash
  ## To tell the BAM/SAM files are single- or paired-end
  samtools view -c -f 1 input.bam
  # This command counts the matching records in the bam/sam file
  # It returns 0 for single-end sequeing. Otherwise, the input bam/sam file is paired-end.
  
  ## For paired-end sequencing
  dir_sample5=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample5
  bamtofastq filename=$dir_sample5/rawBAM.toGenome.bam inputformat=bam gz=1 F=/your_path/fqRaw_R1.fq.gz F2=/your_path/fqRaw_R2.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  
  ## For single-end sequencing
  dir_sample6=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample6
  bamtofastq filename=$dir_sample6/rawBAM.toTranscriptome.bam inputformat=bam gz=1 S=/your_path/fqRaw.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  ```

<u>**Key outputs:**</u>

- For paired-end sequencing samples, there are two FASTQ files per sample: **`fqRaw_R1.fq.gz`** and **`fqRaw_R2.fq.gz`**
- For single-end sequencing samples, there are one single file per sample: **`fqRaw.fq.gz`**

Now, though these files are standard-in-format, as noted in the file name, they are still "raw": they contain noisy sequences, and the Phred quality control encoding method could vary.

## 2. Adapter Trimming

Adapter trimming analysis trims not only the **adapter sequences**, but also the **sequences of unknown or low-quality bases**. It also discards the reads of **too-short length**. So, even though no significant adapter content was found in quality control analysis, it is still highly recommended to perform this analysis to  remove the low-quality reads.

We previously employed **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)** for adapter trimming. It requires the sequence of adapter(s) and takes ~ two hours for a regular run. Now, we move to the **[fastp](https://github.com/OpenGene/fastp#adapters)** which can automatically detect the adapter sequence(s) and trim much faster.

``` bash
## For paired-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw_R1.fq.gz -I /your_path/fqRaw_R2.fq.gz -o /your_path/fqClean_R1.fq.gz -O /your_path/fqClean_R2.fq.gz

## For single-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw.fq.gz -o /your_path/fqClean.fq.gz
```

**<u>Key arguments:</u>**

* **-6/--phred64**: enable it if the input is using Phred+64 encoding. If enabled, fastp will automatically convert the Thread scores from Phread+64 to Phread+33. So the outputs of fastp are always encoded by Phread+33.

  > Not sure about the answer? These two guidelines can help you determine the correct Phred encoding method:
  >
  > - **Phred+64** was retired in late 2011. Data genrated after that time should use **Phred+33**.
  >
  > - You can use ***FastQC*** to identify the Phred encoding of your input files: 
  >
  >   ``` bash
  >   ## To tell the Phred quality score encoding method in FASTQ/BAM/SAM files
  >   fastqc input.fq.gz # for FASTQ files
  >   fastqc input.bam # for BAM/SAM files
  >   # This command generates a html report. In the "Basic Statistics" section, there is a measure called "Endcoding":
  >   # "Sanger / Illumina 1.9" indicates Phred33, while "Illumina 1.5 or lower" indicates Phred64.
  >   ```

* **-w/--thread**: number of threads to use concurrently.

* **-l/--length_required**: the trimmed reads shorter than this value will be discarded. The deault is 15, but 30 is recommended. The shorter reads tend to have multiple alignments, which may affect the quantification accuracy.

* **-q/--qualified_quality_phred:** the quality value that a base is qualified. The default is 15, but 20 is recommended.

* **-n/--n_base_limit**: the read/pair with more N bases will be discarded. The default is 5.

<u>**Key outputs:**</u>

* **fqClean.fq.gz**: FASTQ files with clean reads. And the Phred quality score encoding method is Phred+33.
* **fastp.html**: an HTML report which summarizes some key matrices, e.g. read counts before and after trimming, insert size, base quality distribution, et. al.
* **fastp.json**: an json report with same matrices. This file provides the number of reads before and after adapter trimming used in the final QC report.

Now, the FASTQ files are both standard-in-format and clean-in-sequence. They are ready for subsequent quantification analysis.