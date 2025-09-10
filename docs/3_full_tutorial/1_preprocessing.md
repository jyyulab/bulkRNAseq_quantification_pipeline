---
title: Preprocessing
layout: default
nav_order: 1
parent: Full Tutorial
---

# I: Preprocessing

---

## Why is Preprocessing necessary?

In any pipeline, **standardizing input data** is crucial for ensuring smooth and accurate downstream analysis. In practice, input data can **vary widely** in terms of **file format** (e.g., FASTQ, BAM/SAM, FASTA), **Phred quality score encoding method** (e.g., Phred+33, Phred+64), and other attributes. Additionally, most RNA-seq quantification tools, such as quantifiers and aligners, require input data to be free from **adapter contamination** and **low-quality sequences**. So, we need to preprocess the input data to generate **standard-in-format**, **clean-in-sequence** FASTQ files which can be directly proceed to subsequent quantification analysis.

## 1. Data format standardization

Typically, paired-end sequencing samples have two FASTQ files (R1 and R2, e.g., test data: sample1), while single-end samples have one (e.g., test data: sample2). If your data matches this format,  you can proceed directly to **Adapter Trimming**.

> ***NOTE:***
>
> In ***exceedingly rare*** cases, your input data may be an **interleaved FASTQ** file (both **mate1** and **mate2** reads combined in a single file). In this case, although there is only one file per sample, the library type is paired-end. You can split interleaved FASTQ files using the command below:
>
> ``` bash
> ## To split an interleaved FASTQ file
> fastp --interleaved_in --in1 interleaved.fq --out1 fqRaw_R1.fq.gz --out2 fqRaw_R2.fq.gz
> ```
>

However, if your data does not match the standard format, you will need to generate the raw FASTQ files by yourself:

* If you start with multiple FASTQ files for each mate, usually generated in different lanes, you need to **merge** them:

  ```bash
  ## For paired-end sequencing
  dir_sample3=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample3 # path to the test data: sample3
  cat $dir_sample3/fqRaw_R1_L001.fq.gz $dir_sample3/fqRaw_R1_L002.fq.gz $dir_sample3/fqRaw_R1_L003.fq.gz $dir_sample3/fqRaw_R1_L004.fq.gz > /path-to-save-outputs/sample3/fqRaw_R1.fq.gz
  cat $dir_sample3/fqRaw_R2_L001.fq.gz $dir_sample3/fqRaw_R2_L002.fq.gz $dir_sample3/fqRaw_R2_L003.fq.gz $dir_sample3/fqRaw_R2_L004.fq.gz > /path-to-save-outputs/sample3/your_path/fqRaw_R2.fq.gz
  
  ## For single-end sequencing
  dir_sample4=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample4 # path to the test data: sample4
  cat $dir_sample4/fqRaw_L001.fq.gz $dir_sample4/fqRaw_L002.fq.gz $dir_sample4/fqRaw_L003.fq.gz $dir_sample4/fqRaw_L004.fq.gz > /path-to-save-outputs/sample4/fqRaw.fq.gz
  ```

  >  ***NOTE:***
  >
  >  For paired-end data, the files of lanes **MUST BE IN THE SAME ORDER** between **mate1** and **mate2**.

* If you start with BAM/SAM files, usually collected from other sources, you need to **convert** them:

  We previously used **`bedtools bamtofastq`** to convert BAM/SAM to FASTQ files, in which the BAM/SAM files MUST BE SORTED BY NAME. As recommended by GDC, we now move to **`bamtofastq`** integrated in the toolset named **`biobambam2`**. This command does not need the input BAM/SAM files sorted. **The only thing you need to figure out before running it is that the input BAM/SAM file was generated from single-end or paired-end sequencing.** This could be done by **`samtools view -c -f 1 input.bam`** which counts the matching records. It returns 0 for single-end sequencing or a non-zero integer for paired-end sequencing.

  ```bash
  ## To tell the BAM/SAM files are single- or paired-end
  samtools view -c -f 1 input.bam
  # This command counts the matching records in the bam/sam file
  # It returns 0 for single-end; non-zero for paired-end.
  
  ## For paired-end sequencing
  dir_sample5=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample5 # path to the test data: sample5
  bamtofastq filename=$dir_sample5/rawBAM.toGenome.bam inputformat=bam gz=1 F=/your_path/fqRaw_R1.fq.gz F2=/your_path/fqRaw_R2.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  
  ## For single-end sequencing
  dir_sample6=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata/sample6 # path to the test data: sample6
  bamtofastq filename=$dir_sample6/rawBAM.toTranscriptome.bam inputformat=bam gz=1 S=/your_path/fqRaw.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  ```

<u>**Key outputs:**</u>

- **`fqRaw_R1.fq.gz`** and **`fqRaw_R2.fq.gz`** for paired-end samples.
- **`fqRaw.fq.gz`** for single-end samples

These files are now standard in format, but may still contain adapters and low-quality reads.

## 2. Adapter Trimming

Adapter trimming analysis removes **adapter sequences**, **unknown or low-quality bases**, and **short reads**. It also discards the reads of **too-short length**. So, even though no significant adapter contaminations was found any quality control analysis (e.g., [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) analysis), trimming is recommended to improve data quality.

We previously employed **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)** for adapter trimming. It requires the sequence of adapter(s) and takes ~ two hours for a regular run. Now, we move to **[fastp](https://github.com/OpenGene/fastp#adapters)** which automatically detects the adapter sequence(s) and trims much faster.

``` bash
## For paired-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw_R1.fq.gz -I /your_path/fqRaw_R2.fq.gz -o /your_path/fqClean_R1.fq.gz -O /your_path/fqClean_R2.fq.gz

## For single-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw.fq.gz -o /your_path/fqClean.fq.gz
```

**<u>Key arguments:</u>**

* **`-6/--phred64`**: Enable it if the input is using **`Phred+64`** encoding. Otherwise, **`Phred+33`** would be used by default.

  If this is enabled, **fastp** will automatically convert to **`Phread+33`**. So, the outputs of **fastp** are **ALWAYS** in **`Phread+33`**.

  > ***TIPS:***
  >
  > Not sure about the answer? These two guidelines can help you determine the correct Phred encoding method:
  >
  > - **`Phred+64`** was retired in late 2011. Data genrated after that time should use **`Phred+33`**.
  >
  > - You can use ***[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*** to identify the Phred encoding of your input files: 
  >
  >   ``` bash
  >   ## To tell the Phred quality score encoding method in FASTQ/BAM/SAM files
  >   fastqc input.fq.gz # for FASTQ files
  >   fastqc input.bam # for BAM/SAM files
  >   # This command generates a html report. In the "Basic Statistics" section, there is a measure called "Endcoding":
  >   # "Sanger / Illumina 1.9" indicates Phred33, while "Illumina 1.5 or lower" indicates Phred64.
  >   ```

* **`-w/--thread`**: Number of threads to use concurrently.

* **`-l/--length_required`**: Minimum read length required. Any trimmed reads shorter than this value would be discarded.

  > ***NOTE:***
  >
  > The default value by **fastp** is 15, but **we use 30 in this pipeline**. The shorter reads tend to have multiple alignments, which may affect the quantification accuracy.

* **`-q/--qualified_quality_phred`**: Minimum base quality required. The default value by **fastp** is 15, but **we use 20 in this pipeline**.

* **`-n/--n_base_limit`**: Maximum unknown (**`N`**) bases allowed per read. Any reads with more **`N`** bases would be discarded. The default value is **5**.

<u>**Key outputs:**</u>

* FASTQ files containing clean reads only: **`fqClean_R1.fq.gz`** and **`fqClean_R2.fq.gz`** for paired-end samples, and **`fqClean.fq.gz`** for single-end samples. The Phred quality score encoding method for these files is **`Phred+33`**.
* **`fastp.html`**: an HTML report summarizing key matrices of adapter trimming and others (e.g., read counts before and after trimming, insert size, base quality distribution)
* **`fastp.json`**: an JSON report with same matrices as fastp.html. This file is used in Summarization Report generation.

After completing these two steps, your FASTQ files are both standard-in-format and clean-in-sequence. They are ready for subsequent quantification analysis.