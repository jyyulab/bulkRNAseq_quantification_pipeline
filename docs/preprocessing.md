---
title: Data Preprocessing
layout: default
nav_order: 2
---

**Why is preprocessing required?** The RNA-Seq data you start with could **vary in format** (e.g., FASTQ, BAM, FASTA) and usually **contain noisy sequences** (e.g., adapters leftovers, poor quality bases and other contaminations). So, we need to pre-process these data to generate the standard-in-format, clean-in-sequence FASTQ files which can be directly proceed to quantification analysis.

### 1. Prepare raw FASTQ files

Usually, you have two FASTQ files (R1, R2) for each paired-end sequencing sample (e.g. Sample1), or one single FASTQ file for each single-end sequencing sample (e.g. Sample2). If so, you are good to move forward.

However, if this is not the case, you will need to generate the raw FASTQ files by yourself:

* If you start with multiple FASTQ files for each mate, usually generated in different lanes, you need to **merge** them:

  ```bash
  ## For paired-end sequencing
  dir_sample3=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample3
  cat $dir_sample3/fqRaw_R1_L001.fq.gz $dir_sample3/fqRaw_R1_L002.fq.gz $dir_sample3/fqRaw_R1_L003.fq.gz $dir_sample3/fqRaw_R1_L004.fq.gz > /your_path/fqRaw_R1.fq.gz
  cat $dir_sample3/fqRaw_R2_L001.fq.gz $dir_sample3/fqRaw_R2_L002.fq.gz $dir_sample3/fqRaw_R2_L003.fq.gz $dir_sample3/fqRaw_R2_L004.fq.gz > /your_path/fqRaw_R2.fq.gz
  
  ## For single-end sequencing
  dir_sample4=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample4
  cat $dir_sample4/fqRaw_L001.fq.gz $dir_sample4/fqRaw_L002.fq.gz $dir_sample4/fqRaw_L003.fq.gz $dir_sample4/fqRaw_L004.fq.gz > /your_path/fqRaw.fq.gz
  ```

  > **NOTE:** The files of lanes **MUST BE IN THE SAME ORDER** between mate 1 and mate 2.

* If you start with BAM/SAM files, usually collected from other sources, you need to **convert** them:

  We previously used `**bedtools bamtofastq**` to convert BAM/SAM to FASTQ files, in which the BAM/SAM files MUST BE SORTED BY NAME. As recommended by GDC, we now move to `**bamtofastq**` integrated in the toolset named biobambam2. This command does not need the input BAM/SAM files sorted. The only thing you need to figure out before running it is that the input BAM/SAM file was generated from single-end or paired-end sequencing. This could be done by `samtools view -c -f 1 input.bam` which counts the matching records. It returns 0 for single-end sequencing or a non-zero integer for paired-end sequencing.

  ```bash
  ## to tell the BAM/SAM files are single- or paired-end
  samtools view -c -f 1 input.bam
  # This command counts the matching records in the bam/sam file, and returns 0 for single-end sequeing. Otherwise, the input bam/sam file is paired-end.
  
  ## For paired-end sequencing
  dir_sample5=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample5
  bamtofastq filename=$dir_sample5/rawBAM.toGenome.bam inputformat=bam gz=1 F=/your_path/fqRaw_R1.fq.gz F2=/your_path/fqRaw_R2.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  
  ## For single-end sequencing
  dir_sample6=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample6
  bamtofastq filename=$dir_sample6/rawBAM.toTranscriptome.bam inputformat=bam gz=1 S=/your_path/fqRaw.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  ```

After the step, you should have two FASTQ files for each paired-end sequencing sample, or one single FASTQ file for each single-end sequencing sample. If so, you are good to move forward.

### 2. Quality control

With the raw FASTQ file(s) prepared, a quality control analysis by FastQC is suggested to: 1) generate a html report which will give you a whole view of the FASTQ file quality; and 2) figure out some key informatioins which will be used in sebsequent analysis, e.g. encoding method of Phred quality score.

```bash
## For paired-end sequencing
fastqc /your_path/fqRaw_R1.fq.gz -o /your_path
fastqc /your_path/fqRaw_R2.fq.gz -o /your_path

## For single-end sequencing
fastqc /your_path/fqRaw.fq.gz -o /your_path
```

This analysis will generate a .html file with all metrics integerated. For more details of this report, please check out here: https://rtsf.natsci.msu.edu/sites/_rtsf/assets/File/FastQC_TutorialAndFAQ_080717.pdf

Here are some key infomations you may pay attention to:

* **Encoding of Phred Quality Score**:  it defines how the quality score of base calling is encoded. There are two options: Phred+33 or **Phred+64**. Phred+33 is denoted as **Sanger / Illumina 1.9** in the html report. This should always be the case for those sequenced in recent years. You will only find Phred+64 encoding (denoted as **Illumina 1.5 or lower**) on older data, which was sequenced several years ago. For more details: https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)
* **Sequenc Length**: The most common values are 46/45, 76/75, 101/100 or 151/150. The alignment of reads shorter than 50 could be tricky.
* **Adapter Content**: FastQC could detect several widely-used (**not all**) adapters, e.g. Illumina Universal Adapter (AGATCGGAAGAG), Illumina Small RNA 3' Adapter (TGGAATTCTCGG), Illumina Small RNA 5' Adapter (GATCGTCGGACT), Nextera Transposase Sequence (CTGTCTCTTATA) and SOLID Small RNA Adapter(CGCCTTGGCCGT). 

### 3. Adapter Trimming

Adapter trimming analysis trims not only the **adapter sequences**, but also the **sequences of unknown or low-quality bases**. It also discards the reads of **too-short length**. So, even though no significant adapter content was found in quality control analysis, it is still highly recommended to perform this analysis to  remove the low-quality reads.

We previous employed **Cutadapt** (https://cutadapt.readthedocs.io/en/stable/guide.html) for adapter trimming. It requires the sequence of adapter(s) and takes ~ two hours for a regular run. Now, we move to the **fastp** (https://github.com/OpenGene/fastp#adapters) which can automatically detect the adapter sequence(s) and trim much faster.

``` bash
## For paired-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw_R1.fq.gz -I /your_path/fqRaw_R2.fq.gz -o /your_path/fqClean_R1.fq.gz -O /your_path/fqClean_R2.fq.gz

## For single-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw.fq.gz -o /your_path/fqClean.fq.gz
```

**<u>Key arguments:</u>**

* <u>**-6/--phred64**: enable it if the input is using Phred+64 encoding. If enabled, fastp will automatically convert the Thread scores from Phread+64 to Phread+33. So the outputs of fastp are always encoded by Phread+33. </u>
* **-w/--thread**: number of threads to use concurrently.
* **-l/--length_required**: the trimmed reads shorter than this value will be discarded. The deault is 15, but 30 is recommended. The shorter reads tend to have multiple alignments, which may affect the quantification accuracy.
* **-q/--qualified_quality_phred:** the quality value that a base is qualified. The default is 15, but 20 is recommended.
* **-n/--n_base_limit**: the read/pair with more N bases will be discarded. The default is 5.

<u>**Key outputs:**</u>

* **fqClean.fq.gz**: FASTQ files with trimmed reads.
* **fastp.html**: an htmal report which summarizes some key matrices, e.g. read counts before and after trimming, insert size, base quality distribution, et. al.
* **fastp.json**: an json report with same matrices. This file provides the number of reads before and after adapter trimming used in the final QC report.

After the adapter timing, your FASTQ files are clean-in-sequence and can be directly proceed to quantification analysis.