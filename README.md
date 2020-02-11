# RNASeq_Quantification_Pipelines

## Overview
![Picture1](https://user-images.githubusercontent.com/33663247/66729168-3ef9f700-ee0f-11e9-866a-d4a621466c9d.png)


### 0. Preprocess to prepare the RAW FASTQ

* If you have multiple fastq files for each sample, usually generated from different lanes, you need to merge them into one.
* If you start from BAM files, usually downloaded from some databases, you need to convert them into FASTQ files.


### 1. Quality Control of RAQ FASTQ Files by FastQC

Note: From this analysis, we need to figure out:
* a. **Encoding of Phred Scores**: Phred+33 is listed as Illumina 1.9/Sanger, while Phred+64 encoding as illumina 1.5 or lower. (Find more details here: https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)
* b. **Sequence Length**: The most common values are 46/45, 76/75, 101/100 or 151/150.
* c. **Adapter Type**: Illumina Universal Adapter(AGATCGGAAGAG), Illumina Small RNA 3' Adapter(TGGAATTCTCGG), Illumina Small RNA 5' Adapter(GATCGTCGGACT), Nextera Transposase Sequence(CTGTCTCTTATA) and SOLID Small RNA Adapter(CGCCTTGGCCGT).

Here is a nice tutorial for FastQC: https://www.youtube.com/watch?v=bz93ReOv87Y

If no adaptor is found in the RAW FASTQ files, we are done for this step, and use the RAW FASTQ files in subsequent analysis. Otherwise, we have to trim the adaptors.

### 2. Quantification by Salmon

* Salmon is alignment-free and hence ultra-fast!  
* Salmon is easy to use: there is only a few options you need to specify. Salmon could figure some options out by itself.
* Just provide the mapping file of transcripts to genes, then it will generate the quantification results of both transcripts and genes.
* Salmon predicts the library type by default. If you don't the library type of your samples, you could use it to figure out.

### 3. Quantification by RSEM

* RSEM is a well-accepted gold standard for RNA-Seq quantification.
* RSEM is a alignment-based quantification method, which makes it a little bit more complicated to use: you have to specify some options.
* RSEM generates the BAM file of transcriptomic alignment by default, and could also generate the one of genomic alignment by specify the corresponding arguments.

### 4. Quantifcation by STAR-HTSeq Strategy

* STAR-HTSeq strategy is recommended by GDC.
* The 2-pass STAR alignment is famous for its speed and accuracy.
* HTSeq is very popular to quantify the expression of genes.

### 5. Gene Body Coverage Analysis

* This analysis plots the distributions of reads along the tanscripts/gens. It is used to evaluate the quality of sample library, especially the status of RNA degradation.

### 6. Quantification Summary

* The expression matrix of all samples under a certain quantification method is generated
* The correlation between quantification methods is calculated.
