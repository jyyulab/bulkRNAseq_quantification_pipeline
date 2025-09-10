---
title: III-Test Data Collection
layout: default
nav_order: 3
parent: Pipeline Setup

---

## Part III: Test Data Collection (optional)

We orginally curated this test data, contaiining 6 samples of different conditions (see table below), for the purpose of **testing and troubleshooting** this pipeline.   We also use them for **showcasing the input formats** accepted by this pipeline, and  **benchmarking pipeline performance** (e.g., Run Time, CPU Time and Max Memory Usage), though these six samples were prepared by downsampling from real data.

|         | Format                     | Library Type | Phred Encoding | Species | Number of Reads | Preparation                |
| ------- | -------------------------- | ------------ | -------------- | ------- | --------------- | -------------------------- |
| sample1 | FASTQ                      | Paired-end   | Phred+33       | human   | 19,410,373 * 2  | Downsampled from real data |
| sample2 | FASTQ                      | Single-end   | Phred+33       | mouse   | 16,005,450      | Downsampled from real data |
| sample3 | FASTQ, multiple lanes      | Paired-end   | Phred+33       | human   | 19,410,373 * 2  | Split by lane from sample1 |
| sample4 | FASTQ, multiple lanes      | Single-end   | Phred+33       | mouse   | 16,005,450      | Split by lane from sample2 |
| sample5 | BAM/SAM (to genome)        | Paired-end   | Phred+33       | human   | 13,856,075 * 2  | Downsampled from real data |
| sample6 | BAM/SAM (to transcriptome) | Single-end   | Phred+33       | mouse   | 13,533,162      | Downsampled from real data |

This is only one situation that you need collect the test data: **<u>you would like to walk through this pipeline but do not have your own bulk RAN-seq available</u>**. In this case, you can download the test data for one or more samples to get started.

- For those who can access to the pre-built conda environment, the test data is available at: `/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/pipeline/testdata`. 
- For other users, please download them from [Zenodo](https://zenodo.org/records/17088549).