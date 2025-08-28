---
title: II-Database Preparation
layout: default
nav_order: 2
parent: Pipeline Setup
---

### Part II: Reference preparation

In addition to the tools, you will also need to prepare **reference genomes** for alignment, quantification and QC assessment. Below is a summary of reference preparation (use hg38 as an example):

![Picture](referencePreparation.png)

#### 1. Data collection

The reference preparation stars with FOUR files that can be  directly downloaded from websites:

- **Gene Annotation file in [GTF](https://biocorecrg.github.io/PhD_course_genomics_format_2021/gtf_format.html) (Gene Transfer Format) format**: e.g., /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.primary_assembly.annotation.gtf

- **Genome sequence file in [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) format**: e.g., /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa

- **Transcriptome sequence file in [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) format**: e.g., `/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.transcripts.fa`.

  ***<u>NOTE:</u>*** For the three files above, they are usually available at open-sourced websites. For human and mouse, we recommend [GENCODE](https://www.gencodegenes.org/) to collect them, while for other species, we recommend [Ensembl](https://useast.ensembl.org/info/data/ftp/index.html). 

- **HouseKeeping gene list**: the housekeeping genes defined by [this study](https://www.sciencedirect.com/science/article/pii/S0168952513000899?via%3Dihub) (N = 3804), e.g., /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/housekeeping_genes.human.txt.

  ***<u>NOTE:</u>*** For other species, you can generate the housekeeping gene list by gene homology conversion using BiomaRt or other tools. Below is the codes I used to generate the housekeeping genes in mouse:

  ``` R
  library(NetBID2)
  
  HK_hg <- read.table("/Volumes/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/housekeeping_genes.human.txt")
  HK_mm <- get_IDtransfer_betweenSpecies(
    from_spe = "human", to_spe = "mouse", from_type = "hgnc_symbol", to_type = "mgi_symbol",
    use_genes = unique(HK_hg$V1));
  colnames(HK_mm) <- paste0("#", colnames(HK_mm))
  write.table(HK_mm[,c(2,1)], ## Please make sure the mouse gene symbols are in the FIRST column
              file = "/Volumes/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM37/housekeeping_genes.mouse.txt", col.names = T, row.names = F, sep = "\t", quote = F)
  ```



#### 2. Parsing annotation file

In this step, we will parse the gene annotation file and generate four files that required in downstream analysis:

- gencode.v48.primary_assembly.annotation.gene2transcript.txt and gencode.v48.primary_assembly.annotation.transcript2gene.txt: the mappings between transcripts and genes. They are required in gene-level quantification.
- gencode.v48.primary_assembly.annotation.geneAnnotation.txt and gencode.v48.primary_assembly.annotation.transcriptAnnotation.txt: the gene or transcript annotation file. They are required in the final gene expression matrix generation.

To make the parsing analysis easiler, we created a script, `parseAnnotation.pl`, that you can easily generate the four files wit h it:

``` bash
## parse the gene anotation file
## This command will generate the four files in the same folder as gencode.v48.primary_assembly.annotation.gtf.
## Only ONE argument is need: gene annotation file in GTF format.
perl /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/git_repo/scripts/setup/parseAnnotation.pl gencode.v48.primary_assembly.annotation.gtf
```



#### 3. Creating gene body bins

In this step, we will create a bin list for the longest transcript of each gene, with 100 bins per transcript by default. This list is required in genebody coverage statistics - an important QC metrics that indicates the extent of RNA degradation.

Two files will be generated:

- ./bulkRNAseq/genebodyBins/genebodyBins_allTranscripts.txt: bins of the longest transcripts of all genes. This one is the most reliable solution since it calculates the gene body coverage across all genes (N = 46,402 for human). However, it's much slower.
- ./bulkRNAseq/genebodyBins/genebodyBins_HouseKeepingTranscripts.txt: bins of the longest transcripts of precurated housekeeping genes. Though it only considers the housekeeping genes (N = 3,515 in human), based on our tests across 30+ datasets, no significant difference of gene coverate statistics was observed compared to the all-transcript version. And it's way faster. This is widely-used in many pipelines, including the [RseQC](https://rseqc.sourceforge.net/#genebody-coverage-py). So, we set it as the default in this pipeline.

You can easily generate these two files with the command below:

``` bash
## create the gene body bins
## This command will generate the two files containing the bin list of the longest transcript of all genes and housekeeping genes.
## Three arguments are needed: transcriptome sequence file in FASTA format, a txt file containiing housekeeping genes in the first column, and a directory to save the output files.
perl /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/git_repo/scripts/setup/createBins.pl gencode.v48.transcripts.fa housekeeping_genes.human.txt ./bulkRNAseq/genebodyBins
```



#### 4. Create genome index files for RSEM

``` bash
#BSUB -P buildIndex
#BSUB -n 8
#BSUB -M 8000
#BSUB -oo 01_buildIndex.out -eo 01_buildIndex.err
#BSUB -J buildIndex
#BSUB -q priority

rsem_index=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/bin/rsem-prepare-reference

## Bowtie2-RSEM
mkdir /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/RSEM/index_bowtie2
$rsem_index --gtf /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.primary_assembly.annotation.gtf --bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/bin --num-threads 8 /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/RSEM/index_bowtie2/hg38

## STAR-RSEM
mkdir /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/RSEM/index_star
$rsem_index --gtf /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.primary_assembly.annotation.gtf --star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025/bin --num-threads 8 --star-sjdboverhang 100 /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/RSEM/index_star/hg38
```



#### 5. Create genome index files for Salmon

``` bash
#BSUB -P salmonIndex
#BSUB -n 8
#BSUB -M 8000
#BSUB -oo 01_buildIndex.out -eo 01_buildIndex.err
#BSUB -J buildIndex
#BSUB -q standard

## generate a decoy-aware transcriptome
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 | cut -d ">"     -f 2 > /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/Salmon/decoys.txt

cat /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.transcripts.fa /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa > /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/Salmon/gentrome.fa

salmon index -t /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/Salmon/gentrome.fa -d /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/Salmon/decoys.txt -p 8 -i index_decoy --gencode -k 31
```



#### 6. Create genome index files for STAR

``` bash
#BSUB -P STAR_Index
#BSUB -n 8
#BSUB -M 8000
#BSUB -oo 01_buildIndex.out -eo 01_buildIndex.err
#BSUB -J buildIndex
#BSUB -q standard

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/bulkRNAseq/STAR/index_overhang100 --genomeFastaFiles /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48/gencode.v48.primary_assembly.annotation.gtf --sjdbOverhang 100
```

