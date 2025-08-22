# RNAseq pipeline

Collection of wrapper scripts for RNA-seq pipeline that runs all preprocessing (trimming, mapping, filtering bam files,  QC plots followed by counting reads at both gene and isoform level)

## Dependency

The pipeline automatically loads all the necessary modules in HPCF and submit jobs for all the samples in parallel.  
```bash
trimmomatic/0.36
star/2.5.3a
samtools/1.9
rsem/1.3.1
python/3.7.0
salmon
```
## Config file for mapping indexed reference genomes: 
```
[hg19]
star_path = /research/rgs01/applications/hpcf/apps/star/install/2.5.3a/bin
rsem = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/RSEM/index_star/hg19
salmon_path = /research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/Salmon/bin/salmon
salmon_index = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/Salmon/index_quasi
salmon_annotation = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/idmap.tr2gene.txt
[hg38]
star_path = /research/rgs01/applications/hpcf/apps/star/install/2.5.3a/bin
rsem = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/RSEM/index_star/hg38
salmon_path = /research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/Salmon/bin/salmon
salmon_index = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/Salmon/index_quasi
salmon_annotation = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/idmap.tr2gene.txt
[mm10]
star_path = /research/rgs01/applications/hpcf/apps/star/install/2.5.3a/bin
rsem = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/RSEM/index_star/mm10
salmon_path = /research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/Salmon/bin/salmon
salmon_index = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/Salmon/index_quasi
salmon_annotation = /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/idmap.tr2gene.txt
```

## Step 1

1. Merged fastq.gz files either single end or paired end for all the samples in one directory. 
eg
```
1830518_R1_001.fastq.gz
1830518_R2_001.fastq.gz
1830519_R1_001.fastq.gz
1830519_R2_001.fastq.gz
```
```bash
#Example script to merge fastq files for same sample from multiple lanes

#!/bin/bash
mkdir names
cd names
for file in path_to_fastq/*
do
 touch `echo $file | awk -F/ '{print $NF}'|awk -F_ '{print $1}'` 
done
mkdir ../merged
cd ..
for i in names/*
do
 file=$(basename "$i")
 echo "merging file $file"
 cat path_to_fastq/"$file"_*_*_R1_*.fastq.gz > merged/"$file"_R1_001.fastq.gz
 cat path_to_fastq/"$file"_*_*_R2_*.fastq.gz > merged/"$file"_R2_001.fastq.gz
done
```
2. After you have merged the fastq files, run the main script job_runRNAseq.py and follow --help for documentation: by default, salmon is set as method for transcript quantification

```bash
python /home/ydhungan/pipelines/rnaseq_pipeline/bin/job_runRNAseq.py --help
usage: Wrapper for BSUB job submission for RNA-seq data. [-h]
                                                         [--path-to-fastqs PATH_TO_FASTQS]
                                                         [--paired-end]
                                                         [--configfile CONFIGFILE]
                                                         [--adapters ADAPTERS]
                                                         [--genome {mm10,hg19,hg38}]
                                                         [--method {salmon,rsem}]
                                                         [--strandedness {none,forward,reverse}]
                                                         [--trim-status {1,0}] # 1 = perform trimming; 0 = no trimming
                                                         [--memory MEMORY]
                                                         [--queue QUEUE]
                                                         [--out-dir OUT_DIR]
                                                         [--log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}]

optional arguments:
  -h, --help            show this help message and exit
  --path-to-fastqs PATH_TO_FASTQS
                        Path to FASTQ files.
  --paired-end          Paired-end fastqs.
  --configfile CONFIGFILE
                        Configuration file with information for genome
                        database
  --adapters ADAPTERS   Fasta file of adapter sequences
  --genome {mm10,hg19,hg38}
                        Genome to use for quantification
  --method {salmon,rsem}
                        method of RNA-seq quantification
  --strandedness {none,forward,reverse}
                        method of RNA-seq quantification
  --trim-status {1,0}   Option to trim or not to trim the input reads
  --memory MEMORY       Memory requested to run the analysis.
  --queue QUEUE         Queue to submit the job in HPCF (use bqueues to
                        choose).
  --out-dir OUT_DIR     Output Directory.
  --log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}
                        Log level



```
This script will create analysis directory under provided output directory using first part of the fastq.gz file name eg 1830518_R1_001.fastq.gz and dump all the intermediate files in that directory. The quantification at gene and isoform level can be imported for all the samples using tximport() function in R. 


