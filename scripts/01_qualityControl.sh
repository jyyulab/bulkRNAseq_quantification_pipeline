#BSUB -P RNASeq_pipelines
#BSUB -n 1
#BSUB -M 2000
#BSUB -oo 01_qualityControl.out -eo 01_qualityControl.err
#BSUB -J 01_qualityControl
#BSUB -q standard

## 1. Quality Control of RAQ FASTQ Files by FastQC

fastqc=/research/rgs01/applications/hpcf/apps/fastqc/install/0.11.5/fastqc

# For Single-end Reads
$fastqc sample.raw.fq.gz -o output_dir

# For Paired-end Reads
$fastqc sample_R1.raw.fq.gz -o output_dir
$fastqc sample_R2.raw.fq.gz -o output_dir

# Note: From the results, you need to figure out:
# a. Encoding of Phred Scores: Phred+33 is listed as Illumina 1.9/Sanger, while Phred+64 encoding as illumina 1.5 or lower. (Find more details here: https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)
# b. Sequence Length: The most common values are 46/45, 76/75, 101/100 or 151/150.
# c. Adapter Type: Illumina Universal Adapter(AGATCGGAAGAG), Illumina Small RNA 3' Adapter(TGGAATTCTCGG), Illumina Small RNA 5' Adapter(GATCGTCGGACT), Nextera Transposase Sequence(CTGTCTCTTATA) and SOLID Small RNA Adapter(CGCCTTGGCCGT).

########## If no adaptor is found in the RAW FASTQ files, we are done for this step, and use the RAW FASTQ files in subsequent analysis. Else, we have to trim the adaptors. ##########

## 2. Adaptor Trimming

cutadapt=/hpcf/apps/python/install/3.6.1/bin/cutadapt
fastqc=/research/rgs01/applications/hpcf/apps/fastqc/install/0.11.5/fastqc

# For Single-end Reads
$cutadapt -a AGATCGGAAGAG --trim-n --max-n=0.5 --quality-base=33 -q 30 -m 30 -o sample.clean.fq.gz sample.raw.fq.gz
$fastqc sample.clean.fq.gz -o output_dir

# For Paired-end Reads
$cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --trim-n --max-n=0.5 --quality-base=33 -q 30,20 -m 30 -o sample_R1.clean.fq.gz -p sample_R2.clean.fq.gz sample_R1.raw.fq.gz sample_R2.raw.fq.gz
$fastqc sample_R1.clean.fq.gz -o output_dir
$fastqc sample_R2.clean.fq.gz -o output_dir

# Note:
# 1) -m refers to the minimum length. Trimmed reads that are shorter than this value will be discarded. Usually, we use 20 for raw sequences of <=50nt, and 30 for those of >= 50. Adjust the value to suit your needs.
# 2) -a/-A are used to trim the 3' adapters.
# 3) --trim-n and --max-n are used to trim the sequence of unknown bases. Disable it if necessary.
# 4) --quality-base and -q are used to trim the sequence of low quality. Disable it if necessary.
# 5) For more information about cutadapt, check out here: https://cutadapt.readthedocs.io/en/stable/guide.html.
# 6) We usually run another round of FastQC on CLEAN FASTQ files to make sure the quality is good.
