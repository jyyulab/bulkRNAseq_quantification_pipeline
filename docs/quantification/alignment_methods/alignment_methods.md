---
title: Alignment-based methods
layout: default
nav_order: 2
permalink: /docs/quantification/alignment_methods/alignment_methods
parent: Quantification
---

RSEM (RNA-Seq by Expectation Maximization) is an accurate and user-friendly software tool for quantifying transcript and gene abundances from RNA-seq data. Here is the detailed introduction of RSEM: https://github.com/bli25/RSEM_tutorial.

* RSEM is a quantifier, not an aligner. RSEM can directly take the FASTQ files, but it does not align the reads by itself. Instead, it empolys the Bowtie2 (by default) or STAR/HISAT2 (optional) for reads alignment.
* RSEM doesn't rely on the existence of a reference genome. It quantifies from the alignments to reference transcriptome.
* RSEM is famouse for its ability to **effectively use ambiguously-mapping reads**. This is the main reason for its high quantification accuracy.

### 1. Bowtie2-RSEM (default)

In this pipeline, Bowtie2-RSEM is set as the default. 

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/RSEM
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/RSEM
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/RSEM
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/RSEM

## Transcript bins
bins_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/genebodyBins
bins_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/genebodyBins
bins_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/genebodyBins
bins_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/genebodyBins

## For paired-end sequencing
rsem-calculate-expression --num-threads 8 \
--bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--bowtie2-sensitivity-level sensitive \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
--paired-end /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz \
$index_hg38v43/index_bowtie2 /your_path/quantRSEM

## For single-end sequencing
rsem-calculate-expression --num-threads 8 \
--bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--bowtie2-sensitivity-level sensitive \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
/your_path/fqClean.fq.gz \
$index_hg38v43/index_bowtie2 /your_path/quantRSEM

## Prepare gene body coverage statistics
bedtools multicov -bams /your_path/quantRSEM.transcript.sorted.bam -bed $bins_hg38v43/genebodyBins_HouseKeepingTranscripts.txt > /your_path/genebodyCoverage.txt
```

**<u>Key arguments:</u>**

There is only one key argument you need to manually specify, `**--strandedness [none|forward|reverse]**`. You can figure it out from the **lib_format_counts.json** returned by Salmon, as mentioned before.

**<u>Key outputs</u>**: (for more details: https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)

* **quantRSEM.isoforms.results**: contains transcript-level abundance estimates, raw counts, TPM and FPKM.
* **quantRSEM.genes.results**: contains the aggregated gene-level abundance estimates, raw counts, TPM and FPKM.
* **quantRSEM.transcript.bam/quantRSEM.transcript.sorted.bam**: BAM files with alignments to reference transcriptome. These files will be used to calculate gene body coverage and others.
* **quantRSEM.stat/quantRSEM.cnt**: contains alignment statistics based purely on the alignment results obtained from aligners. This file provides the number of totally-mapped and uniquely-mapped reads used in the final QC report.

### 2. STAR-RSEM (optional)

Following the RSEM tutorial, we also provide the codes for STAR-RSEM **as optional**. STAR (Spliced Transcripts Alignment to a Reference) was uniquely designed for RNA-Seq alignment. It's ultrafast and it does support the spliced alignment. For more details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/.

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/RSEM
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/RSEM
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/RSEM
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/RSEM

## For paired-end sequencing
rsem-calculate-expression -num--threads 8 \
--star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--star-gzipped-read-file \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
--paired-end /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz \
$index_hg38v43/index_star /your_path/quantRSEM

## For single-end sequencing
rsem-calculate-expression -num--threads 8 \
--star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--star-gzipped-read-file \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
/your_path/fqClean.fq.gz \
$index_hg38v43/index_star /your_path/quantRSEM
```

There are two key arguments you need to specify: 1) `**--strandedness [none|forward|reverse]**`, which is same with Bowtie2-RSEM; and 2) `**--star-gzipped-read-file**`, which should be enabled if the input FASTQ files are gzipped.

The outputs are exactly same with those from Bowtie2-RSEM.

### 3. STAR-HTSeq (optional)

In addtion RSEM, HTSeq (https://htseq.readthedocs.io/en/master/) is another popular quantifier for RNA-Seq data. Compared with RSEM which employed the Expectaton Maximization-algorithsm to handle the ambiguously-mapped reads, the way HTSeq counts reads is simpler and more transparent: as shown in the figure below, **HTSeq does NOT count ambiguously-mapped reads and it does NOT weight the reads by either their overlaps with refence genes or mapping scores.**

![_images/count_modes.png](https://htseq.readthedocs.io/en/master/_images/count_modes.png)

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/STAR/index_overhang100
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/STAR/index_overhang100

## Alignment by STAR: for paired-end sequencing
STAR --genomeDir $index_hg38oh100 --readFilesIn /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz --readFilesCommand zcat \
--runThreadN 8 --twopassMode Basic --quantMode TranscriptomeSAM \
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
--outSAMattributes All --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outFileNamePrefix /your_path/quantSTAR

## Alignment by STAR: for single-end sequencing
STAR --genomeDir $index_hg38oh100 --readFilesIn /your_path/fqClean.fq.gz --readFilesCommand zcat \
--runThreadN 8 --twopassMode Basic --quantMode TranscriptomeSAM \
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
--outSAMattributes All --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outFileNamePrefix /your_path/quantSTAR

## Quantification by HTSeq
htseq-count -f bam -r pos -s reverse -i gene_id -m intersection-nonempty -n 8 /your_path/quantSTAR/Aligned.out.bam $gtf_hg38 > /your_path/quantSTAR/htseq_counts.txt
```

**<u>Key arguments</u>**:

In this pipeline, STAR was employed as the aligner, with the arguments recommended by GDC (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-transformation). The only key argument for STAR you need to specify is the `**--outSAMattrRGline**`.

For HTSeq, there are two key arguments to be specified: 1) `**-r/--order [pos|name]**`, which is the way BAM files were sorted; and 2) `**-s/--stranded [no|yes|reverse]**`, which is the library type.

**<u>Key outputs:</u>**

* **htseq_counts.txt**: contains gene-level abundance estimates, raw counts only.
* **Aligned.sortedByCoord.out.bam**: BAM file with alignments to reference genome, sorted by coordinates.
* **Aligned.toTranscriptome.out.bam**: BAM file with alignments to reference transcriptome.
* **Log.final.out**: contains alignment statistics generated by STAR. This file provides the number of totally-mapped and uniquely-mapped reads used in the final QC report.