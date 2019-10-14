#BSUB -P RNASeq_pipelines
#BSUB -n 8
#BSUB -M 5000
#BSUB -oo 04_quantHTSeq.out -eo 04_quantHTSeq.err
#BSUB -J 04_quantHTSeq
#BSUB -q standard

# software
star=/research/rgs01/applications/hpcf/apps/star/install/2.5.3a/bin/STAR
htseq=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env/bin/htseq-count
bedtools=/research/rgs01/applications/hpcf/apps/bedtools/install/2.25.0/bin/bedtools-bed
samtools=/research/rgs01/applications/hpcf/apps/samtools/install/1.2/bin/samtools
genebody=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/RNASeq_pipelines/05_genebodyCoverage.R

# database
index_hg38oh100=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/STAR/index_overhang100
index_hg38oh150=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/STAR/index_overhang150
index_hg19oh100=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/STAR/index_overhang100
index_hg19oh150=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/STAR/index_overhang150
index_mm10oh100=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/STAR/index_overhang100
index_mm10oh150=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/STAR/index_overhang150
binlist_hg38=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/binlist_150.txt
binlist_hg19=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/binlist_150.txt
binlist_mm10=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/binlist_150.txt
gtf_hg38=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/annotation.gtf
gtf_hg19=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/annotation.gtf
gtf_mm10=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/annotation.gtf

# I/O
indir=
outdir=

## For single-end sequencing
# STAR Alignment
$star --runMode alignReads --runThreadN 8 --twopassMode Basic \
    --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/sample/ --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outSAMattributes All --outSAMunmapped Within \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignMatesGapMax 1000000 \
    --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimSegmentMin 15 --limitSjdbInsertNsj 1200000 --chimOutType SeparateSAMold SoftClip \
    --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outSAMstrandField intronMotif \
    --readFilesCommand zcat --genomeDir $index_hg38oh100 --readFilesIn $indir/sample.clean.fq.gz
# HTSeq Quantification
$htseq -f bam -r pos -s reverse -a 10 -t exon -i gene_id -m intersection-nonempty --nonunique none --secondary-alignments score --supplementary-alignments score $outdir/sample/Aligned.out.bam $gtf_hg38 > $outdir/sample/htseq_counts.txt

# Gene Body Coverage
$samtools sort $outdir/sample/Aligned.toTranscriptome.out.bam $outdir/sample/Aligned.toTranscriptome.out.sorted && $samtools index $outdir/sample/Aligned.toTranscriptome.out.sorted.bam
$bedtools multicov -bams $outdir/sample/Aligned.toTranscriptome.out.sorted.bam -bed $binlist_hg38 > $outdir/sample/readsDistribution.txt
$genebody $outdir/sample/readsDistribution.txt $outdir/sample/genebodyCoverage

## For paired-end sequencing
# STAR Alignment
$star --runMode alignReads --runThreadN 8 --twopassMode Basic \
    --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/sample/ --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outSAMattributes All --outSAMunmapped Within \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignMatesGapMax 1000000 \
    --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimSegmentMin 15 --limitSjdbInsertNsj 1200000 --chimOutType SeparateSAMold SoftClip \
    --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outSAMstrandField intronMotif \
    --readFilesCommand zcat --genomeDir $index_hg38oh100 --readFilesIn $indir/sample_R1.clean.fq.gz $indir/sample_R2.clean.fq.gz
# HTSeq Quantification
$htseq-count -f bam -r pos -s reverse -a 10 -t exon -i gene_id -m intersection-nonempty --nonunique none --secondary-alignments score --supplementary-alignments score $outdir/sample/Aligned.out.bam $gtf_hg38 > $outdir/sample/htseq_counts.txt

# Gene Body Coverage
$samtools sort $outdir/sample/Aligned.toTranscriptome.out.bam $outdir/sample/Aligned.toTranscriptome.out.sorted && $samtools index $outdir/sample/Aligned.toTranscriptome.out.sorted.bam
$bedtools multicov -bams $outdir/sample/Aligned.toTranscriptome.out.sorted.bam -bed $binlist_hg38 > $outdir/sample/readsDistribution.txt
$genebody $outdir/sample/readsDistribution.txt $outdir/sample/genebodyCoverage

# Options to set
# 1) --quantMode TranscriptomeSAM: use it to generate the alignment to transcriptome. Highly recommended.
# 2) --outSAMattrRGline: edit it everytime.
# 3) --readFilesCommand: zcat for .gz files, bzcat for .bz2 files.
# 4) Here we follow the GDC mRNA quantification analysis pipeline: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/.
# 5) For more information about STAR, please check out here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf.
# 6) For more information about HTSeq, please check out here: https://htseq.readthedocs.io/en/release_0.11.1/count.html.

# Key Outputs:
# 1) Aligned.sortedByCoord.out.bam: Alignments to reference genome, sorted by coordinates.
# 2) Aligned.toTranscriptome.out.bam: Alignments to reference transcriptome.
# 3) htseq_counts.txt: HTSeq gene expression quantification, COUNTS.


