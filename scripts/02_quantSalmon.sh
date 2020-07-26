#BSUB -P RNASeq_pipelines
#BSUB -n 8
#BSUB -M 2000
#BSUB -oo 02_quantSalmon.out -eo 02_quantSalmon.err
#BSUB -J 02_quantSalmon
#BSUB -q standard

# software
salmon=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/Salmon/bin/salmon

# database
index_hg38=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/Salmon/index_quasi
index_hg19=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/Salmon/index_quasi
index_mm10=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/Salmon/index_quasi
tr2gene_hg38=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/idmap.tr2gene.txt
tr2gene_hg19=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release32/idmap.tr2gene.txt
tr2gene_mm10=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM23/idmap.tr2gene.txt

# I/O
indir=
outdir=

# For single-end sequencing
$salmon quant -i $index -l A -p 8 -g $tr2gene -r $indir/sample.clean.fq.gz -o $outdir/sample

# For paired-end sequencing
$salmon quant -i $index -l A -p 8 -g $tr2gene -1 $indir/sample_R1.clean.fq.gz -2 $indir/sample_R2.clean.fq.gz -o $outdir/sample

# Key Outputs:
# 1) quant.sf: quantification results of transcripts/isoforms. Both TPM and COUNTS are provided.
# 2) quant.genes.sf: quantification results of genes. Both TPM and COUNTS are provided.
# 3) lib_format_counts.json: library type is predicted: IOM(inward, outward, matching) + SU(stranded, unstranded) + FR(Forward, Reverse)
# 4) For more information, please check out here: https://salmon.readthedocs.io/en/latest/salmon.html.

# NOTE
# 1) The default index files were generated under k=31. This is recommanded by Salmon, and works well for reads of 75bp or longer. If your reads are shorter than 31bp, no hits will be matched for them. In this case, you have to generate the index files by your self. The script to generate the index files is here: /research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release32/Salmon/00_buildIndex_quasi.sh.



