#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);

my $conda_prefix = $ENV{'CONDA_PREFIX'};
if (defined $conda_prefix) {
    print "The CONDA_PREFIX is: $conda_prefix\n";
} else {
    print "CONDA_PREFIX is not set.\n";
}

my $dir_bin = "$conda_prefix/bin";
my $f_input = $ARGV[0];

print "Job Started: $f_input is loading...\n";
open (INPUT, $f_input) or die;
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $dir_preProcess = "$F[5]/$F[0]/preProcessing";
    my $dir_quantSTAR_HTSeq = "$F[5]/$F[0]/quantSTAR_HTSeq";
    make_path($dir_quantSTAR_HTSeq);

    my $libtype = "";
    if (-e "$F[5]/$F[0]/quantSalmon/lib_format_counts.json") {
        open (LT, "$F[5]/$F[0]/quantSalmon/lib_format_counts.json") or die;
        while (<LT>) {
            chomp;
            next unless ($_ =~ /\"expected_format\": \"(\w+)\",/);
            if ($1 =~ /U$/) {
                $libtype = "no";
            }
            elsif ($1 =~ /R$/) {
                $libtype = "reverse";
            }
            elsif ($1 =~ /F$/) {
                $libtype = "yes";
            }
            else {
                print "The putative library type of $F[0], $1, was not recoganized.\n";
            }
        }
        close LT;
    }
    else {
        print "The lib_format_counts.json file of $F[0] was not found. Please check if the Salmon quantification was completely done.\n";
    }

    open (OUT, "> $dir_quantSTAR_HTSeq/quantSTAR_HTSeq.sh") or die;
    print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 8000\n#BSUB -oo $dir_quantSTAR_HTSeq/quantSTAR_HTSeq.out -eo $dir_quantSTAR_HTSeq/quantSTAR_HTSeq.err\n#BSUB -J quantSTAR_HTSeq_$F[0]\n#BSUB -q superdome\n\n";

    if ($F[1] eq "PE") {
        if ((-e "$dir_preProcess/fqClean_R1.fq.gz") and (-e "$dir_preProcess/fqClean_R2.fq.gz")) {
            print OUT "$dir_bin/STAR --readFilesIn $dir_preProcess/fqClean_R1.fq.gz $dir_preProcess/fqClean_R2.fq.gz --outFileNamePrefix $dir_quantSTAR_HTSeq/ --outSAMattrRGline ID:$F[0] SM:$F[0] LB:Illumina PL:Illumina PU:Illumina --genomeDir $F[3]/bulkRNAseq/STAR/index_overhang100 --readFilesCommand zcat --runThreadN 8 --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --genomeLoad NoSharedMemory --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimOutJunctionFormat 1 --chimMainSegmentMultNmax 1 --outSAMattributes NH HI AS nM NM ch\n\nhtseq-count -f bam -r pos -s reverse -i gene_id -m intersection-nonempty -n 8 $dir_quantSTAR_HTSeq/Aligned.sortedByCoord.out.bam $F[3]/annotation.gtf > $dir_quantSTAR_HTSeq/htseq_counts.genes.txt\n\nhtseq-count -f bam -r pos -s reverse -i transcript_id -m intersection-nonempty -n 8 $dir_quantSTAR_HTSeq/Aligned.sortedByCoord.out.bam $F[3]/annotation.gtf > $dir_quantSTAR_HTSeq/htseq_counts.transcripts.txt";
        }
        else {
            die "The inputs files, fqClean_R1.fq.gz and fqClean_R1.fq.gz, were not found in $dir_preProcess.\n";
        }
    }
    elsif ($F[1] eq "SE") {
        if (-e "$dir_preProcess/fqClean.fq.gz") {
            print OUT "$dir_bin/STAR --readFilesIn $dir_preProcess/fqClean.fq.gz --outFileNamePrefix $dir_quantSTAR_HTSeq/ --outSAMattrRGline ID:$F[0] SM:$F[0] LB:Illumina PL:Illumina PU:Illumina --genomeDir $F[3]/bulkRNAseq/STAR/index_overhang100 --readFilesCommand zcat --runThreadN 8 --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --genomeLoad NoSharedMemory --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimOutJunctionFormat 1 --chimMainSegmentMultNmax 1 --outSAMattributes NH HI AS nM NM ch\n\nhtseq-count -f bam -r pos -s reverse -i gene_id -m intersection-nonempty -n 8 $dir_quantSTAR_HTSeq/Aligned.sortedByCoord.out.bam $F[3]/annotation.gtf > $dir_quantSTAR_HTSeq/htseq_counts.genes.txt\n\nhtseq-count -f bam -r pos -s reverse -i transcript_id -m intersection-nonempty -n 8 $dir_quantSTAR_HTSeq/Aligned.sortedByCoord.out.bam $F[3]/annotation.gtf > $dir_quantSTAR_HTSeq/htseq_counts.transcripts.txt";
        }
        else {
            die "The input file, fqClean.fq.gz, was not found in $dir_preProcess.\n";
        }
    }
    close OUT;

    `bsub < $dir_quantSTAR_HTSeq/quantSTAR_HTSeq.sh`;
    print "\tSTAR quantification jobs for $F[0] have been submitted.\n";
}
close INPUT;
