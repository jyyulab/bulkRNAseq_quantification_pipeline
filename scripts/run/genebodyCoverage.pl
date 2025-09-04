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

print "Job Started: $f_input is loading ...\n";
open (INPUT, $f_input) or die;
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $dir = "$F[5]/$F[0]/quantRSEM_STAR";
    my $bam = "$F[5]/$F[0]/quantRSEM_STAR/quant.transcript.sorted.bam";

    if (-e $bam) {
        open (OUT, "> $dir/genebodyCoverage.sh") or die;
        print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 4000\n#BSUB -oo $dir/genebodyCoverage.out -eo $dir/genebodyCoverage.err\n#BSUB -J genebodyCoverage_$F[0]\n#BSUB -q superdome\n\n";

        print OUT "$dir_bin/bedtools multicov -bams $bam -bed $F[3]/bulkRNAseq/genebodyBins/genebodyBins_housekeeping.txt > $dir/genebodyCoverage.txt";
        close OUT;

        `bsub < $dir/genebodyCoverage.sh`;
        print "Genebody coverage jobs for $F[0] have been submitted.\n";
    }
    else {
        print "The input file, $bam, was not found. Please check and retry.\n";
        next;
    }
}
close INPUT;
