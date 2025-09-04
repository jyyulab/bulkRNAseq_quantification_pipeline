#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);

my $f_input = $ARGV[0];

print "Job Started:\n\tWe are loading $f_input...\n";
open (INPUT, $f_input) or die;
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $dir_preProcess = "$F[5]/$F[0]/preProcessing";
    make_path($dir_preProcess);

    open (OUT, "> $dir_preProcess/adapterTrimming.sh") or die;
    print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 4000\n#BSUB -oo $dir_preProcess/adapterTrimming.out -eo $dir_preProcess/adapterTrimming.err\n#BSUB -J adapterTrimming_$F[0]\n#BSUB -q superdome\n\n";

    if ($F[1] eq "PE") {
        if ((-e "$dir_preProcess/fqRaw_R1.fq.gz") and (-e "$dir_preProcess/fqRaw_R2.fq.gz")) {
            print OUT "fastp -w 8 -l 30 -q 20 -n 5 -h $dir_preProcess/adapterTrimming.html -j $dir_preProcess/adapterTrimming.json -i $dir_preProcess/fqRaw_R1.fq.gz -I $dir_preProcess/fqRaw_R2.fq.gz -o $dir_preProcess/fqClean_R1.fq.gz -O $dir_preProcess/fqClean_R2.fq.gz";
        }
        else {
            die "The inputs files, fqRaw_R1.fq.gz and fqRaw_R1.fq.gz, were not found in $dir_preProcess\n";
        }
    }
    elsif ($F[1] eq "SE") {
        if (-e "$dir_preProcess/fqRaw.fq.gz") {
            print OUT "fastp -w 8 -l 30 -q 20 -n 5 -h $dir_preProcess/adapterTrimming.html -j $dir_preProcess/adapterTrimming.json -i $dir_preProcess/fqRaw.fq.gz -o $dir_preProcess/fqClean.fq.gz";
        }
        else {
            die "The inputs file, fqRaw.fq.gz, was not found in $dir_preProcess\n";
        }
    }
    close OUT;

    `bsub < $dir_preProcess/adapterTrimming.sh`;
    print "Adapter trimming jobs for $F[0] have been submitted.\n";
}
close INPUT;





