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

print "Job Start: we are loading $f_input...\n";
open (INPUT, $f_input) or die;
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $dir_preProcess = "$F[5]/$F[0]/preProcessing";
    my $dir_quantSalmon = "$F[5]/$F[0]/quantSalmon";
    make_path($dir_quantSalmon);

    open (OUT, "> $dir_quantSalmon/quantSalmon.sh") or die;
    print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 4000\n#BSUB -oo $dir_quantSalmon/quantSalmon.out -eo $dir_quantSalmon/quantSalmon.err\n#BSUB -J quantSalmon_$F[0]\n#BSUB -q superdome\n\n";

    if ($F[1] eq "PE") {
        if ((-e "$dir_preProcess/fqClean_R1.fq.gz") and (-e "$dir_preProcess/fqClean_R2.fq.gz")) {
            print OUT "$dir_bin/salmon quant -i $F[3]/bulkRNAseq/Salmon/index_decoy -l A -p 8 -g $F[3]/annotation.gtf -1 $dir_preProcess/fqClean_R1.fq.gz -2 $dir_preProcess/fqClean_R2.fq.gz --validateMappings -o $dir_quantSalmon";
        }
        else {
            die "The inputs files, fqClean_R1.fq.gz and fqClean_R1.fq.gz, were not found in $dir_preProcess\n";
        }
    }
    elsif ($F[1] eq "SE") {
        if (-e "$dir_preProcess/fqClean.fq.gz") {
            print OUT "$dir_bin/salmon quant -i $F[3]/bulkRNAseq/Salmon/index_decoy -l A -p 8 -g $F[3]/annotation.gtf -r $dir_preProcess/fqClean.fq.gz --validateMappings -o $dir_quantSalmon";
        }
        else {
            die "The input file, fqClean.fq.gz, was not found in $dir_preProcess\n";
        }
    }
    close OUT;

    `bsub < $dir_quantSalmon/quantSalmon.sh`;
    print "\tSalmon quantification jobs for $F[0] have been submitted.\n";
}
close INPUT;

