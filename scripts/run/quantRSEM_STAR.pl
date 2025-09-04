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
    my $dir_quantRSEM = "$F[5]/$F[0]/quantRSEM_STAR";
    make_path($dir_quantRSEM);

    my $libtype = "";
    if (-e "$F[5]/$F[0]/quantSalmon/lib_format_counts.json") {
        open (LT, "$F[5]/$F[0]/quantSalmon/lib_format_counts.json") or die;
        while (<LT>) {
            chomp;
            next unless ($_ =~ /\"expected_format\": \"(\w+)\",/);
            if ($1 =~ /U$/) {
                $libtype = "none";
            }
            elsif ($1 =~ /R$/) {
                $libtype = "reverse";
            }
            elsif ($1 =~ /F$/) {
                $libtype = "forward";
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

    open (OUT, "> $dir_quantRSEM/quantRSEM_STAR.sh") or die;
    print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 4000\n#BSUB -oo $dir_quantRSEM/quantRSEM_STAR.out -eo $dir_quantRSEM/quantRSEM_STAR.err\n#BSUB -J quantRSEM_STAR_$F[0]\n#BSUB -q superdome\n\n";

    if ($F[1] eq "PE") {
        if ((-e "$dir_preProcess/fqClean_R1.fq.gz") and (-e "$dir_preProcess/fqClean_R2.fq.gz")) {
            print OUT "$dir_bin/rsem-calculate-expression --num-threads 8 --star --star-path $dir_bin --star-gzipped-read-file --strandedness $libtype --phred33-quals --sort-bam-by-coordinate --paired-end $dir_preProcess/fqClean_R1.fq.gz $dir_preProcess/fqClean_R2.fq.gz $F[3]/bulkRNAseq/RSEM/index_star/assembly $dir_quantRSEM/quant";
        }
        else {
            die "The inputs files, fqClean_R1.fq.gz and fqClean_R1.fq.gz, were not found in $dir_preProcess.\n";
        }
    }
    elsif ($F[1] eq "SE") {
        if (-e "$dir_preProcess/fqClean.fq.gz") {
            print OUT "$dir_bin/rsem-calculate-expression --num-threads 8 --star --star-path $dir_bin --star-gzipped-read-file --strandedness $libtype --phred33-quals --sort-bam-by-coordinate $dir_preProcess/fqClean.fq.gz $F[3]/bulkRNAseq/RSEM/index_star/assembly $dir_quantRSEM/quant";
        }
        else {
            die "The input file, fqClean.fq.gz, was not found in $dir_preProcess.\n";
        }
    }
    close OUT;

    `bsub < $dir_quantRSEM/quantRSEM_STAR.sh`;
    print "\tRSEM quantification jobs for $F[0] have been submitted.\n";
}
close INPUT;
