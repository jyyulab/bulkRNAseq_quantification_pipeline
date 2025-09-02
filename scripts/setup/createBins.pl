#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);

## We take the longest transcript to present each gene. The genes of which the longest transcript is less than 300nt were excluded.

if (@ARGV == 3) {
    if (-e $ARGV[0]) { unless ($ARGV[0] =~ /\.fa$/ig) { print "ERROR: The transcriptome sequence file, $ARGV[0], must be in FASTA format, and the file name much be ended with '.fa'.\n"; die; }
    } else { print "ERROR: The transcriptome sequence file, $ARGV[0], doesn't exist. Please check and retry.\n"; die; }
    unless (-e $ARGV[1]) { print "ERROR: The housekeeping gene file, $ARGV[1], doesn't exist. Please check and retry.\n"; die; }
    unless (-d $ARGV[2]) { print "The output directory, $ARGV[2], doesn't exist. We are generating it...\n"; make_path($ARGV[2]); }
} else {
    print "ERROR: Please use THREE arguments to specify 1) transcriptome sequence file, 2) housekeeping gene file, 3) directory to save output files, respectively.\n\n\tperl prepareBins.pl xxx.transcripts.fa housekeeping_genes.xxx.txt ./bulkRNAseq/genebodyBins\n\n"; die;
}

my $input_gtf = $ARGV[0];
my $input_hk = $ARGV[1];
my $outout_dir = $ARGV[2];

open (REF, $input_gtf) or die;
print "We are paring $input_gtf ...\n";
my %ref = ();
while (<REF>) {
    chomp;
    next unless ($_ =~ /^>ENS/);
    my @F = split(/\|/, $_);
    next if ($F[6] < 300);
    $F[0] =~ /^>(.+)/; my $isoform = $1;
    if (exists $ref{$F[5]}) {
        for my $x (sort keys %{$ref{$F[5]}}) {
            if ($ref{$F[5]}{$x} > $F[6]) {
                next;
            }
            else {
                delete $ref{$F[5]}{$x};
                $ref{$F[5]}{$isoform} = $F[6];
            }
        }
    }
    else {
        $ref{$F[5]}{$isoform} = $F[6];
    }
}
close REF;

my %hk_genes = ();
open (HK, $input_hk) or die;
print "We are collecting housekeeping genes from $input_hk ...\n";
while (<HK>) {
    chomp;
    next if ($_ =~ /^#/);
    my @HK = split(/\t/, $_);
    $hk_genes{$HK[0]}++;
}
close HK;

open (OUT1, "> $outout_dir/genebodyBins_all.txt") or die;
open (OUT2, "> $outout_dir/genebodyBins_housekeeping.txt") or die;
print "We are writing the outputs ...\n";
for my $x (sort keys %ref) {
    for my $y (sort keys %{$ref{$x}}) {
        my $length = $ref{$x}{$y};
        my $width = $length / 100;

        for my $i (0 .. 99) {
            my $start = int($width * $i) + 1;
            my $end = int($width * ($i + 1));
            $end = $length if ($end > $length);
            my $k = $i + 1;
            print OUT1 "$y\t$start\t$end\t$x\|Bin$k\n";
            print OUT2 "$y\t$start\t$end\t$x\|Bin$k\n" if (exists($hk_genes{$x}));
        }
    }
}
close OUT1; close OUT2;
print "Done!\n";
