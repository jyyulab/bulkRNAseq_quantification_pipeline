#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV == 1) {
    unless (-e $ARGV[0]) {
        print "ERROR: The INPUT file doesn't exist. Please check and retry.\n"; die;
    }
} elsif (@ARGV < 1) {
    print "ERROR: Please specify the gene annotation file to parse!\n\n\tperl parseAnnotation.pl xxx.annotation.gtf\n\n"; die;
} else {
    print "ERROR: Please use only ONE argument which specifies the gene annotation file to parse!\n\n\tperl parseAnnotation.pl xxx.annotation.gtf\n\n"; die;
}

my $f_input = $ARGV[0];
if ($f_input =~ /\.gtf$/ig) {
    print "We are parsing $f_input ...\n";

    open (IN, $f_input) or die;
    my %tr = (); my %gene = (); my %tr2gene = (); my %gene2tr = (); my %hash = ();
    while (<IN>) {
        chomp;
        next if ($_ =~ /^#/);
        my @F = split(/\t/, $_);
        next unless ($F[2] eq "transcript");
        $F[8] =~ /gene_id "(ENS.+)"; transcript_id "(ENS.+)"; gene_type "(\w+)";.* gene_name "(.+)"; transcript_type "(\w+)";.* transcript_name "(.+)"; level /;
        $tr{$2} = "$2\t$6\t$5\t$1\t$4";
        $gene{$1} = "$1\t$4\t$3";
        $tr2gene{$2} = "$2\t$1"; $gene2tr{$2} = "$1\t$2"; $hash{$2} = "$2\t$6\t$1\t$4";
    }
    close IN;

    print "We are writing the outputs ...\n";
    (my $f_output = $f_input) =~ s/\.gtf$//g;
    open (OUT1, "> $f_output.transcriptAnnotation.txt") or die;
    print OUT1 "transcript_id\ttranscript_name\ttranscript_type\tgene_id\tgene_name\n";
    for my $x (sort keys %tr) {
        print OUT1 "$tr{$x}\n";
    }
    close OUT1;

    open (OUT2, "> $f_output.geneAnnotation.txt") or die;
    print OUT2 "gene_id\tgene_name\tgene_type\n";
    for my $x (sort keys %gene) {
        print OUT2 "$gene{$x}\n";
    }
    close OUT2;

    open (OUT3, "> $f_output.transcript2gene.txt") or die;
    open (OUT4, "> $f_output.gene2transcript.txt") or die;
    for my $x (sort keys %hash) {
        print OUT3 "$tr2gene{$x}\n";
        print OUT4 "$gene2tr{$x}\n";
    }
    close OUT3; close OUT4;
    print "Done!\n";
}
else {
    print "ERROR: The INPUT file must be in GTF(Gene Transfer Format) format, and the file name much be ended with '.gtf' or '.GTF'.\n";
}

