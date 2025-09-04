#!/usr/bin/perl

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
    #mkdir $dir_preProcess unless (-e $dir_preProcess);

    open (OUT, "> $dir_preProcess/all2Fastq.sh") or die;
    print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 4000\n#BSUB -oo $dir_preProcess/all2Fastq.out -eo $dir_preProcess/all2Fastq.err\n#BSUB -J all2Fastq_$F[0]\n#BSUB -q standard\n\n";

    if ($F[4] =~ /(bam)$/i or $F[4] =~ /(sam)$/i) {
        my $format = lc($1);
        print OUT "bamtofastq filename=$F[4] inputformat=$format gz=1 F=$dir_preProcess/fqRaw_R1.fq.gz F2=$dir_preProcess/fqRaw_R2.fq.gz\n" if ($F[1] eq "PE");
        print OUT "bamtofastq filename=$F[4] inputformat=$format gz=1 S=$dir_preProcess/fqRaw.fq.gz\n" if ($F[1] eq "SE");
    }
    elsif ($F[4] =~ /fq$|fastq$|fq.gz$|fastq.gz$/i) {
        if ($F[4] =~ /;/) {
            my @mates = split(/;/, $F[4]);
            $mates[0] =~ s/\s+//g; $mates[0] =~ s/,/ /g;
            $mates[1] =~ s/\s+//g; $mates[1] =~ s/,/ /g;
            print OUT "cat $mates[0] > $dir_preProcess/fqRaw_R1.fq.gz\ncat $mates[1] > $dir_preProcess/fqRaw_R2.fq.gz\n";
        }
        else {
            if ($F[4] =~ /,/) {
                $F[4] =~ s/\s+//g; $F[4] =~ s/,/ /g;
            }
            print OUT "cat $F[4] > $dir_preProcess/fqRaw.fq.gz\n";
        }
    }
    else {
        print "The format of $F[4] can not be recoganized.\n";
    }
    close OUT;

    `bsub < $dir_preProcess/all2Fastq.sh`;
    print "\tThe all2Fastq job for $F[0] has been sumitted.\n";
}
close INPUT;
