#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);

my $conda_prefix = $ENV{'CONDA_PREFIX'};
if (defined $conda_prefix) {
    print "The CONDA_PREFIX is: $conda_prefix\n";
} else {
    print "CONDA_PREFIX is not set.\n";
}

my $dir_rmd = "$conda_prefix/pipeline/scripts/run";
my $f_input = $ARGV[0];
my $dir_output = $ARGV[1];

print "Job Started: we are parsing $f_input ...\n";
open (INPUT, $f_input) or die;
my %ref_sample = ();
my %ref_path = ();
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $str_reference = $F[3];
    $str_reference =~ s/\//_/g;
    unless (exists $ref_sample{$str_reference}) { make_path("$dir_output/$str_reference"); }
    unless (exists $ref_path{$str_reference}) { $ref_path{$str_reference} = $F[3]; }

    my $dir_summarization = "$F[5]/$F[0]/summarization";
    if ((-e "$dir_summarization/quant.transcripts.txt") and (-e "$dir_summarization/quant.genes.txt") and (-e "$dir_summarization/01_alignmentStatistics.txt") and (-e "$dir_summarization/02_quantificationStatistics.txt") and (-e "$dir_summarization/03_biotypeDistribution.transcripts.txt") and (-e "$dir_summarization/03_biotypeDistribution.genes.txt") and (-e "$dir_summarization/05_genebodyCoverage.txt")) {
        push @{ $ref_sample{$str_reference} }, $_;
        next;
    } else {
        print "\t\t$dir_summarization/quant.transcripts.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/quant.transcripts.txt");
        print "\t\t$dir_summarization/quant.genes.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/quant.genes.txt");
        print "\t\t$dir_summarization/01_alignmentStatistics.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/01_alignmentStatistics.txt");
        print "\t\t$dir_summarization/02_quantificationStatistics.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/02_quantificationStatistics.txt");
        print "\t\t$dir_summarization/03_biotypeDistribution.transcripts.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/03_biotypeDistribution.transcripts.txt");
        print "\t\t$dir_summarization/03_biotypeDistribution.genes.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/03_biotypeDistribution.genes.txt");
        print "\t\t$dir_summarization/05_genebodyCoverage.txt is not found. Sample: $F[0] skipped...\n" unless (-e "$dir_summarization/05_genebodyCoverage.txt");
        next;
    }
}
close INPUT;

my $ref_count = keys %ref_sample;
if ($ref_count == 1) {print "1 reference genome is found in $f_input...\n";} else {print "$ref_count refrence genomes are found in $f_input...\n";}

for my $x (sort keys %ref_sample) {
    open (OUT1, "> $dir_output/$x/sampleTable.txt") or die;
    print OUT1 "sampleID\tlibraryType\tphredMethod\treference\tinput\toutput\n";
    my $y = join("\n", @{$ref_sample{$x}});
    print OUT1 "$y\n";
    close OUT1;

    open (OUT2, "> $dir_output/$x/summarizationMultiple.sh") or die;
    print OUT2 "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 8000\n#BSUB -oo $dir_output/$x/summarizationMultiple.out -eo $dir_output/$x/summarizationMultiple.err\n#BSUB -J summarizationMultiple\n#BSUB -q superdome\n\n";
    print OUT2 "Rscript -e \"rmarkdown::render(input = '$dir_rmd/summarizationMultiple.Rmd', clean = TRUE, quiet = FALSE, output_format = 'html_document', output_file = 'summarizationMultiple.html', output_dir = '$dir_output/$x', params = list(sampleTable = '$dir_output/$x/sampleTable.txt', dir_anno = '$ref_path{$x}', dir_output = '$dir_output/$x'))\"\n";
    close OUT2;

    `bsub < $dir_output/$x/summarizationMultiple.sh`;
    print "\t$dir_output/$x/summarizationMultiple.sh has been submitted.\n";
}
