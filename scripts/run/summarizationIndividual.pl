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

print "Job Started: $f_input is loading ...\n";
open (INPUT, $f_input) or die;
while (<INPUT>) {
    chomp;
    next if ($_ =~ /^sampleID/);
    next if ($_ =~ /^#/);

    my @F = split(/\t/, $_);
    my $dir_summarization = "$F[5]/$F[0]/summarization";
    make_path($dir_summarization);

    if ((-e "$F[5]/$F[0]/quantRSEM_STAR/quant.isoforms.results") and (-e "$F[5]/$F[0]/quantRSEM_STAR/quant.genes.results") and (-e "$F[5]/$F[0]/quantSalmon/quant.sf") and (-e "$F[5]/$F[0]/quantSalmon/quant.genes.sf") and (-e "$F[5]/$F[0]/preProcessing/adapterTrimming.json") and (-e "$F[5]/$F[0]/quantRSEM_STAR/quant.stat/quant.cnt") and (-e "$F[5]/$F[0]/quantRSEM_STAR/genebodyCoverage.txt")) {
        if (-e "$F[5]/$F[0]/summarization/summarization.html") {
            print "\tSummarization jobs have completed for $F[0].\n";
        }
        else {
            unless (-e "$F[5]/$F[0]/summarization/summarization.sh") {
                open (OUT, "> $dir_summarization/summarization.sh") or die;
                print OUT "#BSUB -P bulkRNAseqQuantification\n#BSUB -n 8\n#BSUB -M 8000\n#BSUB -oo $dir_summarization/summarization.out -eo $dir_summarization/summarization.err\n#BSUB -J summarization_$F[0]\n#BSUB -q superdome\n\n";
                print OUT "Rscript -e \"rmarkdown::render(input = '$dir_rmd/summarizationIndividual.Rmd', clean = TRUE, quiet = F, output_format = 'html_document', output_file = 'summarization.html', output_dir = '$F[5]/$F[0]/summarization', params = list(sampleName = '$F[0]', dir_quant = '$F[5]/$F[0]', dir_anno = '$F[3]'))\"\n";
                close OUT;
            }

            `bsub < $dir_summarization/summarization.sh`;
            print "\tSummarization jobs for $F[0] have been submitted.\n";
        }
    }
    else {
        print "\tSummarization jobs for $F[0] FAILED:\n";
        print "\t\t$F[5]/$F[0]/quantRSEM_STAR/quant.isoforms.results was not found.\n" unless (-e "$F[5]/$F[0]/quantRSEM_STAR/quant.isoforms.results");
        print "\t\t$F[5]/$F[0]/quantRSEM_STAR/quant.genes.results was not found.\n" unless (-e "$F[5]/$F[0]/quantRSEM_STAR/quant.genes.results");
        print "\t\t$F[5]/$F[0]/quantSalmon/quant.sf was not found.\n" unless (-e "$F[5]/$F[0]/quantSalmon/quant.sf");
        print "\t\t$F[5]/$F[0]/quantSalmon/quant.genes.sf was not found.\n" unless (-e "$F[5]/$F[0]/quantSalmon/quant.genes.sf");
        print "\t\t$F[5]/$F[0]/preProcessing/adapterTrimming.json was not found.\n" unless (-e "$F[5]/$F[0]/preProcessing/adapterTrimming.json");
        print "\t\t$F[5]/$F[0]/quantRSEM_STAR/quant.stat/quant.cnt was not found.\n" unless (-e "$F[5]/$F[0]/quantRSEM_STAR/quant.stat/quant.cnt");
        print "\t\t$F[5]/$F[0]/quantRSEM_STAR/genebodyCoverage.txt was not found.\n" unless (-e "$F[5]/$F[0]/quantRSEM_STAR/genebodyCoverage.txt");
        next;
    }
}
close INPUT;
