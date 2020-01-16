#!/usr/bin/env Rscript

####################
## Usage: Rscript 05_genebodyCoverage.R readsDistribution.txt prefix_of_output
##
## Inputs: Two Arguments
##      1) readsDistribution.txt: this file is generated in previous step. It's the reads counts mapped to each bin for each TRANSCRIPT.
##      2) prefix_of_output: The is the prefix of the outputs, including prefix_of_output.txt and prefix_of_output.pdf
##
## Outputs:
##      1) prefix_of_output.txt: it's the reads counts mapped to each bin for each GENE.
##      2) prefix_of_output.pdf: visulization of reads distribution along GENEs.
####################


args = commandArgs(trailingOnly=TRUE)
paste0(args[1], " is in process...")

input <- read.table(args[1], header = F, sep = "\t")
sel <- input[,4:5];
sel$V4 <- gsub(".+Bin(\\d+)", "\\1", sel$V4)
x.sel <- aggregate(. ~ V4, sel, sum, na.rm = T)
colnames(x.sel) <- c("GeneBodyBins", "NumberOfReads")

output1 <- paste0(args[2], ".txt")
write.table(x.sel, file = output1, quote = F, row.names = F, col.names = T, sep = "\t")

output2 <- paste0(args[2], ".pdf")
pdf(file = output2, width = 10, height = 10)
plot(x.sel, type = "h", pch = 19, col = "red", xlab = "5'-Gene Body-3'", ylab = "Number of Reads", lwd = 2, cex.lab = 1.5)
dev.off()
