#!/usr/bin/env Rscript

## 0 configure the inputs
indir_salmon <- "/research/rgs01/project_space/yu3grp/brainTumor_JY/yu3grp/m6A/m6A/04_quantification/01_salmon" # dir of Salmon Outputs
indir_rsem <- "/research/rgs01/project_space/yu3grp/brainTumor_JY/yu3grp/m6A/m6A/04_quantification/02_rsem" # dir of RSEM Outputs
outdir <- "/research/rgs01/project_space/yu3grp/brainTumor_JY/yu3grp/m6A/m6A/04_quantification" # dir for Summary, including Correlation and Master Tables
samples <- c("sample1", "sample2")
cat("The input information has been read!\n")

## 1 collect expression information
# create data frames
masterTable_salmon_g <- data.frame(); masterTable_salmon_t <- data.frame();
masterTable_rsem_g <- data.frame(); masterTable_rsem_t <- data.frame();
# extract quantification information
cat("The quantification results are being collected...\n")
for (i in 1:length(samples)) {
    input_salmon_g <- read.table(paste0(indir_salmon, "/", samples[i], "/quant.genes.sf"), header = T, sep = "\t")
    input_salmon_t <- read.table(paste0(indir_salmon, "/", samples[i], "/quant.sf"), header = T, sep = "\t")
    input_rsem_g <- read.table(paste0(indir_rsem, "/", samples[i], "/quant.genes.results"), header = T, sep = "\t")
    input_rsem_t <- read.table(paste0(indir_rsem, "/", samples[i], "/quant.isoforms.results"), header = T, sep = "\t")

    if (i == 1) {
        masterTable_salmon_g <- data.frame(input_salmon_g); masterTable_salmon_t <- data.frame(input_salmon_t);
        masterTable_rsem_g <- data.frame(input_rsem_g); masterTable_rsem_t <- data.frame(input_rsem_t);
    } else {
        masterTable_salmon_g <- merge(masterTable_salmon_g, input_salmon_g, by = "Name", all = T)
        masterTable_salmon_t <- merge(masterTable_salmon_t, input_salmon_t, by = "Name", all = T)
        masterTable_rsem_g <- merge(masterTable_rsem_g, input_rsem_g, by = "gene_id", all = T)
        masterTable_rsem_t <- merge(masterTable_rsem_t, input_rsem_t, by = "transcript_id", all = T)
    }
}
# adjust the format
row.names(masterTable_salmon_g) <- masterTable_salmon_g$Name; masterTable_salmon_g <- masterTable_salmon_g[, -1];
row.names(masterTable_salmon_t) <- masterTable_salmon_t$Name; masterTable_salmon_t <- masterTable_salmon_t[, -1];
row.names(masterTable_rsem_g) <- masterTable_rsem_g$gene_id; masterTable_rsem_g <- masterTable_rsem_g[, -1];
row.names(masterTable_rsem_t) <- masterTable_rsem_t$transcript_id; masterTable_rsem_t <- masterTable_rsem_t[, -1];
cat("The quantification information have been collected!\n")

## 2 prepare and print expression matrix
# salmon_gene_TPM
output_salmon_g.tpm <- masterTable_salmon_g[, grepl("TPM", names(masterTable_salmon_g))]
names(output_salmon_g.tpm) <- samples; output_salmon_g.tpm <- output_salmon_g.tpm[rowSums(output_salmon_g.tpm) > 0,]
output_salmon_g.tpm.print <- data.frame(row.names = row.names(output_salmon_g.tpm), Gene_ID = row.names(output_salmon_g.tpm), output_salmon_g.tpm)
write.table(output_salmon_g.tpm.print, file = paste0(outdir, "/expressionMatrix_Salmon_geneLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# salmon_gene_count
output_salmon_g.count <- masterTable_salmon_g[, grepl("NumReads", names(masterTable_salmon_g))]
names(output_salmon_g.count) <- samples; output_salmon_g.count <- output_salmon_g.count[rowSums(output_salmon_g.count) > 0,]
output_salmon_g.count.print <- data.frame(row.names = row.names(output_salmon_g.count), Gene_ID = row.names(output_salmon_g.count), output_salmon_g.count)
write.table(output_salmon_g.count.print, file = paste0(outdir, "/expressionMatrix_Salmon_geneLevel_Count.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# salmon_isoform_TPM
output_salmon_t.tpm <- masterTable_salmon_t[, grepl("TPM", names(masterTable_salmon_t))]
names(output_salmon_t.tpm) <- samples; output_salmon_t.tpm <- output_salmon_t.tpm[rowSums(output_salmon_t.tpm) > 0,]
output_salmon_t.tpm.print <- data.frame(row.names = row.names(output_salmon_t.tpm), Gene_ID = row.names(output_salmon_t.tpm), output_salmon_t.tpm)
write.table(output_salmon_t.tpm.print, file = paste0(outdir, "/expressionMatrix_Salmon_isoformLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# salmon_isoform_count
output_salmon_t.count <- masterTable_salmon_t[, grepl("NumReads", names(masterTable_salmon_t))]
names(output_salmon_t.count) <- samples; output_salmon_t.count <- output_salmon_t.count[rowSums(output_salmon_t.count) > 0,]
output_salmon_t.count.print <- data.frame(row.names = row.names(output_salmon_t.count), Gene_ID = row.names(output_salmon_t.count), output_salmon_t.count)
write.table(output_salmon_t.count.print, file = paste0(outdir, "/expressionMatrix_Salmon_isoformLevel_Count.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_gene_FPKM
output_rsem_g.fpkm <- masterTable_rsem_g[, grepl("FPKM", names(masterTable_rsem_g))]
names(output_rsem_g.fpkm) <- samples; output_rsem_g.fpkm <- output_rsem_g.fpkm[rowSums(output_rsem_g.fpkm) > 0,]
output_rsem_g.fpkm.print <- data.frame(row.names = row.names(output_rsem_g.fpkm), Gene_ID = row.names(output_rsem_g.fpkm), output_rsem_g.fpkm)
write.table(output_rsem_g.fpkm.print, file = paste0(outdir, "/expressionMatrix_RSEM_geneLevel_FPKM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_gene_TPM
output_rsem_g.tpm <- masterTable_rsem_g[, grepl("TPM", names(masterTable_rsem_g))]
names(output_rsem_g.tpm) <- samples; output_rsem_g.tpm <- output_rsem_g.tpm[rowSums(output_rsem_g.tpm) > 0,]
output_rsem_g.tpm.print <- data.frame(row.names = row.names(output_rsem_g.tpm), Gene_ID = row.names(output_rsem_g.tpm), output_rsem_g.tpm)
write.table(output_rsem_g.tpm.print, file = paste0(outdir, "/expressionMatrix_RSEM_geneLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_gene_count
output_rsem_g.count <- masterTable_rsem_g[, grepl("expected_count", names(masterTable_rsem_g))]
names(output_rsem_g.count) <- samples; output_rsem_g.count <- output_rsem_g.count[rowSums(output_rsem_g.count) > 0,]
output_rsem_g.count.print <- data.frame(row.names = row.names(output_rsem_g.count), Gene_ID = row.names(output_rsem_g.count), output_rsem_g.count)
write.table(output_rsem_g.count.print, file = paste0(outdir, "/expressionMatrix_RSEM_geneLevel_Count.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_isoform_FPKM
output_rsem_t.fpkm <- masterTable_rsem_t[, grepl("FPKM", names(masterTable_rsem_t))]
names(output_rsem_t.fpkm) <- samples; output_rsem_t.fpkm <- output_rsem_t.fpkm[rowSums(output_rsem_t.fpkm) > 0,]
output_rsem_t.fpkm.print <- data.frame(row.names = row.names(output_rsem_t.fpkm), Gene_ID = row.names(output_rsem_t.fpkm), output_rsem_t.fpkm)
write.table(output_rsem_t.fpkm.print, file = paste0(outdir, "/expressionMatrix_RSEM_isoformLevel_FPKM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_isoform_TPM
output_rsem_t.tpm <- masterTable_rsem_t[, grepl("TPM", names(masterTable_rsem_t))]
names(output_rsem_t.tpm) <- samples; output_rsem_t.tpm <- output_rsem_t.tpm[rowSums(output_rsem_t.tpm) > 0,]
output_rsem_t.tpm.print <- data.frame(row.names = row.names(output_rsem_t.tpm), Gene_ID = row.names(output_rsem_t.tpm), output_rsem_t.tpm)
write.table(output_rsem_t.tpm.print, file = paste0(outdir, "/expressionMatrix_RSEM_isoformLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# rsem_isoform_count
output_rsem_t.count <- masterTable_rsem_t[, grepl("expected_count", names(masterTable_rsem_t))]
names(output_rsem_t.count) <- samples; output_rsem_t.count <- output_rsem_t.count[rowSums(output_rsem_t.count) > 0,]
output_rsem_t.count.print <- data.frame(row.names = row.names(output_rsem_t.count), Gene_ID = row.names(output_rsem_t.count), output_rsem_t.count)
write.table(output_rsem_t.count.print, file = paste0(outdir, "/expressionMatrix_RSEM_isoformLevel_Count.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
cat("The expression master tables have been created and printed!\n")

## 3 comparison of quant methods
# TPM_geneLevel
col_1 <- ""; col_2 <- ""; col_3 <- ""; col_4 <- ""; col_5 <- ""; col_6 <- ""; col_7 <- ""; col_8 <- "";
for (i in 1:length(samples)) {
    df_salmon <- data.frame(row.names = row.names(output_salmon_g.tpm), Salmon = output_salmon_g.tpm[,i])
    df_rsem <- data.frame(row.names = row.names(output_rsem_g.tpm), RSEM = output_rsem_g.tpm[,i])
    ident.salmon <- row.names(df_salmon)[df_salmon$Salmon > 0]
    ident.rsem <- row.names(df_rsem)[df_rsem$RSEM > 0]
    ident.overlap <- intersect(ident.salmon, ident.rsem)

    df <- data.frame(row.names = ident.overlap, IDs = ident.overlap)
    df <- merge(df, df_salmon, by = "row.names", all.x = T); row.names(df) <- df$Row.names; df <- df[, -1]
    df <- merge(df, df_rsem, by = "row.names", all.x = T); row.names(df) <- df$Row.names; df <- df[, -1]

    cor.pearson <- cor.test(df$Salmon, df$RSEM, method = "pearson", exact = T)
    cor.spearman <- cor.test(df$Salmon, df$RSEM, method = "spearman", exact = T)

    col_1[i] <- samples[i];
    col_2[i] <- length(ident.salmon); col_3[i] <- length(ident.rsem); col_4[i] <- length(ident.overlap)
    col_5[i] <- as.numeric(cor.pearson$estimate); col_6[i] <- as.numeric(cor.pearson$p.value)
    col_7[i] <- as.numeric(cor.spearman$estimate); col_8[i] <- as.numeric(cor.spearman$p.value)
}

comparison.table <- matrix(c("SampleName", col_1, "# IdentifiedBySalmon", col_2, "# IdentifiedByRSEM", col_3, "# IdentifiedByBoth", col_4,
                             "coef_Pearson", col_5, "pvalue_Pearson", col_6, "rho_Spearman", col_7, "pvalue_Spearman", col_8), ncol = 8, byrow = F)
write.table(comparison.table, file = paste0(outdir, "/03_comparisonOFmethods_geneLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
cat("The gene-level comparison of quantification methods is done!\n")

# TPM_isoformLevel
col_1 <- ""; col_2 <- ""; col_3 <- ""; col_4 <- ""; col_5 <- ""; col_6 <- ""; col_7 <- ""; col_8 <- "";
for (i in 1:length(samples)) {
    df_salmon <- data.frame(row.names = row.names(output_salmon_t.tpm), Salmon = output_salmon_t.tpm[,i])
    df_rsem <- data.frame(row.names = row.names(output_rsem_t.tpm), RSEM = output_rsem_t.tpm[,i])
    ident.salmon <- row.names(df_salmon)[df_salmon$Salmon > 0]
    ident.rsem <- row.names(df_rsem)[df_rsem$RSEM > 0]
    ident.overlap <- intersect(ident.salmon, ident.rsem)

    df <- data.frame(row.names = ident.overlap, IDs = ident.overlap)
    df <- merge(df, df_salmon, by = "row.names", all.x = T); row.names(df) <- df$Row.names; df <- df[, -1]
    df <- merge(df, df_rsem, by = "row.names", all.x = T); row.names(df) <- df$Row.names; df <- df[, -1]

    cor.pearson <- cor.test(df$Salmon, df$RSEM, method = "pearson", exact = T)
    cor.spearman <- cor.test(df$Salmon, df$RSEM, method = "spearman", exact = T)

    col_1[i] <- samples[i];
    col_2[i] <- length(ident.salmon); col_3[i] <- length(ident.rsem); col_4[i] <- length(ident.overlap)
    col_5[i] <- as.numeric(cor.pearson$estimate); col_6[i] <- as.numeric(cor.pearson$p.value)
    col_7[i] <- as.numeric(cor.spearman$estimate); col_8[i] <- as.numeric(cor.spearman$p.value)
}

comparison.table <- matrix(c("SampleName", col_1, "# IdentifiedBySalmon", col_2, "# IdentifiedByRSEM", col_3, "# IdentifiedByBoth", col_4,
                             "coef_Pearson", col_5, "pvalue_Pearson", col_6, "rho_Spearman", col_7, "pvalue_Spearman", col_8), ncol = 8, byrow = F)
write.table(comparison.table, file = paste0(outdir, "/comparisonOFmethods_isoformLevel_TPM.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
cat("The isoform-level comparison of quantification methods is done!\n")
