if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("apeglm")
BiocManager::install("pheatmap")
library(tximport)
library(apeglm)
library(pheatmap)
library(DESeq2)
setwd('D:/sky_project/isoforms_result/')

# comparison I and II conditions
#files <- list.files()[c(10, 14:19, 29, 33:38)]
files <- list.files()[c(1:21, 25:29, 33:38)]
file_names = c()
for (file in files){
  file_names <- c(file_names, substr(file, 1, nchar(file) - 4))
}


names(files) <- file_names
txi_rsem <- tximport(files, type = 'rsem', txIn = T, txOut = T)
head(txi_rsem$length, n = 10)
head(txi_rsem$counts, n = 10)
head(txi_rsem$abundance, n = 10)
txi_rsem$counts <- round(txi_rsem$counts)


datacol <- matrix(c(rep('I', times = 19), rep('II', times = 13)), nrow = 32, ncol = 1)
colnames(datacol) <- c('condition')
rownames(datacol) <- colnames(txi_rsem$counts)
datacol <- data.frame(datacol)
datacol$condition <- factor(datacol$condition)

dds <- DESeqDataSetFromMatrix(countData = txi_rsem$counts,
                              colData = datacol,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
#----------------------------------------------------------------
res_without_na <- res[!is.na(res$padj), ]
res_padj <- res_without_na[!is.na(res_without_na$padj), ][res_without_na[!is.na(res_without_na$padj), ]$padj < 0.05,]
upper_transcript_ids <- rownames(res_padj[res_padj$log2FoldChange > 0, ])
down_transcript_ids <- rownames(res_padj[res_padj$log2FoldChange < 0, ])
#write(upper_transcript_ids, file = 'I_III_transcripts_up_fc.txt')
#write(down_transcript_ids, file = 'I_III_transcripts_dw_fc.txt')
resLFC <- lfcShrink(dds, coef="condition_II_vs_I", type="apeglm")
resLFC
#
#
#
#
#
# comparison I and III conditions
files2 <- list.files()[c(1:19, 39:57)]
file_names2 = c()
for (file in files2){
  file_names2 <- c(file_names2, substr(file, 1, nchar(file) - 4))
}


names(files2) <- file_names2
txi_rsem2 <- tximport(files2, type = 'rsem', txIn = T, txOut = T)
head(txi_rsem3$length, n = 10)
head(txi_rsem3$counts, n = 10)
head(txi_rsem3$abundance, n = 10)
txi_rsem2$counts <- round(txi_rsem2$counts)


datacol2 <- matrix(c(rep('I', times = 19), rep('III', times = 19)), nrow = 38, ncol = 1)
colnames(datacol2) <- c('condition')
rownames(datacol2) <- colnames(txi_rsem2$counts)
datacol2 <- data.frame(datacol2)
datacol2$condition <- factor(datacol2$condition)

dds2 <- DESeqDataSetFromMatrix(countData = txi_rsem2$counts,
                              colData = datacol2,
                              design = ~ condition)
#keep <- rowSums(counts(dds2)) >= 10
#dds2 <- dds2[keep,]
dds2 <- DESeq(dds2)
res2 <- results(dds2)


resLFC2 <- lfcShrink(dds2, coef="condition_III_vs_I", type="apeglm")
resLFC2
#----------------------------------------------------------------
res_without_na2 <- res2[!is.na(res2$padj), ]
res_padj2 <- res_without_na2[!is.na(res_without_na2$padj), ][res_without_na2[!is.na(res_without_na2$padj), ]$padj < 0.05,]
upper_transcript_ids2 <- rownames(res_padj2[res_padj2$log2FoldChange > 0, ])
down_transcript_ids2 <- rownames(res_padj2[res_padj2$log2FoldChange < 0, ])

#write(upper_transcript_ids2, file = 'I_III_transcripts_up_fc.txt')
#write(down_transcript_ids2, file = 'I_III_transcripts_dw_fc.txt')
#
#
#
#
#
# comparison II and III condition
files3 <- list.files()[c(20:21, 25:29, 33:57)]
file_names3 = c()
for (file in files3){
  file_names3 <- c(file_names3, substr(file, 1, nchar(file) - 4))
}


names(files3) <- file_names3
txi_rsem3 <- tximport(files3, type = 'rsem', txIn = T, txOut = T)
head(txi_rsem3$length, n = 10)
head(txi_rsem3$counts, n = 10)
head(txi_rsem3$abundance, n = 10)
txi_rsem3$counts <- round(txi_rsem3$counts)


datacol3 <- matrix(c(rep('II', times = 13), rep('III', times = 19)), nrow = 32, ncol = 1)
colnames(datacol3) <- c('condition')
rownames(datacol3) <- colnames(txi_rsem3$counts)
datacol3 <- data.frame(datacol3)
datacol3$condition <- factor(datacol3$condition)

dds3 <- DESeqDataSetFromMatrix(countData = txi_rsem3$counts,
                               colData = datacol3,
                               design = ~ condition)
keep <- rowSums(counts(dds3)) >= 10
dds3 <- dds3[keep,]
dds3 <- DESeq(dds3)
res3 <- results(dds3)
#----------------------------------------------------------------
res_without_na3 <- res3[!is.na(res3$padj), ]
res_padj3 <- res_without_na3[!is.na(res_without_na3$padj), ][res_without_na3[!is.na(res_without_na3$padj), ]$padj < 0.05,]
upper_transcript_ids3 <- rownames(res_padj3[res_padj3$log2FoldChange > 0, ])
down_transcript_ids3 <- rownames(res_padj3[res_padj3$log2FoldChange < 0, ])

#write(upper_transcript_ids3, file = 'II_III_transcripts_up_fc.txt')
#write(down_transcript_ids3, file = 'II_III_transcripts_dw_fc.txt')

resLFC3 <- lfcShrink(dds3, coef="condition_III_vs_II", type="apeglm")
resLFC3



#---------------------------------------------------------------
# NEXT STAGE
vsd <- vst(dds3, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:10000]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)


pcaData1_3 <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar1_3 <- round(100 * attr(pcaData1_3, "percentVar"))
ggplot(pcaData1_3, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar1_3[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1_3[2],"% variance")) + 
  coord_fixed() + 
  geom_text(aes(label=name))
#ggsave('PCA_IvsIII_isoforms.png', plot = last_plot())
