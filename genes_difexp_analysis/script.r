library("DESeq2")
library("dplyr")

# Считываем данные счетных матриц по генам и изоформам

I_df <- read.csv("I_count_mtr_genes.csv", sep = "\t")
II_df <- read.csv("II_count_mtr_genes.csv", sep = "\t")
III_df <- read.csv("III_count_mtr_genes.csv", sep = "\t")
I_iso_df <- read.csv("I_count_mtr_iso.csv", sep = "\t")
II_iso_df <- read.csv("II_count_mtr_iso.csv", sep = "\t")
III_iso_df <- read.csv("III_count_mtr_iso.csv", sep = "\t")


# Подготовка массивов имен для столбцов счетных матриц и метаданных

names_I <- c("gene_id", "I_69", "I_70", "I_71", "I_72","I_73", "I_74", "I_75",
             "I_76", "I_77", "I_78", "I_79",
             "I_80", "I_81", "I_82", "I_83", "I_84", "I_85", "I_86", "I_87")
names_II <- c("gene_id", "II_69", "II_70", "II_71", "II_72","II_73", "II_74", "II_75",
              "II_76", "II_77", "II_78", "II_79",
              "II_80", "II_81", "II_82", "II_83", "II_84", "II_85", "II_86", "II_87")
names_III <- c("gene_id", "III_69", "III_70", "III_71", "III_72","III_73", "III_74", "III_75",
               "III_76", "III_77", "III_78", "III_79",
               "III_80", "III_81", "III_82", "III_83", "III_84", "III_85", "III_86", "III_87")

names(I_df) <- names_I
names(II_df) <- names_II
names(III_df) <- names_III
names(I_iso_df) <- names_I
names(II_iso_df) <- names_II
names(III_iso_df) <- names_III

cond1 <- rep("I", 19)
cond2 <- rep("II", 19)
cond3 <- rep("III", 19)
cond2_ref <- rep("II", 13)

# Создание df для анализа
#  Тут случилась печаль, можно посмотреть на II_df и заметить какие-то проблемы с данными (71-73, 79-81)
# Складывается ощущение, что были какие-то проблемы с секвениррованием, в дальнейшем анализе мы их не используем

II_df_reformed <- II_df[, c(1,2,3,7:11, 15:20)]

I_II_genes_ref <- data.frame(I_df, II_df_reformed[c(2:14)])
df <- mutate_all(I_II_genes_ref[, c(2:33)], as.integer)
I_II_genes_ref <- data.frame(I_II_genes_ref[,1], df, row.names = "I_II_genes_ref...1.")

I_II_genes <- data.frame(I_df, II_df[c(2:20)])
df <- mutate_all(I_II_genes[, c(2:39)], as.integer)
I_II_genes <- data.frame(I_II_genes[,1], df, row.names = "I_II_genes...1.")

all_genes <- data.frame(I_df, II_df[c(2:20)], III_df[c(2:20)])
df <- mutate_all(all_genes[, c(2:58)], as.integer)
all_genes <- data.frame(all_genes[,1], df, row.names = "all_genes...1.")

I_III_genes <- data.frame(I_df, III_df[c(2:20)])
df <- mutate_all(I_III_genes[, c(2:39)], as.integer)
I_III_genes <- data.frame(I_III_genes[,1], df, row.names = "I_III_genes...1.")

II_III_genes <- data.frame(II_df, III_df[c(2:20)])
df <- mutate_all(II_III_genes[, c(2:39)], as.integer)
II_III_genes <- data.frame(II_III_genes[,1], df, row.names = "II_III_genes...1.")

II_III_genes_ref <- data.frame(II_df_reformed, III_df[c(2:20)])
df <- mutate_all(II_III_genes_ref[, c(2:33)], as.integer)
II_III_genes_ref <- data.frame(II_III_genes_ref[,1], df, row.names = "II_III_genes_ref...1.")


# Всякие действия с изоформами (в другом коде)

# Подготовка меты

meta_df_I_II <- data.frame(row.names =  c(names(I_II_genes)),
                           condition = as.factor(c(cond1,cond2)))
meta_df_II_III <- data.frame(row.names =  names(II_III_genes),
                             condition = as.factor(c(cond2,cond3)))
meta_df_I_III <- data.frame(row.names =  names(I_III_genes),
                            condition = as.factor(c(cond1,cond3)))

meta_df_I_II_ref <- data.frame(row.names =  c(names(I_II_genes_ref)),
                           condition = as.factor(c(cond1,cond2_ref)))
meta_df_II_III_ref <- data.frame(row.names =  names(II_III_genes_ref),
                                 condition = as.factor(c(cond2_ref,cond3)))

meta_all <- data.frame(row.names =  names(all_genes),
                       condition = as.factor(c(cond1,cond2, cond3)))
# Make DESeq dataset
# Из-за проблем с II_df, здесь пока адекватно работает только dds1_3, поэтому запускать лучше только его
# В итоговом анализе не учитывались 6 образцов второй точки, для этого были созданы исправленные DESeqDataSet (ddsn_n_ref)

dds1_2 <- DESeqDataSetFromMatrix(countData = I_II_genes,
                              colData = meta_df_I_II,
                              design = ~ condition)

dds2_3 <- DESeqDataSetFromMatrix(countData = II_III_genes,
                              colData = meta_df_II_III,
                              design = ~ condition)

dds1_3 <- DESeqDataSetFromMatrix(countData = I_III_genes,
                              colData = meta_df_I_III,
                              design = ~ condition)

dds1_2_ref <- DESeqDataSetFromMatrix(countData = I_II_genes_ref,
                                 colData = meta_df_I_II_ref,
                                 design = ~ condition)
dds2_3_ref <- DESeqDataSetFromMatrix(countData = II_III_genes_ref,
                                 colData = meta_df_II_III_ref,
                                 design = ~ condition)

dds_all <- DESeqDataSetFromMatrix(countData = all_genes,
                                  colData = meta_all,
                                  design = ~ condition)

# Подсчет экспрессии / с аналогичными проблемами для 2-го дф исправленными в _ref

dds1_2 <- DESeq(dds1_2)
res1_2 <- results(dds1_2)

dds2_3 <- DESeq(dds2_3)
res2_3 <- results(dds2_3)

dds1_3 <- DESeq(dds1_3)
res1_3 <- results(dds1_3)

dds_all <- DESeq(dds_all)
res <- results(dds_all)

dds1_2_ref <- DESeq(dds1_2_ref)
res1_2_ref <- results(dds1_2_ref, alpha = 0.05)

dds2_3_ref <- DESeq(dds2_3_ref)
res2_3_ref <- results(dds2_3_ref, alpha = 0.05)
res2_3_ref1 <- results(dds2_3_ref)

res1_3_1 <- results(dds1_3, alpha = 0.05)

# Получаем список всех изменившихся генов
# Далее мы планируем анализировать гены в MGigDB, поэтому для точек с наибольшими различиями 
# (при сравнении 1 и 3 условий), так как значимо изменившихся получилось достаточно много, выбрали 
# из них наиболее изменившие экспрессию

p_low_res1_3 <- subset(res1_3_1, res1_3_1$padj < 0.05 )
p_low_res1_31 <- subset(res1_3, res1_3$padj < 0.05 & abs(res1_3$log2FoldChange) > 1 )
p_low_res1_3_upper <- subset(res1_3, res1_3$padj < 0.05 & res1_3$log2FoldChange > 0 )
p_low_res1_3_lower <- subset(res1_3, res1_3$padj < 0.05 & res1_3$log2FoldChange < 0 )
p_low_res1_3_upper1 <- subset(res1_3, res1_3$padj < 0.05 & res1_3$log2FoldChange > 1 )
p_low_res1_3_lower1 <- subset(res1_3, res1_3$padj < 0.05 & res1_3$log2FoldChange < -0.5 )
nrow(p_low_res1_3_upper1[ order( p_low_res1_3_upper1$log2FoldChange ), ])
nrow(p_low_res1_3[ order( p_low_res1_3$log2FoldChange ), ]) # 4436
gene_list1_3 <- row.names(p_low_res1_3)
write(gene_list1_3, file = "gene_list1_3.txt")

gene_list1_3_upper1 <- row.names(p_low_res1_3_upper1)
write(gene_list1_3_upper1, file = "gene_list1_3_upper1.txt")
gene_list1_3_lower1 <- row.names(p_low_res1_3_lower1)
write(gene_list1_3_lower1, file = "gene_list1_3_lower1.txt")

gene_list1_31 <- row.names(p_low_res1_31)
write(gene_list1_31, file = "gene_list1_31.txt")

p_low_res1_2_ref <- subset(res1_2_ref, res1_2_ref$padj < 0.05 )
p_low_res1_2_ref_upper <- subset(res1_2_ref, res1_2_ref$padj < 0.05 & res1_2_ref$log2FoldChange > 0 )
p_low_res1_2_ref_lower <- subset(res1_2_ref, res1_2_ref$padj < 0.05 & res1_2_ref$log2FoldChange < 0 )
nrow(p_low_res1_2_ref_lower[ order( p_low_res1_2_ref_lower$log2FoldChange ), ])
nrow(p_low_res1_2_ref[ order( p_low_res1_2_ref$log2FoldChange ), ]) # 1113
gene_list1_2_ref <- row.names(p_low_res1_2_ref)
write(gene_list1_2_ref, file = "gene_list1_2.txt")

gene_list1_2_ref_upper <- row.names(p_low_res1_2_ref_upper)
write(gene_list1_2_ref_upper, file = "gene_list1_2_upper.txt")

gene_list1_2_ref_lower <- row.names(p_low_res1_2_ref_lower)
write(gene_list1_2_ref_lower, file = "gene_list1_2_lower.txt")

p_low_res2_3_ref <- subset(res2_3_ref, res2_3_ref$padj < 0.05 )
p_low_res2_3_ref_upper <- subset(res2_3_ref, res2_3_ref$padj < 0.05 & res2_3_ref$log2FoldChange > 0 )
p_low_res2_3_ref_lower <- subset(res2_3_ref, res2_3_ref$padj < 0.05 & res2_3_ref$log2FoldChange < 0 )
nrow(p_low_res2_3_ref_lower[ order( p_low_res2_3_ref_lower$log2FoldChange ), ])

nrow(p_low_res2_3_ref[ order( p_low_res2_3_ref$log2FoldChange ), ]) # 183
gene_list2_3_ref <- row.names(p_low_res2_3_ref)
write(gene_list2_3_ref, file = "gene_list2_3_a05.txt")

gene_list2_3_ref_upper <- row.names(p_low_res2_3_ref_upper)
write(gene_list2_3_ref_upper, file = "gene_list2_3_upper.txt")

gene_list2_3_ref_lower <- row.names(p_low_res2_3_ref_lower)
write(gene_list2_3_ref_lower, file = "gene_list2_3_lower.txt")

p_low_res2_3_ref1 <- subset(res2_3_ref1, res2_3_ref1$padj < 0.1 )
nrow(p_low_res2_3_ref1[ order( p_low_res2_3_ref1$log2FoldChange ), ]) # 460
gene_list2_3_ref1 <- row.names(p_low_res2_3_ref1)
write(gene_list2_3_ref1, file = "gene_list2_3_a1.txt")

p_low_res_all <- subset(res, res$padj < 0.05 )
nrow(p_low_res_all[ order( p_low_res_all$log2FoldChange ), ]) # 4333
gene_list_all<- row.names(p_low_res_all)
write(gene_list_all, file = "gene_list_all.txt")

# Посмотрели на экспрессию гена с наибольшей разницей между 1 и 3 точками
# Некоторые визуализации
library("ggplot2")
# Можно так:
plotCounts(dds1_3, gene=which.min(res1_3$padj), intgroup="condition")

# Или так:
d <- plotCounts(dds1_3, gene=which.min(res1_3$padj), intgroup="condition", 
                returnData=TRUE)
plt <- ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
plt
# Построили график главных компонент

vsd1_3 <- vst(dds1_3, blind=FALSE)
head(assay(vsd1_3), 3)

# С ggplot
pcaData1_3 <- plotPCA(vsd1_3, intgroup=c("condition"), returnData=TRUE)
percentVar1_3 <- round(100 * attr(pcaData1_3, "percentVar"))
ggplot(pcaData1_3, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar1_3[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1_3[2],"% variance")) + 
  coord_fixed()

# Без ggplot
plotPCA(vsd1_3, intgroup=c("condition"))

sampleDists <- dist(t(assay(vsd1_3)))

library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd1_3$condition, vsd1_3$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Далее мы планируем провести анализ с помощью genequery и msigdb
# Для этого мы получили ncbi id  из готовых ensg
# get stable id

library(clusterProfiler)
library(org.Hs.eg.db)

ensg = gene_list1_2_ref
# убираем номер версии
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_2 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_2 <- ncbi_conv1_2$SYMBOL
write(ncbi_id1_2, 'ncbi1_2.txt')

ensg = gene_list1_3
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_3 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
View(ncbi_id1_2)
ncbi_id1_3 <- ncbi_conv1_3$SYMBOL
write(ncbi_id1_2, 'ncbi1_2.txt')

ensg = gene_list2_3_ref
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv2_3 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
View(ncbi_id1_2)
ncbi_id2_3 <- ncbi_conv2_3$SYMBOL
write(ncbi_id2_3, 'ncbi2_3.txt')

ensg = gene_list1_2_ref_upper
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_2_upper <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_2_upper <- ncbi_conv1_2_upper$SYMBOL
write(ncbi_id1_2_upper, 'ncbi1_2_upper.txt')

ensg = gene_list1_2_ref_lower
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_2_lower <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_2_lower <- ncbi_conv1_2_lower$SYMBOL
write(ncbi_id1_2_lower, 'ncbi1_2_lower.txt')

ensg = gene_list1_3_upper
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_3_upper <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_3_upper <- ncbi_conv1_3_upper$SYMBOL
write(ncbi_id1_3_upper, 'ncbi1_3_upper.txt')

ensg = gene_list1_3_lower
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_3_lower <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_3_lower <- ncbi_conv1_3_lower$SYMBOL
write(ncbi_id1_3_lower, 'ncbi1_3_lower.txt')

ensg = gene_list1_3_upper1
# remove version number
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_3_upper1 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_3_upper1 <- ncbi_conv1_3_upper1$SYMBOL
write(ncbi_id1_3_upper1, 'ncbi1_3_upper1.txt')

ensg = gene_list1_3_lower1
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_3_lower1 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_3_lower1 <- ncbi_conv1_3_lower1$SYMBOL
write(ncbi_id1_3_lower1, 'ncbi1_3_lower1.txt')

ensg = gene_list2_3_ref_upper
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv2_3_upper <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id2_3_upper <- ncbi_conv2_3_upper$SYMBOL
write(ncbi_id2_3_upper, 'ncbi2_3_upper.txt')

ensg = gene_list2_3_ref_lower
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv2_3_lower <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id2_3_lower <- ncbi_conv2_3_lower$SYMBOL
write(ncbi_id2_3_lower, 'ncbi2_3_lower.txt')

ensg = gene_list1_31
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
ncbi_conv1_31 <- bitr(ensg.no_version, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
ncbi_id1_31 <- ncbi_conv1_31$SYMBOL
write(ncbi_id1_31, 'ncbi1_31.txt')

# Построили volcanoplot, визуализировали изменившиеся гены при сравнении трех условий 
#  1-3, 2-3, 1-2

ggplot(data=res1_3_plt, aes(x=-log2FoldChange, y=-pvalue)) + geom_point()
res1_3_plt <- as.data.frame(res1_3)
res2_3_plt <- as.data.frame(res2_3_ref)
res1_2_plt <- as.data.frame(res1_2_ref)

alpha <- 0.05 # Threshold on the adjusted p-value
cols1 <- densCols(res1_3_plt$log2FoldChange, -log10(res1_3_plt$pvalue))
pl1 <- plot(res1_3_plt$log2FoldChange, -log10(res1_3_plt$padj), col=cols, panel.first=grid(),
     main="I vs III", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6, xlim=c(-2,2), ylim=c(-0.1,7))
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

cols2 <- densCols(res1_2_plt$log2FoldChange, -log10(res1_2_plt$pvalue))
pl2 <- plot(res1_2_plt$log2FoldChange, -log10(res1_2_plt$padj), col=cols, panel.first=grid(),
     main="I vs II", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6, xlim=c(-2,2), ylim=c(-0.1,7))
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

cols3 <- densCols(res2_3_plt$log2FoldChange, -log10(res2_3_plt$pvalue))
pl3 <- plot(res2_3_plt$log2FoldChange, -log10(res2_3_plt$padj), col=cols, panel.first=grid(),
     main="II vs III", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6, xlim=c(-2,2), ylim=c(-0.1,7))
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res2_3_plt$log2FoldChange) > 2.2 & res2_3_plt$padj < 0.01 
text(res2_3_plt$log2FoldChange[gn.selected],
     -log10(res2_3_plt$padj)[gn.selected],
     lab=rownames(res2_3_plt)[gn.selected ], cex=0.4)

# Далее мы использовали подготовленные ранее файлы со списками генов, изменивших экспрессию,
# при сравнениии трех условий попарно, чтобы определить, подходят ли они под какие-либо уже известные
# функциональные списки генов.
# Результаты анализа были занесены в две таблицы по геенам и изоформам.
# Кроме того с помощью phantasus  был проведен кластерный анализ и визуализация


