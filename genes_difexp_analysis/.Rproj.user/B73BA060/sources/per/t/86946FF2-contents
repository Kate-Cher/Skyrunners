library("DESeq2")
library("dplyr")

# 

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
# Складывается ощущение, что были какие-то проблемы с секвениррованием

II_df_reformed <- II_df[, c(1,2,3,7:11, 15:20)]

I_II_genes_ref <- data.frame(I_df, II_df_reformed[c(2:14)])
df <- mutate_all(I_II_genes_ref[, c(2:33)], as.integer)
I_II_genes_ref <- data.frame(I_II_genes_ref[,1], df, row.names = "I_II_genes_ref...1.")

I_II_genes <- data.frame(I_df, II_df[c(2:20)])
df <- mutate_all(I_II_genes[, c(2:39)], as.integer)
I_II_genes <- data.frame(I_II_genes[,1], df, row.names = "I_II_genes...1.")

I_III_genes <- data.frame(I_df, III_df[c(2:20)])
df <- mutate_all(I_III_genes[, c(2:39)], as.integer)
I_III_genes <- data.frame(I_III_genes[,1], df, row.names = "I_III_genes...1.")

II_III_genes <- data.frame(II_df, III_df[c(2:20)])
df <- mutate_all(II_III_genes[, c(2:39)], as.integer)
II_III_genes <- data.frame(II_III_genes[,1], df, row.names = "II_III_genes...1.")

II_III_genes_ref <- data.frame(II_df_reformed, III_df[c(2:20)])
df <- mutate_all(II_III_genes_ref[, c(2:33)], as.integer)
II_III_genes_ref <- data.frame(II_III_genes_ref[,1], df, row.names = "II_III_genes_ref...1.")


# Всякие действия с iso / надо доделать

# I_II_iso <- data.frame(I_iso_df, II_iso_df[c(2:20)])
# II_III_iso <- data.frame(II_iso_df, III_iso_df[c(2:20)])
# I_III_iso <- data.frame(I_iso_df, III_iso_df[c(2:20)])

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
# Make DESeq dataset
# Из-за проблем с II_df, здесь пока адекватно работает только dds1_3, поэтому запускать лучше только его

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

# Подсчет экспрессии / с аналогичными проблемами для 2-го дф

dds1_2 <- DESeq(dds1_2)
res1_2 <- results(dds1_2)

dds2_3 <- DESeq(dds2_3)
res2_3 <- results(dds2_3)

dds1_3 <- DESeq(dds1_3)
res1_3 <- results(dds1_3)

dds1_2_ref <- DESeq(dds1_2_ref)
res1_2_ref <- results(dds1_2_ref, alpha = 0.05)

dds2_3_ref <- DESeq(dds2_3_ref)
res2_3_ref <- results(dds2_3_ref, alpha = 0.05)
res2_3_ref1 <- results(dds2_3_ref)

res1_3_1 <- results(dds1_3, alpha = 0.05)

# Получаем список всех изменившихся генов
p_low_res1_3 <- subset(res1_3_1, res1_3_1$padj < 0.05 )
nrow(p_low_res1_3[ order( p_low_res1_3$log2FoldChange ), ]) # 4436
gene_list1_3 <- row.names(p_low_res1_3)
write(gene_list1_3, file = "gene_list1_3.txt")

p_low_res1_2_ref <- subset(res1_2_ref, res1_2_ref$padj < 0.05 )
nrow(p_low_res1_2_ref[ order( p_low_res1_2_ref$log2FoldChange ), ]) # 1113
gene_list1_2_ref <- row.names(p_low_res1_2_ref)
write(gene_list1_2_ref, file = "gene_list1_2.txt")

p_low_res2_3_ref <- subset(res2_3_ref, res2_3_ref$padj < 0.05 )
nrow(p_low_res2_3_ref[ order( p_low_res2_3_ref$log2FoldChange ), ]) # 183
gene_list2_3_ref <- row.names(p_low_res2_3_ref)
write(gene_list2_3_ref, file = "gene_list2_3_a05.txt")

p_low_res2_3_ref1 <- subset(res2_3_ref1, res2_3_ref1$padj < 0.1 )
nrow(p_low_res2_3_ref1[ order( p_low_res2_3_ref1$log2FoldChange ), ]) # 460
gene_list2_3_ref1 <- row.names(p_low_res2_3_ref1)
write(gene_list2_3_ref1, file = "gene_list2_3_a1.txt")

# Посмотрели на экспрессию гена с наибольшей разницей между 1 и 3 точками
library("ggplot2")
# Можно так:
plotCounts(dds1_3, gene=which.min(res1_3$padj), intgroup="condition")

# Или так:
d <- plotCounts(dds1_3, gene=which.min(res1_3$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Попробовали построить график главных компонент

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
