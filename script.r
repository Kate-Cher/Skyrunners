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

# Создание df для анализа
#  Тут случилась печаль, можно посмотреть на II_df и заметить какие-то проблемы с данными (71-73, 79-81)
# Складывается ощущение, что были какие-то проблемы с секвениррованием

I_II_genes <- data.frame(I_df, II_df[c(2:20)])
df <- mutate_all(I_II_genes[, c(2:39)], as.integer)
I_II_genes <- data.frame(I_II_genes[,1], df, row.names = "I_II_genes...1.")
names(I_II_genes)

I_III_genes <- data.frame(I_df, III_df[c(2:20)])
df <- mutate_all(I_III_genes[, c(2:39)], as.integer)
I_III_genes <- data.frame(I_III_genes[,1], df, row.names = "I_III_genes...1.")

II_III_genes <- data.frame(II_df, III_df[c(2:20)])
df <- mutate_all(II_III_genes[, c(2:39)], as.integer)
II_III_genes <- data.frame(II_III_genes[,1], df, row.names = "II_III_genes...1.")

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

# Подсчет экспрессии / с аналогичными проблемами для 2-го дф

dds1_2 <- DESeq(dds1_2)
res1_2 <- results(dds1_2)

dds2_3 <- DESeq(dds2_3)
res2_3 <- results(dds2_3)

dds1_3 <- DESeq(dds1_3)
res1_3 <- results(dds1_3)
design(dds1_3)

res1_3_1 <- results(dds1_3, alpha = 0.05)
summary(res1_3_1)

# Получаем список всех изменившихся генов
p_low_res <- subset(res1_3_1, res1_3_1$padj < 0.05 )
nrow(p_low_res[ order( p_low_res$log2FoldChange ), ]) # 4436
gene_llist <- row.names(p_low_res)
save(gene_llist, file = "gene_list.txt",ascii = T)

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
