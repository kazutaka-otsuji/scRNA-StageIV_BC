#------------------------------------------------------------------------------
# Code_for_scRNA-StageIV_BC
#------------------------------------------------------------------------------
library(data.table)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrastr)

### Correlation analysis; inferCNV vs bulk DNA-seq ###

#----- import inferCNV data ----#
sampleName <- c("Pre" = "pb216",
                "Post1" = "pb275bT1-kbr",
                "Post2" = "pb275T2",
                "Meta" = "pb308-1")

metadata_list <- list.files("data/metadata/")
metadata <- lapply(metadata_list, function(x){
  out <- read.csv(paste0("data/metadata/",x), row.names = 1, header = T)
  return(out)
})
names(metadata) <- gsub("_metadata.txt", "", metadata_list)

infercnv_mtx <- mclapply(names(sampleName), function(x){
  df <- read.table(paste0("data/",x,"/infercnv.observations.txt"), sep = " ", header = T, row.names = 1)
  df <- t(df)
  rownames(df) <- gsub("\\.", "-", rownames(df))
  return(df)
}, mc.cores = 4)
names(infercnv_mtx) <- names(sampleName)

#check
infercnv_mtx$Pre[c(1:5),c(1:5)]

#----- import cnvkit data ----#
cnvkit_cns <- mclapply(seq_along(sampleName), function(i){
  df <- fread(paste0("data/cnvkit/", as.character(sampleName)[i], ".cs.rmdup.cns"))
  out <- GRanges(seqnames = df$chromosome, IRanges(start = df$start, end = df$end), log2cov = df$log2)
  return(out)
}, mc.cores = 4)
names(cnvkit_cns) <- names(sampleName)

#segment to gene
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
genes_gr$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes_gr$gene_id, keytype = "ENTREZID", column = "SYMBOL")

cnvkit_geneAmp <- lapply(cnvkit_cns, function(x){
  fo1 <- findOverlaps(genes_gr, x)
  out <- data.frame(gene = genes_gr$SYMBOL[queryHits(fo1)], log2cov = x$log2cov[subjectHits(fo1)])
  return(out)
})

#----- compare infercnv and cnvkit ----#
#cancer cell only
idx1 <- list(Pre = intersect(rownames(infercnv_mtx$Pre), rownames(metadata$Pre_cancer)),
             Post1 = intersect(rownames(infercnv_mtx$Post1), rownames(metadata$Post1_cancer)),
             Post2 = intersect(rownames(infercnv_mtx$Post2), rownames(metadata$Post2[which(metadata$Post2$Celltype == "Cancer-cell"),])),
             Meta = intersect(rownames(infercnv_mtx$Meta), rownames(metadata$Meta[which(metadata$Meta$Celltype == "Cancer-cell"),])))

corr_df <- lapply(names(sampleName), function(i){
  df1 <- infercnv_mtx[[i]][idx1[[i]], ] %>% colMeans %>% as.data.frame()
  
  df2 <- cnvkit_geneAmp[[i]]
  df2 <- df2[!duplicated(df2$gene),]
  df2 <- df2[!is.na(df2$gene),]
  rownames(df2) <- df2$gene
  
  idx2 <- intersect(rownames(df1), rownames(df2)) %>% sort(.,)
  out <- data.frame(gene = idx2, infercnv = df1[idx2,1], cnvkit = df2[idx2,2])
  return(out)
})
names(corr_df) <- names(sampleName)

lapply(corr_df, function(x) cor.test(x$infercnv, x$cnvkit) %>% .$estimate) %>% unlist
#Pre.cor Post1.cor Post2.cor  Meta.cor 
#0.6363692 0.7161570 0.6955246 0.7108559 
#all test: p-value < 2.2e-16

p1 <- lapply(names(corr_df), function(i) 
  ggplot(corr_df[[i]], aes(x = infercnv, y = cnvkit)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
  scale_fill_viridis_c(trans = "log") + labs(x = "InferCNV derived score [Modified Expression]", y = "bulkDNA-seq derived CN [log2 coverage depth]") + 
  xlim(0.9, 1.1) + ylim(-2,2) + geom_smooth(method = "lm", color = "red", size = 0.5) + ggtitle(i) +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 1, lty = "dotted")
  )
pdf("output/Plots/CorrelationPlot_InferCNV_vs_CNVkit.pdf", width = 5, height = 6)
p1
dev.off()
lapply(names(corr_df), function(i) write.csv(corr_df[[i]], paste0("output/Tables/CorrelationTable_InferCNV_vs_CNVkit_", i, ".csv")))

#----- within Pre sample (C9R vs rest of cells) ----#
table(metadata$Pre_cancer$sample.original.clusters)
# Pre-RNA-C9R_Cancer 34 cells

#mean correlation
idx1 <- rownames(metadata$Pre_cancer)[which(metadata$Pre_cancer$sample.original.clusters == "Pre-RNA-C9R_Cancer")] #34 cells
idx2 <- rownames(metadata$Pre_cancer)[which(metadata$Pre_cancer$sample.original.clusters != "Pre-RNA-C9R_Cancer")] #2426 cells

idx1 <- intersect(idx1, rownames(infercnv_mtx$Pre))
idx2 <- intersect(idx2, rownames(infercnv_mtx$Pre))

Pre_df1 <- data.frame(C9RMean = colMeans(infercnv_mtx$Pre[idx1,]), OthersMean = colMeans(infercnv_mtx$Pre[idx2,]))
p2 <- ggplot(Pre_df1, aes(x = C9RMean, y = OthersMean)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
  scale_fill_viridis_c(trans = "log") + labs(x = "Pre-C9R InferCNV score", y = "Pre-Others InferCNV score") + 
  geom_smooth(method = "lm", color = "red", size = 0.5) +
  geom_hline(yintercept = 1, lty = "dotted") + geom_vline(xintercept = 1, lty = "dotted")

PreMeanDiff1 <- c(Pre_df1$C9RMean - Pre_df1$OthersMean)
names(PreMeanDiff1) <- rownames(Pre_df1)
plot(sort(PreMeanDiff1))

#statistical test, wilcox rank sum test
Pre_Pval <- mclapply(c(1:ncol(infercnv_mtx$Pre)), function(i){
  tmp <- wilcox.test(infercnv_mtx$Pre[idx1,i], infercnv_mtx$Pre[idx2,i])
  out <- tmp$p.value
  return(out)
}, mc.cores = 10)
Pre_Pval <- unlist(Pre_Pval)
Pre_FDR <- p.adjust(Pre_Pval, method = "fdr")

PreStatTable <- data.frame(row.names = colnames(infercnv_mtx$Pre),
                           MeanDiff = PreMeanDiff1,
                           pvalue = Pre_Pval,
                           FDR = Pre_FDR)
PreStatTable$minuslog10FDR <- -log10(PreStatTable$FDR)
p3 <- ggplot(PreStatTable, aes(x = MeanDiff, y = minuslog10FDR)) + geom_point_rast(size = 0.4) + 
  ArchR::theme_ArchR() + labs(x = "MeanDiff (Pre-C9R - Pre-Others)", y = "-log10(FDR)") +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted")
  # geom_hline(yintercept = 2, lty = "dotted") + geom_vline(xintercept = c(-0.1, 0.1), lty = "dotted") + 

pdf("output/Plots/CompareC9RvsOthers.pdf", width = 5, height = 5)
p2
p3
dev.off()
write.csv(PreStatTable, "output/Tables/CompareC9RvsOthers_StatTable.csv")

#----- within Pre sample (Post-RNA-C6 vs rest of Post1 cells) ----#
table(metadata$Post1_cancer$Post1_RNA1)
# Post1-RNA-C6 366 cells

#mean correlation
idx3 <- rownames(metadata$Post1_cancer)[which(metadata$Post1_cancer$Post1_RNA1 == "Post1-RNA-C6")] #366 cells
idx4 <- rownames(metadata$Post1_cancer)[which(metadata$Post1_cancer$Post1_RNA1 != "Post1-RNA-C6")] #1671 cells

idx3 <- intersect(idx3, rownames(infercnv_mtx$Post1))
idx4 <- intersect(idx4, rownames(infercnv_mtx$Post1))

Post1_df1 <- data.frame(C6Mean = colMeans(infercnv_mtx$Post1[idx3,]), OthersMean = colMeans(infercnv_mtx$Post1[idx4,]))
p4 <- ggplot(Post1_df1, aes(x = C6Mean, y = OthersMean)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
  scale_fill_viridis_c(trans = "log") + labs(x = "Post1-RNA-C6 InferCNV score", y = "Post1-Others InferCNV score") + 
  geom_smooth(method = "lm", color = "red", size = 0.5) +
  geom_hline(yintercept = 1, lty = "dotted") + geom_vline(xintercept = 1, lty = "dotted")

Post1MeanDiff1 <- c(Post1_df1$C6Mean - Post1_df1$OthersMean)
names(Post1MeanDiff1) <- rownames(Post1_df1)
plot(sort(Post1MeanDiff1))

#statistical test, wilcox rank sum test
Post1_Pval <- mclapply(c(1:ncol(infercnv_mtx$Post1)), function(i){
  tmp <- wilcox.test(infercnv_mtx$Post1[idx3,i], infercnv_mtx$Post1[idx4,i])
  out <- tmp$p.value
  return(out)
}, mc.cores = 10)
Post1_Pval <- unlist(Post1_Pval)
Post1_FDR <- p.adjust(Post1_Pval, method = "fdr")

Post1StatTable <- data.frame(row.names = colnames(infercnv_mtx$Post1),
                             MeanDiff = Post1MeanDiff1,
                             pvalue = Post1_Pval,
                             FDR = Post1_FDR)
Post1StatTable$minuslog10FDR <- -log10(Post1StatTable$FDR)
p5 <- ggplot(Post1StatTable, aes(x = MeanDiff, y = minuslog10FDR)) + geom_point_rast(size = 0.4) + 
  ArchR::theme_ArchR() + labs(x = "MeanDiff (Post1-RNA-C6 - Post1-Others)", y = "-log10(FDR)") +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted")
# geom_hline(yintercept = 2, lty = "dotted") + geom_vline(xintercept = c(-0.1, 0.1), lty = "dotted") + 

pdf("output/Plots/ComparePost1C6vsOthers.pdf", width = 5, height = 5)
p4
p5
dev.off()
write.csv(Post1StatTable, "output/Tables/ComparePost1C6vsOther_StatTable.csv")
