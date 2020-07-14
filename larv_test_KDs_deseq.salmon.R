library(GenomicRanges)
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)
library(gplots)
library(tximport)
library(edgeR)
library(reshape2)
library(openxlsx)
library(dm3)

# Load output from salmon (our RNA-seq analysis from larvae)
load("RData/tximport.RData", verbose = T)
ddsTxi <- DESeqDataSetFromTximport(txi.lartest, colData = samples, design =~ conditions)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]
ddsTxi$conditions <- relevel(ddsTxi$conditions, ref = "wt")

dds <- DESeq(ddsTxi)
res.elys <- results(dds, contrast = c("conditions", "elysKD", "wt"))
res.lam <- results(dds, contrast = c("conditions", "lamKD", "wt"))
resultsNames(dds)

resLFC <- lfcShrink(dds, coef = "conditions_lamKD_vs_wt", type = "apeglm",
                    lfcThreshold = 1)
pdf("./plots/MA.plot.pdf")
DESeq2::plotMA(resLFC, alpha = 0.05, ylim=c(-6, 6))
dev.off()



res.df.lam <- as.data.frame(res.lam) %>% mutate(id = rownames(res.lam))
res.df.lam <- merge(genes, res.df.lam, by = "id", all.y = T)

save(res.df.lam, file = "RData/lartest.DESeq.RData")
