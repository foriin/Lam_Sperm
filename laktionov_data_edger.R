library(dplyr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(tximport)
library(edgeR)

load("RData/txi.lakt.wt.bam.RData", verbose = T)


dgList <- DGEList(counts = txi.bam.wt$counts, genes = row.names(txi.bam.wt$counts))
dgList
head(dgList$counts)

dg.cpm <- cpm(dgList)
summary(dg.cpm)
# Filtering
ccheck <- dg.cpm > 1
summary(ccheck)

keep <- which(rowSums(ccheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))
# Normalise 
dgList <- calcNormFactors(dgList, method = "TMM")
x <- removeBatchEffect(dgList)

plotMDS(x, dim.plot = c(1,3))
# Set up model
sampleType <- c("BAM", "BAM", "WT", "WT")
sampleReplicate <- paste("S", rep(1:2,times = 2), sep="")

designMat <- model.matrix(~sampleReplicate + sampleType)
designMat

# Estimate dispersions
dgList <- estimateDisp(dgList, design = designMat)
dgList <- estimateGLMCommonDisp(dgList, design = designMat)
dgList <- estimateGLMTrendedDisp(dgList, design = designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design = designMat)

cor(txi.bam.wt$counts[,1], txi.bam.wt$counts[,2])

plotBCV(dgList)

# Differential expression
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef = 2)
# ?glmLRT
res <- topTags(lrt, n = 5000)
res

deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)


