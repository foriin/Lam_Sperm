library(dplyr)
library(gplots)
library(data.table)

# load data from DamID and RNA-seq which will be used to calculate correlations
load("RData/corr.RData", verbose = T)

MakeCorMat <- function(data, cormethod, use.vals = "complete.obs"){
  cormat <- as.matrix(data) %>%
    cor(method = cormethod, use = use.vals) %>% 
    round(digits = 2)
}

MainCorrelations <- function(data, use.opt="everything", corr.desc, createPDF=T,
                             file.p, ...) {  
  lapply(c("spearman", "pearson"), function(meth){
    cors <- MakeCorMat(data, meth)
    assign(paste0(corr.desc, ".", meth, ".cor"), cors)
    if (createPDF == T){
      options(warn=-1)
      pdf(file=file.path(file.p, paste0(corr.desc,  "_", meth, "_correlation_heatmap", ".pdf")), width=12, height=12)
      heatmap.2(cors, col=greenred(200), breaks=seq(from=-1, to=1, by=0.01),
                Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=cors,
                notecol="white", notecex = 1.5, margins=c(7,7),
                main=paste(meth, "'s correlation coefficients and hierarchical clustering for\n", "'",
                           sub("\\d+\\.(.*)", "\\1", corr.desc), "'", sep=""),
                cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2),
                lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="density")
      options(warn=0)
      dev.off()
    }  
  })
}

raw.spg <- raw.damid.spc.spg %>% select(5,7,6,8) %>% setNames(c(
  "DAM.1", "DAM.2", "DAM-LAM.1", "DAM-LAM.2"
))
raw.spc <- raw.damid.spc.spg %>% select(9:12) %>% setNames(c(
  "DAM.1", "DAM.2", "DAM-LAM.1", "DAM-LAM.2"))

MainCorrelations(raw.spg, corr.desc = "spermatogonia", file.p = "plots")
MainCorrelations(raw.spc, corr.desc = "spermatocytes", file.p = "plots")

dd.spg <- dist(scale(t(raw.spg)), method = "euclidean")
hc.spg <- hclust(dd.spg, method = "ward.D2")
plot(as.dendrogram(hc.spg), type = "rectangle", main = "Spermatogonia", xlab = "")

dd.spc <- dist(scale(t(raw.spc)), method = "euclidean")
hc.spc <- hclust(dd.spc, method = "ward.D2")
plot(as.dendrogram(hc.spc), type = "rectangle", main = "Spermatocytes", xlab = "")

pdf("plots/spg.spc.dendro.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(as.dendrogram(hc.spg), type = "rectangle", main = "Spermatogonia", xlab = "")
plot(as.dendrogram(hc.spc), type = "rectangle", main = "Spermatocytes", xlab = "")
dev.off()

# RNA-seq



dd.rna <- dist(scale(rnaseq.rpm), method = "euclidean")
hc.rna <- hclust(dd.rna, method = "ward.D2")
pdf("plots/rnaseq.dendro.2.pdf")
plot(as.dendrogram(hc.rna), type = "rectangle", main = "RNA-seq in larvae",
     xlab = "")
dev.off()

MainCorrelations(rnaseq.rpm, corr.desc = "larva_RNA-seq", file.p = "plots")
