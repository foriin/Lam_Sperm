library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(dm3)
library(ggplot2)
library(reshape2)

load("RData/lamin.damid.granges.RData", verbose = T)
load("RData/tpms.analysis.RData", verbose = T)
load("RData/tpms.RData", verbose = T)
load("RData/spc.tss.tpm.RData", verbose = T)
load("RData/gene.lists.RData", verbose = T)

# Here we use SpC-specific genes that have negative or positive log2 
# Dam-Lam/Dam values in TSS in SpC and bam-expressed genes
# Both types must meet the criteria of having no expressed (TPM > 1) genes
# in a 3 kb vicinity in a given tissue

# Spermatocyte-specific genes that have negative lod2 Dam-Lam/Dam 
# values
spc.tss.gr <- dm3.tss.gr[dm3.tss.gr$id %in%
                               spc.tss.tpm$id]
# Ubuquitously expressed genes
ub.tss.gr <- dm3.tss.gr[dm3.tss.gr$id %in% ubiqs]
# Bam-expressed genes
bam.tss.gr <- dm3.tss.gr[dm3.tss.gr$id %in% bam.exp]
# Genes that have TPM > 1 in RNA-seq in SpC, excluding SpC-specific genes
spc.tpm.m.1.tss.gr <- dm3.tss.gr[dm3.tss.gr$id %in%
                                lartest.tpms$id[lartest.tpms$TPM_SpC > 1] &
                                !(dm3.tss.gr$id %in% spc.spec)]
# Genes that have TPM > 1 in RNA-seq in SpG, except for bam-expressed
spg.tpm.m.1.gr <- dm3.tss.gr[dm3.tss.gr$id %in%
                                   bam.tpms$id[bam.tpms$TPM_bam > 1] &
                               !(dm3.tss.gr$id %in% bam.exp)]
# Filter so there are only genes that do not "cluster", i.e. locate closer
# than 3kb to each other
spc.tss.gr.f <- spc.tss.gr[abs(start(spc.tss.gr) -
                                 start(spc.tss.gr)[nearest(spc.tss.gr,
                                                        ignore.strand = T)]) >
                             3000]
# spc.tss.gr.f <- subsetByOverlaps(spc.tss.gr.f, ub.tss.gr,
#                                  maxgap = 5000, ignore.strand = T,
#                                  invert = T)
bam.tss.gr.f <- bam.tss.gr[abs(start(bam.tss.gr) -
                                 start(bam.tss.gr)[nearest(bam.tss.gr,
                                                           ignore.strand = T)]) >
                             3000]
# Ubiquitous genes
ub.tss.gr.f <-  ub.tss.gr[abs(start(ub.tss.gr) -
                                 start(ub.tss.gr)[nearest(ub.tss.gr,
                                                           ignore.strand = T)]) >
                             3000]
ub.tss.gr.f <- keepSeqlevels(ub.tss.gr.f,
                             euc.chroms,
                             pruning.mode = "coarse")
ub.tss.gr.1 <- subsetByOverlaps(ub.tss.gr.f,
                                spc.tpm.m.1.tss.gr[!(spc.tpm.m.1.tss.gr$id %in%
                                                        ub.tss.gr.f$id)],
                                maxgap = 3000, ignore.strand = T,
                                invert = T)
ub.tss.gr.2 <- subsetByOverlaps(ub.tss.gr.f,
                                spg.tpm.m.1.gr[!(spg.tpm.m.1.gr$id %in%
                                                       ub.tss.gr.f$id)],
                                maxgap = 3000, ignore.strand = T,
                                invert = T)

# Filter out SpC-spec genes that are closer than 3kb to genes that have TPM > 1
# in SpC
spc.tss.gr.1 <- subsetByOverlaps(spc.tss.gr.f,
                                 spc.tpm.m.1.tss.gr[!(spc.tpm.m.1.tss.gr$id %in%
                                                        spc.tss.gr.f$id)],
                                 maxgap = 3000, ignore.strand = T,
                                 invert = T)
spc.tss.gr.2 <- subsetByOverlaps(spc.tss.gr.f,
                                 spg.tpm.m.1.gr,
                                 maxgap = 3000, ignore.strand = T,
                                 invert = T)

# The same for SpG
bam.tss.gr.1 <- subsetByOverlaps(bam.tss.gr.f,spg.tpm.m.1.gr,
                                 maxgap = 3000, ignore.strand = T,
                                 invert = T)

bam.tss.gr.2 <- subsetByOverlaps(bam.tss.gr.f,
                spc.tpm.m.1.tss.gr[!(spc.tpm.m.1.tss.gr$id %in% 
                                       bam.tss.gr.f$id)],
                                 maxgap = 3000, ignore.strand = T,
                                 invert = T)

# Function to get damid signals around TSSs
tss.surr <- function(tss.gr, signal.gr, sig.col, vic = 5){
ovs <- findOverlaps(tss.gr, signal.gr, ignore.strand = T)
ovs <- ovs[ovs@to > vic]
# print(mcols(signal.gr)[[sig.col]][1:10])
sapply(1:length(ovs), function(x){
  if(strand(tss.gr[ovs@from[x]]) %>% as.vector() == "-"){
    rev(mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)])
  } else{
    mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
  }
})

# Don't take the strandness into account

# sapply(1:length(ovs), function(x){
#     mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
#     
#   })
}

spc.sp.spc.lam.surr <- tss.surr(spc.tss.gr.1,
                                scl.lam.pr,vic = 10, sig.col = "log2damid")
plot(rowMeans(spc.sp.spc.lam.surr), type = "l")
# load Lamin DamID profile from neurons
load("RData/nrn.lam.pr.RData", verbose = T)


# Filter data; leave only euc chroms

spc.tss.gr.1 <- keepSeqlevels(spc.tss.gr.1,
                              euc.chroms)
spc.tss.gr.2 <- keepSeqlevels(spc.tss.gr.2,
                              euc.chroms)

bam.tss.gr.1 <- keepSeqlevels(bam.tss.gr.1,
                              euc.chroms)
bam.tss.gr.2 <- keepSeqlevels(bam.tss.gr.2,
                              euc.chroms)

scl.lam.pr.f <- keepSeqlevels(scl.lam.pr,
                              euc.chroms,
                              pruning.mode = "coarse")
spg.lam.pr.f <- keepSeqlevels(spg.lam.pr,
                              euc.chroms,
                              pruning.mode = "coarse")
nrn.lam.pr.gr.f <- keepSeqlevels(nrn.lam.pr.gr,
                               euc.chroms,
                               pruning.mode = "coarse")
# Generate profiles for SpC spec genes  at TSS
tss.prof.bc.1 <- mapply(tss.surr, 
       split(spc.tss.gr.1, seqnames(spc.tss.gr.1)),
       split(scl.lam.pr.f, seqnames(scl.lam.pr.f)),
       MoreArgs = list(vic = 10, sig.col = "log2damid"))

x1 <- do.call(cbind,tss.prof.bc.1) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x1.m <- melt(x1, id.var = "kb")
ggplot(x1.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()
require(msir)
test <- loess.sd(x1.m$kb, x1.m$value, span = 0.4, degree = 2)

plot(x1.m[,-2])
plot(test$x, test$y, type = "l", ylim = c(-2, 2))
lines(test$x, test$upper)
lines(test$x, test$lower)

tss.prof.1 <- do.call(cbind, tss.prof.bc.1) %>% rowMeans()%>% 
                as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))


p1 <- ggplot(tss.prof.1 %>% 
         setNames(c("log2(Dam/Dam-Lam)", "kb")),
       aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x1.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x1.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("SpC-specific genes\nLam DamID in SpC")

# Add loess smoothing

tss.prof.loess.1 <- loess.smooth(tss.prof.1$kb, tss.prof.1$., degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.1) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

# The same for spc-spec genes in SpG

tss.prof.bc.2 <- mapply(tss.surr, 
                                     split(spc.tss.gr.2, seqnames(spc.tss.gr.2)),
                                     split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                                     MoreArgs = list(vic = 10, sig.col = "log2damid"))

x2 <- do.call(cbind,tss.prof.bc.2) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x2.m <- melt(x2, id.var = "kb")
ggplot(x2.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()

tss.prof.2 <- do.call(cbind, tss.prof.bc.2) %>% rowMeans()%>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))


p2 <- ggplot(tss.prof.2 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x2.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x2.m, aes(x = kb, y = value),
            method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
            col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("SpC-specific genes\nLam DamID in SpG")

# Add loess smoothing

tss.prof.loess.2 <- loess.smooth(tss.prof.2$kb, tss.prof.2$., degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.2) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

# Spermatogonia

# Bam-expressed genes/SpG DamID
tss.prof.bc.3 <- mapply(tss.surr, 
                               split(bam.tss.gr.1, seqnames(bam.tss.gr.1)),
                               split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                               MoreArgs = list(vic = 10, sig.col = "log2damid"))
x3 <- do.call(cbind,tss.prof.bc.3) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x3.m <- melt(x3, id.var = "kb")
ggplot(x3.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.3 <- do.call(cbind, tss.prof.bc.3) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.3 <- loess.smooth(tss.prof.3$kb, tss.prof.3$., degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.3) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p3 <- ggplot(tss.prof.3 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x3.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x3.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Bam-expressed genes\nLam DamID in SpG")

# SpC-spec genes with no expressed genes in SpG closer than 3kb

tss.prof.bc.4 <- mapply(tss.surr, 
                        split(spc.tss.gr.1, seqnames(spc.tss.gr.1)),
                        split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                        MoreArgs = list(vic = 10, sig.col = "log2damid"))
x4 <- do.call(cbind,tss.prof.bc.4) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x4.m <- melt(x4, id.var = "kb")
ggplot(x4.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.4 <- do.call(cbind, tss.prof.bc.4) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
# tss.prof.4 <- apply(do.call(cbind, tss.prof.bc.4), 1, median) %>% 
#   as.data.frame()  %>% 
#   mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.4 <- loess.smooth(tss.prof.4$kb, tss.prof.4$., degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.4) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p4 <- ggplot(tss.prof.4 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x4.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x4.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("SpC-specific genes\nw/o Lam enrichment in TSS\nLam DamID in SpG")

# SpC-spec genes with negative DamID score in SpC in TSS/DamID in neurons

tss.prof.bc.5 <- mapply(tss.surr, 
                        split(spc.tss.gr.1, seqnames(spc.tss.gr.1)),
                        split(nrn.lam.pr.gr.f, seqnames(nrn.lam.pr.gr.f)),
                        MoreArgs = list(vic = 6, sig.col = "score"))
x5 <- do.call(cbind, tss.prof.bc.5) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))
x5.m <- melt(x5, id.var = "kb")
ggplot(x5.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.5 <- do.call(cbind, tss.prof.bc.5) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))

tss.prof.loess.5 <- loess.smooth(tss.prof.5$kb, tss.prof.5$., degree = 2, span = 0.5)
ggplot(do.call(cbind, tss.prof.loess.5) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p5 <- ggplot(tss.prof.5 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x5.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x5.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("SpC-specific genes\nw/o Lam enrichment in TSS\nLam DamID in Neurons")

# Use old microarray data for Kc167 cells
load("RData/kc.lam.profile.ma.RData")

tss.prof.bc.6 <- mapply(tss.surr, 
                        split(spc.tss.gr.1, seqnames(spc.tss.gr.1)),
                        split(kc.bins.lam.pr.f, seqnames(kc.bins.lam.pr.f)),
                        MoreArgs = list(vic = 6, sig.col = "score"))
x6 <- do.call(cbind, tss.prof.bc.6) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))
x6.m <- melt(x6, id.var = "kb")
ggplot(x6.m, aes(x = kb, y = value))+
  geom_smooth(method = "loess")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.6 <- do.call(cbind, tss.prof.bc.6) %>% rowMeans(na.rm = T)  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))

tss.prof.loess.6 <- loess.smooth(tss.prof.6$kb, tss.prof.6$., degree = 2, span = 0.5)
ggplot(do.call(cbind, tss.prof.loess.6) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p6 <- ggplot(tss.prof.6 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x6.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x6.m, aes(x = kb, y = value),
              method = "loess", se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("SpC-specific genes\nw/o Lam enrichment in TSS\nLam DamID in Kc167")

# Bam-expressed genes/SpC DamID
tss.prof.bc.7 <- mapply(tss.surr, 
                        split(bam.tss.gr.2, seqnames(bam.tss.gr.2)),
                        split(scl.lam.pr.f, seqnames(scl.lam.pr.f)),
                        MoreArgs = list(vic = 10, sig.col = "log2damid"))
x7 <- do.call(cbind,tss.prof.bc.7) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x7.m <- melt(x7, id.var = "kb")
ggplot(x7.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.7 <- do.call(cbind, tss.prof.bc.7) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.7 <- loess.smooth(tss.prof.7$kb, tss.prof.7$.,
                                 degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.7) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p7 <- ggplot(tss.prof.7 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x7.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x7.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Bam-expressed genes\nLam DamID in SpC")

# Ubiquitous genes in SpC
tss.prof.bc.8 <- mapply(tss.surr, 
                        split(ub.tss.gr.1, seqnames(ub.tss.gr.1)),
                        split(scl.lam.pr.f, seqnames(scl.lam.pr.f)),
                        MoreArgs = list(vic = 10, sig.col = "log2damid"))
x8 <- do.call(cbind,tss.prof.bc.8) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x8.m <- melt(x8, id.var = "kb")
ggplot(x8.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.8 <- do.call(cbind, tss.prof.bc.8) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.8 <- loess.smooth(tss.prof.8$kb, tss.prof.8$.,
                                 degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.8) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p8 <- ggplot(tss.prof.8 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x8.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x8.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Ubiq-expressed genes\nLam DamID in SpC")

# Ubiquitous genes in SpG
tss.prof.bc.9 <- mapply(tss.surr, 
                        split(ub.tss.gr.2, seqnames(ub.tss.gr.2)),
                        split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                        MoreArgs = list(vic = 10, sig.col = "log2damid"))
x9 <- do.call(cbind,tss.prof.bc.9) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x9.m <- melt(x9, id.var = "kb")
ggplot(x9.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.9 <- do.call(cbind, tss.prof.bc.9) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.9 <- loess.smooth(tss.prof.9$kb, tss.prof.9$.,
                                 degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.9) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p9 <- ggplot(tss.prof.9 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x9.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x9.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Ubiq-expressed genes\nLam DamID in SpG")

# Aly dependent spc-spec genes in spc

aly.table <- read.xlsx("tables/SpC-specific and Aly-dependent or independent genes with Lamin-DamID values_03.10.2019.xlsx")
aly.dep <- (aly.table %>% filter(aver.3P.call.in.aly.testes == "A"))$FBgn

spc.aly.dep.gr.1 <- spc.tss.gr.1[spc.tss.gr.1$id %in% aly.dep]

tss.prof.bc.10 <- mapply(tss.surr, 
                        split(spc.aly.dep.gr.1, seqnames(spc.aly.dep.gr)),
                        split(scl.lam.pr.f, seqnames(scl.lam.pr.f)),
                        MoreArgs = list(vic = 10, sig.col = "log2damid"))
x10 <- do.call(cbind,tss.prof.bc.10) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x10.m <- melt(x10, id.var = "kb")
ggplot(x10.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.10 <- do.call(cbind, tss.prof.bc.10) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.10 <- loess.smooth(tss.prof.10$kb, tss.prof.10$.,
                                 degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.10) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p10 <- ggplot(tss.prof.10 %>% 
               setNames(c("log2(Dam/Dam-Lam)", "kb")),
             aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x10.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x10.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly dependent genes\nLam DamID in SpC")

# Aly independent spc-spec genes in spc

aly.indep <- (aly.table %>% filter(aver.3P.call.in.aly.testes == "P"))$FBgn

spc.aly.indep.gr.1 <- spc.tss.gr.1[spc.tss.gr.1$id %in% aly.indep]

save(aly.dep, aly.indep, file = "RData/aly.dep.indep.genes.id.RData")
save(spc.aly.indep.gr.1, file = "RData/spc.spec.aly.indep.genes.gr.RData")

tss.prof.bc.11 <- mapply(tss.surr, 
                         split(spc.aly.indep.gr.1, seqnames(spc.aly.indep.gr.1)),
                         split(scl.lam.pr.f, seqnames(scl.lam.pr.f)),
                         MoreArgs = list(vic = 10, sig.col = "log2damid"))
x11 <- do.call(cbind,tss.prof.bc.11) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x11.m <- melt(x11, id.var = "kb")
ggplot(x11.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.11 <- do.call(cbind, tss.prof.bc.11) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.11 <- loess.smooth(tss.prof.11$kb, tss.prof.11$.,
                                  degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.11) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p11 <- ggplot(tss.prof.11 %>% 
                setNames(c("log2(Dam/Dam-Lam)", "kb")),
              aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x11.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x11.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly independent genes\nLam DamID in SpC")



# Aly dependent spc-spec genes in spG

spc.aly.dep.gr.2 <- spc.tss.gr.2[spc.tss.gr.2$id %in% aly.dep]

tss.prof.bc.12 <- mapply(tss.surr, 
                         split(spc.aly.dep.gr.2, seqnames(spc.aly.dep.gr.2)),
                         split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                         MoreArgs = list(vic = 10, sig.col = "log2damid"))
x12 <- do.call(cbind,tss.prof.bc.12) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x12.m <- melt(x12, id.var = "kb")
ggplot(x12.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.12 <- do.call(cbind, tss.prof.bc.12) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.12 <- loess.smooth(tss.prof.12$kb, tss.prof.12$.,
                                  degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.12) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p12 <- ggplot(tss.prof.12 %>% 
                setNames(c("log2(Dam/Dam-Lam)", "kb")),
              aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x12.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x12.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly dependent genes\nLam DamID in SpG")

# Aly independent spc-spec genes in spc

spc.aly.indep.gr.2 <- spc.tss.gr.2[spc.tss.gr.2$id %in% aly.indep]

tss.prof.bc.13 <- mapply(tss.surr, 
                         split(spc.aly.indep.gr.2, seqnames(spc.aly.indep.gr.2)),
                         split(spg.lam.pr.f, seqnames(spg.lam.pr.f)),
                         MoreArgs = list(vic = 10, sig.col = "log2damid"))
x13 <- do.call(cbind,tss.prof.bc.13) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))
x13.m <- melt(x13, id.var = "kb")
ggplot(x13.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.13 <- do.call(cbind, tss.prof.bc.13) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.3, length.out = 21))

tss.prof.loess.13 <- loess.smooth(tss.prof.13$kb, tss.prof.13$.,
                                  degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.13) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p13 <- ggplot(tss.prof.13 %>% 
                setNames(c("log2(Dam/Dam-Lam)", "kb")),
              aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x13.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x13.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly independent genes\nLam DamID in SpG")

# Aly independent spc-spec genes in spc
spc.tss.gr.f <- keepSeqlevels(spc.tss.gr.f, euc.chroms,
                              pruning.mode = "coarse")

spc.aly.indep.gr <- spc.tss.gr.f[spc.tss.gr.f$id %in% aly.indep]

tss.prof.bc.14 <- mapply(tss.surr, 
                         split(spc.aly.indep.gr, seqnames(spc.aly.indep.gr)),
                         split(nrn.lam.pr.gr.f, seqnames(nrn.lam.pr.gr.f)),
                         MoreArgs = list(vic = 6, sig.col = "score"))
x14 <- do.call(cbind,tss.prof.bc.14) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))
x14.m <- melt(x14, id.var = "kb")
ggplot(x14.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.14 <- do.call(cbind, tss.prof.bc.14) %>% rowMeans()  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))

tss.prof.loess.14 <- loess.smooth(tss.prof.14$kb, tss.prof.14$.,
                                  degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.14) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p14 <- ggplot(tss.prof.14 %>% 
                setNames(c("log2(Dam/Dam-Lam)", "kb")),
              aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x14.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x14.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly independent genes\nLam DamID in neurons")

# Aly independent spc-spec genes in spc
spc.tss.gr.f <- keepSeqlevels(spc.tss.gr.f, euc.chroms,
                              pruning.mode = "coarse")

spc.aly.indep.gr <- spc.tss.gr.f[spc.tss.gr.f$id %in% aly.indep]

tss.prof.bc.15 <- mapply(tss.surr, 
                         split(spc.aly.indep.gr, seqnames(spc.aly.indep.gr)),
                         split(kc.bins.lam.pr.f, seqnames(kc.bins.lam.pr.f)),
                         MoreArgs = list(vic = 6, sig.col = "score"))
x15 <- do.call(cbind,tss.prof.bc.15) %>% as.data.frame() %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))
x15.m <- melt(x15, id.var = "kb")
ggplot(x15.m, aes(x = kb, y = value))+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))


tss.prof.15 <- do.call(cbind, tss.prof.bc.15) %>% rowMeans(na.rm = T)  %>% 
  as.data.frame()  %>% 
  mutate(kb = seq(-3, by = 0.5, length.out = 13))

tss.prof.loess.15 <- loess.smooth(tss.prof.15$kb, tss.prof.15$.,
                                  degree = 2, span = 0.4)
ggplot(do.call(cbind, tss.prof.loess.15) %>% as.data.frame(), aes(x = x, y = y))+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1))

p15 <- ggplot(tss.prof.15 %>% 
                setNames(c("log2(Dam/Dam-Lam)", "kb")),
              aes(x = kb, y = `log2(Dam/Dam-Lam)`))+
  geom_line()+
  geom_smooth(linetype = 0, data = x15.m, aes(x = kb, y = value),
              method = "loess", span = 1/5)+
  stat_smooth(geom = "line",data = x15.m, aes(x = kb, y = value),
              method = "loess",se = T, span = 1/5, alpha = 0.5, size = 1.5,
              col = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Aly independent genes\nLam DamID in Kc167")

# require("gridExtra")
install.packages("ggpubr")
require(ggpubr)
# pdf("plots/lam.surround.tss.diff.groups.3.pdf", width = 7, height = 7)
# grid.arrange(p1 + ylim(c(-0.5, 1.7)),
#              p3 + ylim(c(-2.5, -0.3)),
#              p2 + ylim(c(-0.5, 1.7)),
#              p4 + ylim(c(-2.5, -0.3)),
#              p5 + ylim(c(-0.5, 1.7)),
#              p6 + ylim(c(-2.5, -0.3)),
#              ncol = 2)
mp <- ggarrange(p1, p2+ylab(""), p3, p7+ylab(""), p8, p9 +ylab(""),
             ncol = 2, nrow = 2)
ggexport(mp, filename = "plots/lam.surround.tss.diff.groups.4.pdf")
dev.off()


# Make it on a whole-genome level


zs <- mget(ls(pattern = "^z\\d?$"))
all.tss.surr <- lapply(zs, function(df) data.frame(score = rowMeans(df[,-6]),
                                    kb = df$kb))
all.tss.surr <- all.tss.surr[c("z", "z2", "z1", "z3", "z4", "z5")]
tissues <- rep(c("SpC", "Neurons", "SpG"), each = 2)
gene.class <- rep(c("SpC-specific", "Ubiq expressed"), 3)
limits = rep(list(c(0, 0.7), c(-1.7, 0)), 3)

plots <- lapply(c(1,5,3), function(x){
  ggplot(all.tss.surr[[x]],
         aes(x = kb, y = score))+
    geom_line()+
    theme_bw()+
    scale_x_continuous(breaks = seq(-3, 3, 1))+
    ylab("log2(Dam/Dam-Lam)")+
    # ylim(limits[[x]])+
    ggtitle(paste(gene.class[x], "genes\nLam DamID in", tissues[x]))
})

pdf("plots/spc.spec.lam.surr.spc.spg.nrn.pdf", width = 6, height = 6)
grid.arrange(grobs = plots, ncol = 2)
dev.off()

save(bam.tss.gr, spc.tss.gr, ub.tss.gr, file = "RData/bam.spc.ub.tss.gr.RData")

pdf("plots/bam.exp.spc.tss.prof.pdf")
p7
dev.off()
pdf("plots/ubiq.exp.spc.tss.prof.pdf")
p8
dev.off()
pdf("plots/ubiq.exp.spg.tss.prof.pdf")
p9
dev.off()

pdf("plots/aly.spc.spg.pdf", width = 10, height = 10)
grid.arrange(p10,
             p11+ylab(""),
             p12,
             p13+ylab(""),
             ncol = 2)
dev.off()

pdf("plots/aly.ind.nrn.kc.pdf", width = 10, height = 5)
grid.arrange(p14,
             p15+ylab(""),
             ncol = 2)
dev.off()
