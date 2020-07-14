library(GenomicRanges)
library(dplyr)
library(dm3)
library(vioplotx)

load("RData/bam.spc.ub.tss.gr.RData", verbose = T)
load("RData/lamin.damid.granges.RData", verbose = T)
load("RData/spc.spec.aly.indep.genes.gr.RData", verbose = T)
load("RData/tpms.RData", verbose = T)

# Combine GRanges with DamID scores into one
a <- df.from.GRanges(scl.lam.pr)
a$log2spc <- scl.lam.pr$log2damid
b <- df.from.GRanges(spg.lam.pr)
b$log2spg <- spg.lam.pr$log2damid
ab <- merge(b, a, by = c("chr", "start", "end"), all = T)
ab$start <- ab$start + 1
spg.spc.pr <- makeGRangesFromDataFrame(ab, keep.extra.columns = T)

# Make GRanges from TPM data
tpm.tss.gr <- makeGRangesFromDataFrame(tpms.all %>% 
                                     mutate(start = tss,
                                            end = tss) %>% 
                                     select(-tss),
                                   keep.extra.columns = T)

# Let's do it first for all spc-specific genes

# prepare subset of genes that aren't spc-specific but are expressed in SpC
spc.exp.no.spec <- tpm.tss.gr[(tpm.tss.gr$TPM_SpC > 1 | tpm.tss.gr$TPM_bam > 1) &
                                !(tpm.tss.gr$id %in% spc.tss.gr$id)]
# Filter spc-spec genes so that they don't have genes from that subset closer
# than 3 kb
spc.tss.gr.u <- subsetByOverlaps(spc.tss.gr,
                                 spc.exp.no.spec,
                                 maxgap = 3001,
                                 ignore.strand = T,
                                 invert = T)
# Next step of filtering is throw out genes that are closer than 3 kb to
# each other
spc.spec.d <- distanceToNearest(spc.tss.gr.u)@elementMetadata$distance
spc.tss.gr.u <- spc.tss.gr.u[spc.spec.d > 3000]

spc.spec.ov <- subsetByOverlaps(spg.spc.pr,
                                spc.tss.gr.u, ignore.strand = T)


median(spc.spec.ov$log2spg, na.rm = T)         
median(spc.spec.ov$log2spc, na.rm = T)
spc.spec.ov$diff <- spc.spec.ov$log2spc - spc.spec.ov$log2spg

wilcox.test(spc.spec.ov$diff, mu = 0, alt = "l")
median(spc.spec.ov$diff, na.rm = T)

wilcox.test(spc.spec.ov$log2spg,
            spc.spec.ov$log2spc,
            alt = "g")

# now for aly-independent

# Again, remove all expressing genes that are closer than 3 kb

spc.exp.no.aly <- tpm.tss.gr[(tpm.tss.gr$TPM_SpC > 1 | tpm.tss.gr$TPM_bam > 1) &
                               !(tpm.tss.gr$id %in% spc.aly.indep.gr.1$id)]

spc.aly.indep.gr.u <- subsetByOverlaps(spc.aly.indep.gr.1,
                                       spc.exp.no.aly,
                                       ignore.strand = T, 
                                       maxgap = 3001,
                                       invert = T)
saigu.dtn <- distanceToNearest(spc.aly.indep.gr.u)@elementMetadata$distance
spc.aly.indep.gr.u <- spc.aly.indep.gr.u[saigu.dtn > 3000]

aly.indep.ov <- subsetByOverlaps(spg.spc.pr,
                                 spc.aly.indep.gr.u,
                                 ignore.strand = T)

aly.indep.ov$diff <- aly.indep.ov$log2spc - aly.indep.ov$log2spg

median(aly.indep.ov$log2spg, na.rm = T)         
median(aly.indep.ov$log2spc, na.rm = T)
median(aly.indep.ov$diff, na.rm = T)

wilcox.test(aly.indep.ov$log2spg,
            aly.indep.ov$log2spc,
            alt = "g")
wilcox.test(aly.indep.ov$diff,
            mu = 0, alt = "l")


# For ubiquitously-expressed as control
load("RData/spc.and.ub.tss.grs.RData", verbose = T)

spc.exp.no.ub <- tpm.tss.gr[(tpm.tss.gr$TPM_SpC > 1 | tpm.tss.gr$TPM_bam > 1) &
                              !(tpm.tss.gr$id %in% ub.tss.gr$id)]
ub.tss.gr.u <- subsetByOverlaps(ub.tss.gr,
                                spc.exp.no.ub,
                                maxgap = 3001,
                                ignore.strand = T, 
                                invert = T)
utgu.dtn <- distanceToNearest(ub.tss.gr.u)@elementMetadata$distance
ub.tss.gr.u <- ub.tss.gr.u[utgu.dtn > 3000]


ub.ov <- subsetByOverlaps(spg.spc.pr, ub.tss.gr.u, ignore.strand = T)
ub.ov$diff <- ub.ov$log2spc - ub.ov$log2spg


median(ub.ov$diff, na.rm = T)
vioplotx(ub.ov$diff)

p.v <- sapply(list(spc.spec.ov$diff,
                   aly.indep.ov$diff,
                   ub.ov$diff), function(x){
                     format(wilcox.test(x, mu = 0, alt = "l")$p.value,
                            digits = 2)
                   })
med <- sapply(list(spc.spec.ov$diff,
                   aly.indep.ov$diff,
                   ub.ov$diff), function(y){
                     format(median(y, na.rm = T), digits = 2)
                   })

wilcox.test(ub.ov$diff, mu = 0, alt = "l")

pdf("plots/spc.and.aly.log2.score.diff.pdf")
  vioplotx(spc.spec.ov$diff, aly.indep.ov$diff,
           ub.ov$diff,
           col = c("#66A5A5", "#D4D411", "#069420"),
           names = c("SpC-spec", "Aly-indep", "Ubiq-expressed"))
  title(ylab = "log2(Dam-Lam/Dam)",
        main = "log2 SpC - SpG")
  text(c(0.7,1.7,2.7), y = c(4,4,4),
       paste("p.v. =", p.v))
  text(c(0.8, 1.8, 2.8), y = as.numeric(med),
       med)
  
  
dev.off()
  
pdf("plots/spc.and.aly.log2.score.pdf", height = 8, width = 12)
  par(mfrow = c(1,2))
  vioplotx(spc.spec.ov$log2spg, spc.spec.ov$log2spc,
           col = c("#9999DD", "#555588"), names = c("SpG", "SpC"))
  title(ylab = "log2(Dam-Lam/Dam",
        main = "All spc-specific")
  
  vioplotx(aly.indep.ov$log2spg, aly.indep.ov$log2spc,
           col = c("#9999DD", "#555588"), names = c("SpG", "SpC"),
           ylab = "test")
  title(ylab = "log2(Dam-Lam/Dam",
        main = "Aly-independent spc-specific")
dev.off()


