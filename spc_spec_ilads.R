library(dplyr)
library(GenomicRanges)
library(dm3)

# rm(list=ls())
# Load LADs and Lam profiles in SpG and SpC
load("RData/lamin.damid.granges.RData", verbose = T)

# Load RNA-seq data from SpC
load("RData/tpms.analysis.RData", verbose = T)

tss.tpms.gr <- makeGRangesFromDataFrame(tpms.1 %>%
                                          select(chr, tss, strand, id,
                                                 gene_name, WT) %>% 
                                          mutate(chr = sub("chr", "", chr)),
                                        start.field = "tss", 
                                        end.field = "tss",
                                        keep.extra.columns = T)
# Subset "expressed" genes with TPM values > 1
tss.tpms.m.1 <- tss.tpms.gr[tss.tpms.gr$WT > 1]
# See how many "expressed" genes overlap by their TSS with SpC iLADs
tss.expr.x.iLADs.spc <- length(subsetByOverlaps(tss.tpms.m.1,
                                                gaps(scl.lam.hmm),
                                                ignore.strand = T))
tss.expr.x.LADs.spc <- length(subsetByOverlaps(tss.tpms.m.1,
                                               scl.lam.hmm,
                                               ignore.strand = T))


a <- df.from.GRanges(tss.tpms.m.1)
b <- df.from.GRanges(gaps(scl.lam.hmm))
c <- df.from.GRanges(scl.lam.hmm)

spc.ilam.x.spc.exp.pv <- sum(unlist(mclapply(1:10000,
                function(x) bedTools.shuffle.tss(a, b) >
                                         tss.expr.x.iLADs.spc, mc.cores = 12)))
spc.lam.x.spc.exp.pv <- sum(unlist(mclapply(1:10000,
                                    function(x) bedTools.shuffle.tss(a, c) >
                                     tss.expr.x.LADs.spc, mc.cores = 12)))

# Load testes-specific gene lists
load("RData/gene.lists.RData", verbose = T)

tss.ubiq <- dm3.tss.gr[dm3.tss.gr$id %in% ubiq.genes]
spg.x.ubiq.ilad <- length(subsetByOverlaps(tss.ubiq, 
                                           gaps(spg.lam.hmm)))
spc.x.ubiq.ilad <- length(subsetByOverlaps(tss.ubiq, 
                                           gaps(scl.lam.hmm)))

a <- df.from.GRanges(tss.ubiq)
b <- df.from.GRanges(gaps(spg.lam.hmm))
c <- df.from.GRanges(gaps(scl.lam.hmm))

spg.ilad.x.ubiq.pv <- sum(unlist(mclapply(1:10000,
                                          function(x) bedTools.shuffle.tss(a, b) >
                                            spg.x.ubiq.ilad, mc.cores = 12)))

spc.ilad.x.ubiq.pv <- sum(unlist(mclapply(1:10000,
                                          function(x) bedTools.shuffle.tss(a, c) >
                                            spc.x.ubiq.ilad, mc.cores = 12)))


subsetByOverlaps(spc.spec.gr, scl.lam.hmm)

subsetByOverlaps(spc.spec.gr, spg.lam.hmm)

