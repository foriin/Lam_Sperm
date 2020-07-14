library(stringr)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(dm3)

load("RData/bam.spc.ub.tss.gr.RData", verbose = T)
load("RData/lamin.damid.granges.RData", verbose = T)

# Function to shuffle domains

# !!!!!! Please write full path to BedTools in your system and genome file
# with sizes of chromosomes (consult with 
# https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html, you can
# find my genome file in the 'data' folder)!!!!!!
bedTools.shuffle.tss <- function(query, subject, opt.string="-chrom"){
  
  bed.file.1 <- tempfile()
  # print(bed.file.1)
  bed.file.2 <- tempfile()
  # print(bed.file.2)
  shuf.1 <- tempfile()
  shuf.2 <- tempfile()
  
  
  options(scipen = 99)
  
  write.table(query, file = bed.file.1, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(subject, file = bed.file.2, quote = F, sep = "\t", col.names = F, row.names = F)
  command = paste("/path/to/bedtools shuffle -i", bed.file.1,
                  "-g /path/to/genome/file", opt.string, "|",
                  "/path/to/bedtools sort -i - >", shuf.1, ";",
                  "/path/to/bedtools shuffle -i", bed.file.2,
                  "-g /path/to/genome/file", opt.string, "|",
                  "/path/to/bedtools sort -i - >", shuf.2)
  # cat(command, "\n")
  tryCatch(system(command))
  qu <- import.bed(shuf.1)
  subj <- import.bed(shuf.2)
  q.x.s <- subsetByOverlaps(qu, subj)
  unlink(bed.file.1); unlink(bed.file.2); unlink(shuf.2); unlink(shuf.1)
  return(length(q.x.s))
}

# spc.spec.x.ilads <- length(subsetByOverlaps(makeGRangesFromDataFrame(spc.spec.bed),
                                            # (makeGRangesFromDataFrame(spc.ilads.bed))))

# a <- sapply(1:10000, function(x) bedTools.shuffle.tss(spc.spec.bed, spc.ilads.bed) >
#          spc.spec.x.ilads)
# spc.all.x.ilads <- length(subsetByOverlaps(makeGRangesFromDataFrame(spc.all.bed),
#                                            (makeGRangesFromDataFrame(spc.ilads.bed))))
# b <- sapply(1:11000,
#             function(x) bedTools.shuffle.tss(spc.all.bed, spc.ilads.bed) >
#               spc.all.x.ilads)
# 
# spg.all.x.ilads <- length(subsetByOverlaps(makeGRangesFromDataFrame(spg.all.bed),
#                                            (makeGRangesFromDataFrame(spg.ilads.bed))))
# 
# c <- sapply(1:11000,
#             function(x) bedTools.shuffle.tss(spg.all.bed, spg.ilads.bed) >
#               spg.all.x.ilads)
# sum(c)

bam.exp.bed <- df.from.GRanges(bam.tss.gr) %>% 
  mutate(start = start - 1)
spc.spec.bed <- df.from.GRanges(spc.tss.gr) %>% 
  mutate(start = start - 1)
spg.ilads.bed <- df.from.GRanges(gaps(spg.lam.hmm)) %>% 
  mutate(start = start - 1)
spc.ilads.bed <- df.from.GRanges(gaps(scl.lam.hmm)) %>% 
  mutate(start = start - 1)
spc.spec.spc.ilads.bed <- df.from.GRanges(subsetByOverlaps(spc.tss.gr,
                                                           gaps(scl.lam.hmm),
                                                           ignore.strand = T)) %>% 
  mutate(start = start - 1)

bam.exp.x.ilads <- length(subsetByOverlaps(bam.tss.gr, spg.lam.hmm, invert = T))

a <- sapply(1:10000, function(x) bedTools.shuffle.tss(bam.exp.bed,
                                                      spg.ilads.bed) >
         bam.exp.x.ilads)
sum(a)

spc.spec.x.ilads <- length(subsetByOverlaps(spc.tss.gr, scl.lam.hmm, invert = T))

b <- sapply(1:10000, function(x) bedTools.shuffle.tss(spc.spec.bed,
                                                      spc.ilads.bed) >
              spc.spec.x.ilads)
sum(b)/10000

spc.spec.ilads.x.ilads.spg <- length(subsetByOverlaps(
  subsetByOverlaps(spc.tss.gr,
                   gaps(scl.lam.hmm),
                   ignore.strand = T),
  spg.lam.hmm, invert = T))

c <- sapply(1:10000, function(x) bedTools.shuffle.tss(spc.spec.spc.ilads.bed,
                                                      spg.ilads.bed) >
              spc.spec.ilads.x.ilads.spg)
sum(c)/10000

# Compare TSSs of spc-specific genes that are in iLADs both in SpC and SpG 
# with iLADs in Kc and NRN

load("RData/spc.spec.x.ilads.spg.spc.kc.nrn.RData")

spc.tss.ilads.bed <- df.from.GRanges(spc.tss.s.2) %>% 
  mutate(start = start - 1)
nrn.ilads.bed <- df.from.GRanges(nrn.ilad) %>% 
  mutate(start = start - 1)
kc.ilads.bed <- df.from.GRanges(kc.ilad) %>% 
  mutate(start = start - 1)

d <- mclapply(1:50000, function(x) bedTools.shuffle.tss(spc.tss.ilads.bed,
                                                      nrn.ilads.bed) >
              len.1, mc.cores = 12)
sum(unlist(d))

e <- mclapply(1:10000, function(x) bedTools.shuffle.tss(spc.tss.ilads.bed,
                                                        kc.ilads.bed) >
                len.2, mc.cores = 12)
sum(unlist(e))






