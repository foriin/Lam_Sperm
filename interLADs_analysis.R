library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# Function to shuffle domains and TSS
# !!!! Please provide path to BedTools and genome file to function below
# Example of genome file is in the "data" folder !!!!

bedTools.shuffle.tss <- function(query, subject, opt.string="-chrom"){
  
  bed.file.1 <- tempfile()
  # print(bed.file.1)
  bed.file.2 <- tempfile()
  # print(bed.file.2)
  shuf.1 <- tempfile()
  shuf.2 <- tempfile()
  # print(shuf.1)
  
  
  options(scipen = 99)
  
  write.table(query, file = bed.file.1, quote = F, sep = "\t",
              col.names = F, row.names = F)
  write.table(subject, file = bed.file.2, quote = F, sep = "\t",
              col.names = F, row.names = F)
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

load("RData/lamin.damid.granges.RData")
load("RData/tpms.analysis.RData")

# Make TPMs gr with coordinates of tss of genes that are expressed in SpC
# (TPM > 1) from our SpC RNA-seq data
tpms.1.gr <- makeGRangesFromDataFrame(tpms.1 %>% filter(WT > 1) %>% 
                                        select(-start, -end)%>% 
                                        mutate(chr = sub("chr", "", chr)),
                                      keep.extra.columns = T, seqnames.field = "chr",
                                      start.field = "tss", end.field = "tss")
# Subset TSSs that reside in the interLADs
tpms.mt.1.iLADs <- subsetByOverlaps(tpms.1.gr, gaps(scl.lam.hmm))

scl.lam.pr.low <- scl.lam.pr[scl.lam.pr$log2damid < 0] %>% reduce()

subsetByOverlaps(tpms.1.gr, scl.lam.pr.low)

# Make TPMs gr with coordinates of tss of ScL-specific genes 
# from our SpC RNA-seq data
tpms.2.gr <- makeGRangesFromDataFrame(tpms.2 %>% filter(WT > 1) %>% 
                                        select(-start, -end)%>% 
                                        mutate(chr = sub("chr", "", chr)),
                                      keep.extra.columns = T, seqnames.field = "chr",
                         start.field = "tss", end.field = "tss")
# Subset TSSs that reside in the interLADs
tpms.2.iLADs <- subsetByOverlaps(tpms.2.gr, gaps(scl.lam.hmm), ignore.strand = T)

subsetByOverlaps(tpms.2.gr, scl.lam.pr.low, ignore.strand = T)

# ubiquitous genes

ubiqs <- scan("data/ubiq.chiant.genes.upd.txt", what = character(), sep = "\n")
tpms.ub <- tpms.1.gr[tpms.1.gr$id %in% ubiqs]

subsetByOverlaps(dm3.tss.gr[dm3.tss.gr$id %in% ubiqs], gaps(spg.lam.hmm), ignore.strand = T)
subsetByOverlaps(tpms.ub, scl.lam.pr.low)
subsetByOverlaps(tpms.ub, gaps(spg.lam.hmm), ignore.strand = T)
