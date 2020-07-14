library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(dm3)

# 20.04.20
# This script was used to test non-randomness of association of
# Aly-independent spermatocyte-specific genes with LADs in Spermatogonia
# and interLADs in Spermatocytes

# Function to shuffle domains and TSS
# !!!! Please provide path to BedTools and genome file to function below
# Example of genome file is in the "data" folder !!!!

bedTools.shuffle.gr <- function(gr, opt.string="-chrom"){
  
  bed.file <- tempfile()
  
  shuf <- tempfile()
  
  
  options(scipen = 99)
  
  export.bed(gr, bed.file)
  command = paste("/path/to/bedtools shuffle -i", bed.file,
                  "-g /path/to/genome/file", opt.string, "|",
                  "/path/to/bedtools sort -i - >", shuf)
  # cat(command, "\n")
  tryCatch(system(command))
  gr.shuf <- import.bed(shuf)
  unlink(bed.file); unlink(shuf)
  return(gr.shuf)
}
# Load domain data and modify their starts to help conversion into BED
load("RData/lamin.damid.granges.RData")
start(scl.lam.hmm) <- start(scl.lam.hmm) + 1
start(spg.lam.hmm) <- start(spg.lam.hmm) + 1
# load("RData/tpms.analysis.RData")
# Load Aly-independent subset of spc-specific genes
load("RData/spc.spec.aly.indep.genes.gr.RData")
# Number of TSSs of those genes that are in iLADs in SpC
ix.spc <- length(subsetByOverlaps(spc.aly.indep.gr.1,
                 scl.lam.hmm,
                 ignore.strand = T,
                 invert = T))
# Do the permutation test x10k
perm.test.spc <- mclapply(1:10000, function(somebodyoncetoldmetheworldsgonnarollme){
  length(subsetByOverlaps(bedTools.shuffle.gr(spc.aly.indep.gr.1),
                          scl.lam.hmm,
                          ignore.strand = T, invert = T)) > ix.spc
}, mc.cores = 14)
perm.test.spc <- do.call(c, perm.test.spc)
# receive p-value
sum(perm.test.spc)/10000
# Number of TSSs of those genes that are in LADs in SpG
ix.spg <- length(subsetByOverlaps(spc.aly.indep.gr.1,
                                  spg.lam.hmm,
                                  ignore.strand = T))
# Perm test x10k
perm.test.spg <- mclapply(1:10000, function(somebodyoncetoldmetheworldsgonnarollme){
  length(subsetByOverlaps(bedTools.shuffle.gr(spc.aly.indep.gr.1),
                          spg.lam.hmm,
                          ignore.strand = T)) > ix.spg
}, mc.cores = 14)
perm.test.spg <- do.call(c, perm.test.spg)
# p-value
sum(perm.test.spg)/10000
