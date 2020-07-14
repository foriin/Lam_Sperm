library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
# devtools::install_github("TomKellyGenetics/vioplotx")
library(vioplotx)

rm(list = ls())
chr.len <- data.frame("chr" = c("2L", "2R", "3L", "3R", "4", "X", "2LHet",
                                "2RHet", "3LHet", "3RHet", "XHet", "YHet"),
                      "len" = c(23011544, 21146708, 24543557, 27905053, 1351857,
                                22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038))
chroms <- c("2L", "2R", "3L", "3R", "X")

spg.lam.pr <- import.bedGraph("data/LAM.SG.wt.sum.norm.bedgraph") %>%
  add.chr(rev = T) %>% GRanges()

spg.lam.hmm <- import.bed("data/LAM.SG.300nt.domains.bed") %>% 
  GRanges()

scl.lam.pr <- import.bedGraph("data/LAM.SCL.wt.sum.norm.bedgraph") %>% 
  add.chr(rev = T) %>% GRanges()
names(mcols(scl.lam.pr)) <- "log2damid"

scl.lam.hmm <- import.bed("data/LAM.SCL.300nt.domains.bed") %>%
  GRanges()

spg.dom.len <-  sapply(split(width(spg.lam.hmm), seqnames(spg.lam.hmm)), sum) 
spg.dom.len <- spg.dom.len[names(spg.dom.len) %in% euc.chroms]

chromlen <- chr.lengths[names(chr.lengths) %in% euc.chroms]

scl.dom.len <-  sapply(split(width(scl.lam.hmm), seqnames(scl.lam.hmm)), sum) 
scl.dom.len <- scl.dom.len[names(scl.dom.len) %in% euc.chroms]

sum(spg.dom.len) / sum(chromlen)
sum(scl.dom.len) / sum(chromlen)

save(spg.lam.hmm, spg.lam.pr, scl.lam.hmm, scl.lam.pr, 
     file = "RData/lamin.damid.granges.RData")
