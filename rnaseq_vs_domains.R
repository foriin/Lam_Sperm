library(data.table)
library(dm3)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(tibble)

# Here I use RNA-seq from Laktionov et al. from adult WT and bam(-) 
# testes and compares it with our Lamin DamID data (in spermatogonia and
# spermatocytes)

# Load lists of testis-specific genes and genes expressed in SpG genes
test.spec <- scan("data/testspec.fbgns.2.txt",
                  what = character(), sep = '\n')
test.spec.id <- scan("data/testspec.geneids.2.txt",
                     what = character(), sep = '\n')
spg.spec <- scan("data/spgspec.fbgns.txt",
                 what = character(), sep = '\n')

# Some of IDs were outdated, load table with genes expressed in
# bam mutants in which updated items are shown
bamexp <- fread("data/bamexp.geneids.conv.txt") %>% 
  filter(current_symbol != "-",
         !grepl("D[^\\\\]+\\\\", current_symbol)) %>%
  filter(!(submitted_id %in% submitted_id[duplicated(bamexp$submitted_id)]))
bam.exp <- bamexp$current_id
# set spermatocyte-specific genes as those that are testes-specific, but
# don't express in SpG
scl.spec <- setdiff(test.spec, bamexp$converted_id)


# LOAD EXPRESSION DATA (from Laktionov et al.)
spg.rnaseq <- fread("data/bam_rnaseq.csv", dec = ",")
scl.rnaseq <- fread("data/wt_rnaseq.csv", dec = ",")

spg.exp <- merge(dm3.genes, spg.rnaseq %>% select(-c(1,2,3,4,6,7,8, 9, 10, 17)),
                 by.x = "gene_name", by.y = "gene_id") %>%
  filter(FPKM_status == "OK") %>% select(-FPKM_status) %>% 
  mutate(testspec = ifelse(id %in% test.spec, 1, 0),
         bamspec = ifelse(id %in% spg.spec, 1, 0),
         bamexp = ifelse(id %in% bam.exp, 1, 0))

scl.exp <- merge(dm3.genes, scl.rnaseq %>% select(-c(1,2,3,4,6,7,8, 9, 10,17)),
                 by.x = "gene_name", by.y = "gene_id") %>% 
  filter(FPKM_status == "OK") %>% select(-FPKM_status) %>% 
  mutate(sclspec = ifelse(id %in% scl.spec, 1, 0))

write.table(scl.exp, "tables/scl.fpkm.uniq.csv",quote = F,
            sep = "\t", row.names = F, dec = ",")

# LOAD DOMAINS DATA
spg.lads <- import.bed("data/LAM.SG.300nt.domains.bed")
scl.lads <- import.bed("data/LAM.SCL.300nt.domains.bed")

# LOAD ENRICHMENT DATA
lam.enr <- fread("data/Dam.Normalized.csv")
names(lam.enr) <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", names(lam.enr))

bg <- new.env()

for (i in names(lam.enr)[5:(ncol(lam.enr))]){
  assign(i,
         makeGRangesFromDataFrame(
           data.frame(lam.enr %>% select(c(2:4)), score = lam.enr[[i]]),
           keep.extra.columns = T),
         envir = bg)
}

# Starting with SpG

# subset from bam-mutant RNA-seq those genes that are testes-specific
spg.exp.ts <- spg.exp %>% filter(testspec == 1)
# create GRanges with their TSS as only coordinate
tss.gr.sp.ts <- GRanges(seqnames = Rle(spg.exp.ts$chr),
                     ranges = IRanges(start = spg.exp.ts$tss, width = 1,
                                      names = spg.exp.ts$id),
                     strand = spg.exp.ts$strand)
# Add column which indicates if a given gene is located within a LAD
spg.exp.ts$LAD <- ifelse(spg.exp.ts$id %in% subsetByOverlaps(tss.gr.sp.ts,
                                            spg.lads)@ranges@NAMES, 1, 0)
spg.exp.ts %>% group_by(LAD) %>% summarize(r = n()/nrow(spg.exp.ts), n = n())
# Subset from bam-m. RNA-seq genes that are bam-expressed and create a 
# GRanges with their TSSs
spg.exp.be <- spg.exp %>% filter(bamexp == 1)
tss.gr.sp.be <- GRanges(seqnames = Rle(spg.exp.be$chr),
                        ranges = IRanges(start = spg.exp.be$tss, width = 1,
                                         names = spg.exp.be$id))
# Add a column indicating if gene is within LAD
spg.exp.be$LAD <- ifelse(spg.exp.be$id %in% subsetByOverlaps(tss.gr.sp.be,
                                              spg.lads)@ranges@NAMES, 1, 0)
spg.exp.be %>% group_by(LAD) %>% summarize(r = n()/nrow(spg.exp.be), n = n())

# Spermatocytes

scl.exp.ts <- scl.exp %>% filter(sclspec == 1)
tss.gr.sc.ts <- GRanges(seqnames = Rle(scl.exp.ts$chr),
                        ranges = IRanges(start = scl.exp.ts$tss, width = 1,
                                         names = scl.exp.ts$id))
scl.exp.ts$LAD <- ifelse(scl.exp.ts$id %in% subsetByOverlaps(tss.gr.sc.ts,
                                            scl.lads)@ranges@NAMES, 1, 0)
scl.exp.ts %>% group_by(LAD) %>% summarize(r = n()/nrow(scl.exp.ts), n = n())
scl.nolad <- scl.exp.ts[scl.exp.ts$LAD == 0, ]$id
# Subset those genes from bam-m RNA-seq that are in interLADs in SpC and check
# how are they distributed relative to LADs in SpG
spg.exp.ts %>%  mutate(scl.nLAD = ifelse(id %in% scl.nolad, 1, 0)) %>%
  filter(scl.nLAD == 1) %>% 
  group_by(LAD) %>% summarize(r = n()/ length(scl.nolad), n = n())
# Define GRanges with interLADs in both spg and spc
gaps.spg <- gaps(spg.lads)
gaps.scl <- gaps(scl.lads)

hits <- findOverlaps(gaps.scl, gaps.spg)
# Neat trick to extract hits from findOverlaps as list containing granges
# of all hits to a query
grl <- extractList(gaps.spg, as(hits, "List"))
x <- subsetByOverlaps(gaps.scl, gaps.spg, invert = T, minoverlap = 1)

# Function to get a log2 DamID values at the TSSs of genes
damid.tss <- function(rnaseq, damid, range = 2, fun = median, ...){
  # rnaseq <- data.frame with tissue-specific
  # genes, with TPMs cut in three categories
  # via `cut(rnaseq$TPM, quantile(rnaseq$TPM, c(0, 0.33333, 0.66667, 1)),
  # labels = c("weak", "mid", "high"), include.lowest = T)` command
  # strand information also necessary
  
  tss.gr <- GRanges(
    seqnames = Rle(rnaseq$chr),
    ranges = IRanges(
      start = as.integer(rnaseq$tss),
      width = 1
    )
  )
  
  overlaps <- findOverlaps(tss.gr, damid)
  print(length(unique(overlaps@from)))
  # print(length(which(damid[overlaps@to]$index < 10)))
  vic <- sapply(1:length(overlaps@to), function(ind){
    tss.n <- overlaps@from[ind]
    tss.c <- overlaps@to[ind]
    if (rnaseq[tss.n, 5] == "+"){
      damid$score[-range:range + tss.c]
    } else {
      damid$score[range:-range + tss.c]
    }
  })
  
  names(vic) <- rnaseq$id
  # print(head(vic))
  if (range > 0){
    return(apply(vic, 1, fun, ...))
  }else{
    return(fun(vic, ...))
  }
}


tss.scores.scl.spec <- damid.tss(scl.exp.ts, bg$LAM.SCL, 0, c)
tss.scores.scl.spec <- tss.scores.scl.spec[!is.na(tss.scores.scl.spec)]
scl.spec.lt.0 <- names(tss.scores.scl.spec)[tss.scores.scl.spec < 0]

tss.scores.spg.scl.spec.lt.0 <- damid.tss(spg.exp %>%
                                            filter(id %in% scl.spec.lt.0),
                                          bg$LAM.SG, 0, c)
tss.scores.spg.scl.spec.lt.0 <- tss.scores.spg.scl.spec.lt.0[!is.na(tss.scores.spg.scl.spec.lt.0)]
sum(tss.scores.spg.scl.spec.lt.0 < 0)
# Find genes residing in newly formed interLADs in ScL

new.gaps <- GenomicRanges::setdiff(gaps.scl, gaps.spg)
tss.gr.sc <- GRanges(seqnames = Rle(scl.exp$chr),
                     ranges = IRanges(start = scl.exp$tss, width = 1,
                                      names = scl.exp$id))
# Create a column indicating if genes are in newly formed LADs
scl.exp$newgaps <- ifelse(scl.exp$id %in%
                    subsetByOverlaps(tss.gr.sc, new.gaps)@ranges@NAMES, 1, 0)
scl.exp %>% group_by(newgaps, sclspec) %>%
  summarize(n = n(), FPKM_ave.med = median(as.double(FPKM_ave)))

scl.exp.ts$newgaps <- ifelse(scl.exp.ts$id %in%
                    subsetByOverlaps(tss.gr.sc.ts, new.gaps)@ranges@NAMES, 1, 0)
scl.exp.ts %>% group_by(newgaps) %>% summarize(n = n(),
                                FPKM_ave.med = median(as.double(FPKM_ave)))

wilcox.test(scl.exp.ts$FPKM_ave[scl.exp.ts$newgaps == 0],
            scl.exp.ts$FPKM_ave[scl.exp.ts$newgaps == 1])



scl.exp %>% filter(newgaps == 1) %>% group_by(sclspec) %>% summarize(n = n())


# Same RNA-seq but analysed myself in salmon

load("RData/txi.lakt.wt.bam.RData", verbose = T)
load("RData/gene.lists.RData", verbose = T)

lakt.df <- txi.bam.wt$abundance %>% as.data.frame() %>% rownames_to_column("id")
lakt.df <- merge(lakt.df, dm3.genes, by = "id")
l.d.spg.exp <- lakt.df %>% filter(id %in% bam.exp)
sum(l.d.spg.exp$wt1 > 1)/274
