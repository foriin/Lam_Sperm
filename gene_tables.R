library(dplyr)
library(data.table)
library(GenomicRanges)
library(openxlsx)
library(tibble)
library(dm3)


# Load testis- and bam-specific genes, define spermatocyte-spec genes
test.spec <- scan("data/testspec.fbgns.3.txt",
                  what = character(), sep = '\n')


ubiq.genes <- scan("data/ubiq.chiant.genes.upd.txt",
                   what = character(), sep = '\n')
# Load Laktionov RNA-seq data in bam and wt testes

load("RData/txi.lakt.wt.bam.RData", verbose = T)

# Load DamID profiles and domains

load("RData/lamin.damid.granges.RData", verbose = T)

load("RData/spc.and.ub.tss.grs.RData", verbose = T)

bam.tpms <- data.frame(
  id = rownames(txi.bam.wt$abundance),
  TPM_bam = rowMeans(txi.bam.wt$abundance[, 1:2])
) %>% remove_rownames()

bam.tpms <- merge(dm3.genes, bam.tpms, by = "id")

# Load our RNA-seq in spc

load("RData/tximport.RData", verbose = T)
load("RData/aly.dep.indep.genes.id.RData", verbose = T)

lartest.tpms <- data.frame(
  id = rownames(txi.lartest$abundance),
  TPM_SpC = rowMeans(txi.lartest$abundance[, 1:3],),
  TPM_SpC_ElysKD = rowMeans(txi.lartest$abundance[, 4:6],),
  TPM_SpC_LambcKD = rowMeans(txi.lartest$abundance[, 7:9],)
) %>% remove_rownames()

tpms.all <- merge(bam.tpms, lartest.tpms, by = "id") %>%
  mutate(testspec = ifelse(id %in% test.spec, 1, 0),
         bamexp = ifelse(id %in% test.spec & TPM_bam > 1, 1, 0),
         spcspec = ifelse(id %in% test.spec & TPM_bam < 1 & TPM_SpC > 1, 1, 0),
         ubiq = ifelse(id %in% ubiq.genes & TPM_bam > 1 & TPM_SpC > 1, 1, 0),
         aly_dep = ifelse(id %in% aly.dep & spcspec == 1, 1, 0),
         aly_indep = ifelse(id %in% aly.indep & spcspec == 1, 1, 0))
# write.xlsx(tpms.all %>% select(-10), file = "tables/Table_S3.xlsx")
spc.spec <- tpms.all$id[tpms.all$spcspec == 1]
bam.exp <- tpms.all$id[tpms.all$bamexp == 1]
# tpms.f <- tpms.all %>% filter(testspec == 1, !(bamspec < 1 & spcspec < 1) )
# test.spec <- tpms.f$id
# cat(test.spec, file = "~/IMG/Projects/LAM.SpG.SpC/RNA-seq/testspec.fbgns.3.txt",
#     sep = "\n")


# Load difex data

load("RData/lartest.DESeq.RData")

tpms.all <- tpms.all %>% 
  mutate(ElysKD_difex = ifelse(abs(log2(TPM_SpC_ElysKD/TPM_SpC)) > 1 &
                                 id %in% res.df.elys.sig$id, 1, 0),
         LambcKD_difex = ifelse(abs(log2(TPM_SpC_LambcKD/TPM_SpC)) > 1 &
                                 id %in% res.df.lam.sig$id, 1, 0))




write.xlsx(tpms.all, file = "tables/SpG_and_SpC_tpms_and_spec_genes_subsets.xlsx")

save(ubiq.genes, test.spec, bam.exp, spc.spec, file = "RData/gene.lists.RData")
save(bam.tpms, lartest.tpms, tpms.all, file = "RData/tpms.RData")


# SpC spec + log2

spc.tss.gr <- sort(spc.tss.gr, ignore.strand = T)

spc.tss.tpm <- merge(data.frame(df.from.GRanges(spc.tss.gr),
                          id = spc.tss.gr$id,
                          gene_name = spc.tss.gr$gene_name,
                          stringsAsFactors = F), 
                     lartest.tpms[, 1:2],
                     by = "id") %>% arrange(chr, start)

spc.spec.id.lam.spc <- (subsetByOverlaps(spc.tss.gr, scl.lam.hmm,
                                     ignore.strand = T))$id

spc.tss.ov.lam.spc.pr <- findOverlaps(spc.tss.gr, scl.lam.pr)
spc.tss.tpm$`log2(Dam-Lam/Dam) SpC` <- NA
spc.tss.tpm$`log2(Dam-Lam/Dam) SpC`[spc.tss.ov.lam.spc.pr@from] <- scl.lam.pr$log2damid[spc.tss.ov.lam.spc.pr@to]
spc.tss.tpm$`in Lam domain SpC` <- ifelse(spc.tss.tpm$id %in% spc.spec.id.lam.spc, 1, 0)

spc.spec.id.lam.spg <- (subsetByOverlaps(spc.tss.gr, spg.lam.hmm,
                                         ignore.strand = T))$id

spc.tss.ov.lam.spg.pr <- findOverlaps(spc.tss.gr, spg.lam.pr)
spc.tss.tpm$`log2(Dam-Lam/Dam) SpG` <- NA
spc.tss.tpm$`log2(Dam-Lam/Dam) SpG`[spc.tss.ov.lam.spg.pr@from] <- spg.lam.pr$log2damid[spc.tss.ov.lam.spg.pr@to]
spc.tss.tpm$`in Lam domain SpG` <- ifelse(spc.tss.tpm$id %in% spc.spec.id.lam.spg, 1, 0)

# Bam-expressed + log2 Lam

bam.tss.gr <- sort(bam.tss.gr, ignore.strand = T)

bam.tss.tpm <- merge(data.frame(df.from.GRanges(bam.tss.gr),
                                id = bam.tss.gr$id,
                                gene_name = bam.tss.gr$gene_name,
                                stringsAsFactors = F), 
                     bam.tpms[, c(1,8)],
                     by = "id") %>% arrange(chr, start)

bam.exp.id.lam.spg <- (subsetByOverlaps(bam.tss.gr, spg.lam.hmm,
                                         ignore.strand = T))$id

bam.tss.ov.lam.spg.pr <- findOverlaps(bam.tss.gr, spg.lam.pr)
bam.tss.tpm$`log2(Dam-Lam/Dam) SpG` <- NA
bam.tss.tpm$`log2(Dam-Lam/Dam) SpG`[bam.tss.ov.lam.spg.pr@from] <- spg.lam.pr$log2damid[bam.tss.ov.lam.spg.pr@to]
bam.tss.tpm$`in Lam domain SpG` <- ifelse(bam.tss.tpm$id %in% bam.exp.id.lam.spg,
                                          1, 0)



write.xlsx(spc.tss.tpm, "tables/spc.tss.tpm.log2.xlsx")
write.xlsx(bam.tss.tpm, "tables/bam.tss.tpm.log2.xlsx")
write.xlsx(spc.tss.tpm %>% filter(spc.tss.tpm$`log2(Dam-Lam/Dam) SpC` > 0),
           file = "tables/spc.tss.tpm.pos.log2.xlsx")
write.xlsx(spc.tss.tpm %>% filter(spc.tss.tpm$`log2(Dam-Lam/Dam) SpC` < 0),
           file = "tables/spc.tss.tpm.neg.log2.xlsx")
write.xlsx(tpms.all %>% select(-10), file = "tables/Table_S3.xlsx")

save(spc.tss.tpm, file = "RData/spc.tss.tpm.RData")
cat(spc.tss.tpm$id[spc.tss.tpm$`log2(Dam-Lam/Dam) SpC`< 0],
    file = "data/spc.spec.tss.neg.damid.fbgns.txt", sep = "\n")
cat(spc.tss.tpm$id[spc.tss.tpm$`log2(Dam-Lam/Dam) SpC` > 0],
    file = "data/spc.spec.tss.pos.damid.fbgns.txt", sep = "\n")

