library(dplyr)
library(reshape2)
library(ggplot2)
library(tibble)
library(dm3)

# Analyse RNA-seq of larval testes (WT, Lam B + Lam C KD, ELYS KD)
load("RData/tximport.RData", verbose = T)
load("RData/lartest.DESeq.RData", verbose = T)

# Merge difseq data and TPM data
ab <- txi.lartest$abundance %>% as.data.frame()
tpms <- data.frame(
  id = row.names(ab),
  WT = rowMeans(ab[, 1:3]),
  ElysKD = rowMeans(ab[, 4:6]),
  LamKD = rowMeans(ab[, 7:9])
) 
tpms <- merge(dm3.genes, tpms, by = "id") %>% 
  filter(chr %in% euc.chroms) %>% 
  mutate(chrom = ifelse(chr == "X", "X", "A"),
         chrom = factor(chrom, levels = c("X", "A")))
# Remove LamB, LamC and Elys from analysis:
tpms <- tpms[!(tpms$gene_name %in% c("CG14215", "Lam", "LamC")), ]
# Add pseudo-counts to genes with zero expression
tpms[, 8:10][tpms[, 8:10] == 0] <- min(tpms[, 8:10][tpms[, 8:10] > 0])
tpms.m <- melt(tpms, measure.vars = c("WT", "ElysKD", "LamKD"))
tpms.chrom.sp <- split(tpms.m, list(tpms.m$variable, tpms.m$chrom))
# Wilcox pairwise tests
# WT vs ElysKD, ElysKD vs LamKD, LamKD vs WT for only X chromosome and
# autosomes
mapply(function(x,y) format(wilcox.test(tpms.chrom.sp[[x]]$value,
                                        tpms.chrom.sp[[y]]$value)$p.value,
                            digits = 3),
       1:6,
       c(2,3,1, 5, 6, 4))

# Now for the ratios of expressions

# Ratios of gene expressions on X and autosomes


tpms.1 <- tpms %>% mutate(el.wt = ElysKD/WT,
                          lam.wt = LamKD/WT)


boxplot(el.wt ~ chrom, data = tpms.1,
        outline = F)
tapply(tpms.1$el.wt, tpms.1$chrom, median, na.rm = T)
# This plot was intented to use in an article
boxplot(lam.wt ~ chrom, data = tpms.1,
        outline = F)
abline(h = 1, lty = 2, col = "red")
tapply(tpms.1$lam.wt, tpms.1$chrom, median, na.rm = T)
wilcox.test(lam.wt ~ chrom, data = tpms.1, alt = 'g')
wilcox.test(tpms.1[tpms.1$chrom == "X", ]$lam.wt, mu = 1, alt = "g")
wilcox.test(tpms.1[tpms.1$chrom == "A", ]$lam.wt, mu = 1, alt = "g")

wilcox.test(el.wt ~ chrom, data = tpms.1, alt = "g")
table(tpms.1$chrom)

# Subset genes in elys KD that are differentially expressed via padj < 0.05
# and their TPM ratio between elys KD and WT is more than 2 or less than 1/2
res.df.elys.sig <- res.df.elys %>% filter(padj < 0.05)
tpms.el.s <- tpms %>% filter(abs(log2(ElysKD/WT)) > 1,
                             id %in% res.df.elys.sig$id) %>% 
  mutate(el.wt = ElysKD/WT)
# Melt to plot boxplots more easily
tpms.el.s.m <- melt(tpms.el.s, measure.vars = c("WT", "ElysKD"))
# Check the medians of ElysKD TPM to WT TPM ratios across X and autosomes
tapply(tpms.el.s$el.wt, tpms.el.s$chrom, median)
# See how many of up- and down-regulated genes are on X and on A
tpms.el.s %>% group_by(chrom) %>% 
  summarize(`WT > KD` = sum(WT > ElysKD),
            `WT <= KD` = sum(WT <= ElysKD)) %>% t()

tpms.lam.s %>% group_by(chrom) %>% 
  summarize(`WT > KD` = sum(WT > LamKD),
            `WT <= KD` = sum(WT <= LamKD)) %>% column_to_rownames("chrom") %>% 
  t() %>% prop.table()

tpms %>% group_by(chrom) %>% 
  summarize(`WT > KD` = sum(WT > LamKD),
            `WT <= KD` = sum(WT <= LamKD)) %>% column_to_rownames("chrom") %>% 
  t() %>% chisq.test()



# Do the same for lamins KD
res.df.lam.sig <- res.df.lam %>% filter(padj < 0.05)
tpms.lam.s <- tpms %>% filter(abs(log2(LamKD/WT)) > 1,
                             id %in% res.df.lam.sig$id) %>% 
  mutate(lam.wt = LamKD/WT)
tpms.lam.s.m <- melt(tpms.lam.s, measure.vars = c("WT", "LamKD"))

tpms.lam.s %>% group_by(chrom) %>% 
  summarize(`WT > KD` = sum(WT > LamKD),
            `WT <= KD` = sum(WT <= LamKD)) %>% column_to_rownames("chrom") %>% 
  t() %>% prop.table()

tpms %>% group_by(chrom) %>% 
  summarize(`WT > KD` = sum(WT > LamKD),
            `WT <= KD` = sum(WT <= LamKD)) %>% column_to_rownames("chrom") %>% 
  t() %>% chisq.test()
# Dot plots

# Elys KD


pdf("plots/difex.elys.scatter.pdf")
ggplot(tpms %>% filter(id %in% res.df.elys.sig$id,
                       WT > 0, ElysKD > 0) %>%
         mutate(chrom = ifelse(chr == "X", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(ElysKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()
dev.off()

# Ugly

# Lam KD
pdf("plots/difex.lam.scatter.pdf")
ggplot(tpms %>% filter(id %in% res.df.lam.sig$id, LamKD > 0, WT >0) %>%
         mutate(chrom = ifelse(chr == "X", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(LamKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()
dev.off()


# Check SpC-specific genes


spc.genes <- scan("data/sclspec.fbgns.2.txt",
                  what = "char")


ggplot(tpms %>% filter(id %in% spc.genes, ElysKD > 0, WT >0) %>%
         mutate(chrom = ifelse(chr == "X", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(ElysKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()

boxplot((tpms %>% filter(id %in% spc.genes))$ElysKD,
        (tpms %>% filter(id %in% spc.genes))$WT,
        outline = F)
table((tpms %>% filter(id %in% spc.genes))$chr)

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% c("2L", "2R", "3L", "3R")))$WT,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "X"))$WT,
        outline = F, names = c("A", "X"), ylab = "TPM",
        main = "SpC-specific genes\nWT")

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% euc.chroms[-5]))$ElysKD,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "X"))$ElysKD,
        outline = F, names = c("A", "X"), ylab = "TPM",
        main = "SpC-specific genes\nElys KD")

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% euc.chroms[-5]))$LamKD,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "X"))$LamKD,
        outline = F, names = c("A", "X"), ylab = "TPM",
        main = "SpC-specific genes\nLamin KD")

# Check ratios for spc specific genes

tpms.2 <- tpms %>% filter(id %in% spc.genes) %>% 
  mutate(el.wt = ElysKD/WT,
                          lam.wt = LamKD/WT,
                          chrom = ifelse(chr == "X", "X", "A"),
         chrom = factor(chrom, levels = c("X", "A")))
tpms.2.m <- melt(tpms.2, measure.vars = c("WT", "ElysKD", "LamKD"))

tpms.2.sp <- split(tpms.2.m, list(tpms.2.m$variable, tpms.2.m$chrom))

tapply(tpms.2$el.wt, tpms.2$chrom, median, na.rm = T)

boxplot(lam.wt ~ chrom, data = tpms.2,
        outline = F)
tapply(tpms.2$lam.wt, tpms.2$chrom, median, na.rm = T)

wilcox.test(lam.wt ~ chrom, data = tpms.2, alt = "g")
wilcox.test(el.wt ~ chrom, data = tpms.2, alt = "g")


# Save data
save(tpms, tpms.m, spc.genes,
     tpms.chrom.sp,
     tpms.el.s, tpms.el.s.m,
     tpms.lam.s, tpms.lam.s.m,
     res.df.elys.sig, res.df.lam.sig,
     tpms.1, tpms.2, tpms.2.m, tpms.2.sp,
     file = "RData/tpms.analysis.RData")
