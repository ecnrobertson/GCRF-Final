# 2 June 2017
# Genomic landscape of divergence in south Pacific silvereyes
# R script used to reorder vcf file based on satsuma output

# ------------------------------------------------------------

setwd("~/Desktop/GCRF/GCRF-Final/figures/figures_data/")

# STEP 1: Read in and massage the satsuma and scaffold data

library(tidyverse)
#created in previous script
scaff_ord <- read_csv("GCRF_scaffold_order_from_ZFinch_Un.csv") %>%
  rename(scaffold = sca)
scaff_lengths <- read_tsv("scaffold_lengths_clean",col_names =F) %>% rename(scaffold=X1,length=X2)

missing_lengths <- anti_join(scaff_ord, scaff_lengths)

missing_ords <- anti_join(scaff_lengths, scaff_ord)


scaffs <- scaff_ord %>%
  left_join(scaff_lengths) %>%
  rename(ZFCHROM = chr, CHROM = scaffold)
scaffs

# vcf <- read_tsv("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/lmm/GCRF_bill_PC.lmm.assoc.txt", comment = "##", progress = FALSE) %>% rename(CHROM=chr, POS=ps)
# head(vcf)
vcf <- read_tsv("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/lmm_poplevel/feather_PC_fam2.lmm.assoc.txt", comment = "##", progress = FALSE) %>% rename(CHROM=chr, POS=ps)
head(vcf)

combo <- scaffs %>%
  left_join(vcf)
head(combo)
tail(combo)
combo <- combo[complete.cases(combo),]

#same order as ped file
#combo <- vcf %>% left_join(scaffs)
#tail(combo)
options(scipen = 999)

zf_ified <- combo %>%
  mutate(ZFPOS = {
    ML = floor(mean.loc)  # some temp variables to make it easier to express
    Lo2 = floor(length/2)
    L = length
    ifelse(sca.ori == 1,
           ML - Lo2 + POS,           # forward orientation
           ML - Lo2 + (L - POS))     # reverse orientation
  }) %>%
  dplyr::select(ZFCHROM, ZFPOS, everything()) %>%
  mutate(ZFCHROM = factor(ZFCHROM, levels = unique(ZFCHROM))) %>%   # this is to get them to sort correctly
  arrange(ZFCHROM, ZFPOS)

tail(zf_ified)


#zf_ified$CHROM <- NULL
#zf_ified$POS <- NULL
#zf_ified$mean.loc <- NULL
#zf_ified$sca.ori <- NULL
#zf_ified$length <- NULL

#colnames(zf_ified)[1] <- "CHROM"
#colnames(zf_ified)[2] <- "POS"
#dim(zf_ified)
#zf_ified %>% write.table("~/GCRF_wing_ulmm_snps.map",row.names=F,quote=F,sep="\t")


#Filter so that negative positions and unmapped scaffolds are removed:
zf_ified <- filter(zf_ified, POS >0)
#zf_ified <- na.omit(zf_ified) ##only if there's no NA in dataset except what you introduce
#zf_ified <- zf_ified[!grepl("_Un", zf_ified$`CHROM`),]
head(zf_ified)
dim(zf_ified)
library(pgirmess)
zf_ified %>% 
  write_delim("~/Desktop/GCRF/GCRF-Final/figures/figures_data/feather_PC_fam2_ulmm_snps_ZFchr_assoc.txt", quote = "none" ,delim = "\t", col_names = T)

# Note: To make into a functional vcf file need to add header. This is done by adding the
# following line to the top of the file:
# "##fileformat=VCFv4.2"
# This was done in atom.
