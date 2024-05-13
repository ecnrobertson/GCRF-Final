library(data.table)
library(dplyr)
library(tidyverse)

#This is condensed code to just pull out the high p-value value SNPs from each LMM analysis

# FEATHER ANALYSIS --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/lmm/lmm_feather/")
feather_traits <- c("pen_num", "pen_length", "plum_num", "plum_length", "nodes", "GCRF_feather_PC")
for(trait in feather_traits) {
  width <-read_delim(paste0(trait,"_lmm.assoc.txt"),delim="\t") %>%  
    mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral"))
  out<-width %>% filter(outlier=="outlier")
  if (nrow(out)!= 0) write.table(out, file=paste0(trait,"_outlier.txt"),row.names=F,quote=F,sep="\t")
}

# BEAK ANALYSIS --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/lmm/lmm_bill/")
beak_traits <- c("beak_depth", "beak_width", "culmen_end_length", "nare_length", "GCRF_bill_PC")
for(trait in beak_traits) {
  width <-read_delim(paste0(trait,"_lmm.assoc.txt"),delim="\t") %>%  
    mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral"))
  out<-width %>% filter(outlier=="outlier")
  if (nrow(out)!= 0) write.table(out, file=paste0(trait,"_outlier.txt"),row.names=F,quote=F,sep="\t")
}

# POP LEVEL --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/lmm_poplevel/")
beak_traits <- c("feather_PC_fam1", "feather_PC_fam2")
for(trait in beak_traits) {
  width <-read_delim(paste0(trait,".lmm.assoc.txt"),delim="\t") %>%  
    mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral"))
  out<-width %>% filter(outlier=="outlier")
  if (nrow(out)!= 0) write.table(out, file=paste0(trait,"_outlier.txt"),row.names=F,quote=F,sep="\t")
}
