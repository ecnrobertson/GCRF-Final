library(data.table)
library(dplyr)
library(tidyverse)

#This is condensed code to just pull out the high PIP value SNPs from each BSLMM analysis

# FEATHER ANALYSIS --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/bslmm/bslmm_feather/")
feather_traits <- c("pen_num", "pen_length", "plum_num", "plum_length", "nodes", "GCRF_feather_PC")
for(trait in feather_traits) {
  params<-fread(paste0(trait,".bslmm.param.txt"),header=T,sep="\t", data.table=F)
  params.pipsort<-params[order(-params$gamma),]
  pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
  # variants with effect in 10% MCMC samples or more
  pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
  # variants with effect in 25% MCMC samples or more
  pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
  # variants with effect in 50% MCMC samples or more
  pip50<-params.pipsort[params.pipsort$gamma>=0.50,]
 
  if (nrow(pip01) != 0) write.table(pip01, file=paste0(trait,"_pip01.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip10) != 0) write.table(pip10, file=paste0(trait,"_pip10.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip25) != 0) write.table(pip25, file=paste0(trait,"_pip25.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip50) != 0) write.table(pip50, file=paste0(trait,"_pip50.dsv"), quote=F, row.names=F, sep="\t")
}


# BEAK ANALYSIS --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/bslmm/bslmm_beak/")
beak_traits <- c("beak_depth", "beak_width", "culmen_end_length", "nare_length", "GCRF_bill_PC")
for(trait in beak_traits) {
  params<-fread(paste0(trait,".bslmm.param.txt"),header=T,sep="\t", data.table=F)
  params.pipsort<-params[order(-params$gamma),]
  pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
  # variants with effect in 10% MCMC samples or more
  pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
  # variants with effect in 25% MCMC samples or more
  pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
  # variants with effect in 50% MCMC samples or more
  pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

  if (nrow(pip01) != 0) write.table(pip01, file=paste0(trait,"_pip01.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip10) != 0) write.table(pip10, file=paste0(trait,"_pip10.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip25) != 0) write.table(pip25, file=paste0(trait,"_pip25.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip50) != 0) write.table(pip50, file=paste0(trait,"_pip50.dsv"), quote=F, row.names=F, sep="\t")
}

# POP LEVEL --------------------------------------------------------
setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/bslmm_poplevel/")
beak_traits <- c("feather_PC", "pen_num", "plum_num", "pen_length", "plum_length", "nodes")
for(trait in beak_traits) {
  params<-fread(paste0(trait,"_fam1.bslmm.param.txt"),header=T,sep="\t", data.table=F)
  params.pipsort<-params[order(-params$gamma),]
  pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
  # variants with effect in 10% MCMC samples or more
  pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
  # variants with effect in 25% MCMC samples or more
  pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
  # variants with effect in 50% MCMC samples or more
  pip50<-params.pipsort[params.pipsort$gamma>=0.50,]
  
  if (nrow(pip01) != 0) write.table(pip01, file=paste0(trait,"_pip01_fam1.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip10) != 0) write.table(pip10, file=paste0(trait,"_pip10_fam1.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip25) != 0) write.table(pip25, file=paste0(trait,"_pip25_fam1.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip50) != 0) write.table(pip50, file=paste0(trait,"_pip50_fam1.dsv"), quote=F, row.names=F, sep="\t")
}

setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/bslmm_poplevel/")
beak_traits <- c("feather_PC", "pen_num", "plum_num", "pen_length", "plum_length", "nodes")
for(trait in beak_traits) {
  params<-fread(paste0(trait,"_fam2.bslmm.param.txt"),header=T,sep="\t", data.table=F)
  params.pipsort<-params[order(-params$gamma),]
  pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
  # variants with effect in 10% MCMC samples or more
  pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
  # variants with effect in 25% MCMC samples or more
  pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
  # variants with effect in 50% MCMC samples or more
  pip50<-params.pipsort[params.pipsort$gamma>=0.50,]
  
  if (nrow(pip01) != 0) write.table(pip01, file=paste0(trait,"_pip01_fam2.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip10) != 0) write.table(pip10, file=paste0(trait,"_pip10_fam2.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip25) != 0) write.table(pip25, file=paste0(trait,"_pip25_fam2.dsv"), quote=F, row.names=F, sep="\t")
  if (nrow(pip50) != 0) write.table(pip50, file=paste0(trait,"_pip50_fam2.dsv"), quote=F, row.names=F, sep="\t")
}
