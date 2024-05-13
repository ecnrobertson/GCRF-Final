library(tidyverse)
library(ggplot2)
library(geiger)
library(ggrepel)
library(ggpubr) 
library(data.table)
library("ggnewscale") 
library('ggrastr')

setwd("~/Desktop/GCRF/GCRF-Final/figures/figures_data/")

width <-read_delim("feather_PC_fam2_ulmm_snps_ZFchr_assoc.txt",delim="\t") %>%  
  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width[order(as.numeric(width$ZFCHROM)),]
c <- c %>% mutate(win_num = 1:nrow(c)) #this is important, is basically takes each snp and gives it a separate order number that is relative to all the others instead of relative to the chromosome of scaffold it's on, so it's 1-# of snps. Without this it's pretty tricky to plot

c2 <- c[sample(nrow(c), size=200000, replace = F),] #I take a subsample of the data to actually plot initially because trying to get stuff to load with all 7mil snps just doesn't work, so you can figure out all the details with this sample as the background puoints (actual significant snps plotted separately) and then use the full version when saving as a pdf

frequency <- table(c$ZFCHROM)
start <- 1
lab.pos<-vector()
for (i in seq_along(frequency)){
  size <- frequency[i]
  txtpos <- start+size/2
  lab.pos <- c(lab.pos, txtpos)
  start <- start+size
}
#The above chunk of code is basically to set up where the labels will go, this allows the labels to be centered in the middle of the chromosome block

specific_breaks <- lab.pos

specific_labels <- c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", 
                     "8", "9", "10", "11", "12", "13", "14",
                     " ", "16", " ", " ", "20", " ",
                     " ", " ", " ", " ", " ", " ", "28",
                     "Z", "LGE22", " ") 
#Giving the chromosome labels, note that I leave some as blank because they get so cramped at the end that they overlap, so I removed some. The ticks will still be there but the labels will be readable.

#The next chunk of code is setting up a separate spreadsheet with just the snps that are significant. I found this to be the easiest way to play with colors, control what's rasterized, change shapes and especially at gene labels.
snps_sheet <- read.table("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/lmm_poplevel/feather_PC_fam2_outlier.txt", header=T)
snps_genes <- read.table("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/analysis_output_poplevel/genes/feather_PC_fam2_genes_lmm_100kb.txt", header=T)
snps <- snps_sheet$rs

sig_snps <- c[c$rs %in% snps, ]
#sig_snps <- c[match(snps, c$rs), ]
sig_snps$trait <- "LMM"
sig_snps <- left_join(sig_snps, snps_genes, by="rs")
sig_snps$annotated <- "no"
for (row in 1:nrow(sig_snps)) {
  sig_snps$annotated[row] <- ifelse(!is.na(sig_snps$gene_name[row]), "yes", sig_snps$annotated[row])
}


#Here's the actual plotting!
option1<-ggplot(c, aes(x = win_num, y = -log10(p_wald), color = as.factor(ZFCHROM))) + #make sure you set the chromosome values to a factor
  rasterise(geom_point(alpha = 1, size = 1, shape = 19), dpi=300) + #note that these are the background dots, and what should be rasterized to get anything to work
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", alpha = 0.8) + #adds dashed line for your significance threshold
  scale_color_manual(values = rep(c("grey80", "grey40", "grey20"), 200)) + #changing the colors of the chromosome blocks
  scale_x_continuous(breaks = specific_breaks, labels = specific_labels) + #add the chromosome labels
  scale_y_continuous(expand = c(0, 0), limits = c(0, 9)) + #this can get fiddled with to change the height of the graph, I like them to look a little stouter
  guides(color = "none") + #this keeps it from trying to create a legend with the chromosome info
  new_scale_color()+ #this allows the colors scales to coexist without removing the scale_color_manual above
  geom_point(data = subset(sig_snps, trait == "LMM"), aes(x = win_num, y = -log10(p_wald), color = "LMM"), size = 1.8) +
  geom_label_repel(data=subset(sig_snps,annotated=="yes"), aes(x=win_num, y=-log10(p_wald), label=gene_name), size=3,
                   segment.color = "black",segment.size = .2,
                   min.segment.length = unit(0, 'lines'),nudge_y = .03,max.overlaps = 50)+ #this creates little boxed with lines that have the gene names in them, with a repel to keep them from overlapping
  theme_classic() +
  labs(x = "Chromosome", y = "-log10(p-value)") +
  guides(color = guide_legend(override.aes = list(shape =c(19)))) + #changed the legend to have different shapes from the BSLMM and LMM models
  scale_color_manual(name = "Model",
                     values = c("LMM" = "#C69D44"),
                     breaks = c("LMM"),
                     labels = c("LMM")) #specifying the colors and what the legend will say

ggsave("~/Desktop/GCRF/GCRF-Final/figures/figures_output/manhattan/feather_PC_fam2_lmm.png", plot=option1,width=11.7, height=2.5, units="in")
png("~/Desktop/GCRF/GCRF-Final/figures/figures_output/manhattan/feather_PC_fam2_lmm.png", width=11.7, height=2.5, units="in", res=300)
print(option1)
dev.off()

option2<-ggplot(c2, aes(x = win_num, y = -log10(p_wald), color = as.factor(ZFCHROM))) + #make sure you set the chromosome values to a factor
  rasterise(geom_point(alpha = 1, size = 1, shape = 19), dpi=300) + #note that these are the background dots, and what should be rasterized to get anything to work
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", alpha = 0.8) + #adds dashed line for your significance threshold
  scale_color_manual(values = rep(c("grey80", "grey40", "grey20"), 200)) + #changing the colors of the chromosome blocks
  scale_x_continuous(breaks = specific_breaks, labels = specific_labels) + #add the chromosome labels
  scale_y_continuous(expand = c(0, 0), limits = c(0, 9)) + #this can get fiddled with to change the height of the graph, I like them to look a little stouter
  guides(color = "none") + #this keeps it from trying to create a legend with the chromosome info
  new_scale_color()+ #this allows the colors scales to coexist without removing the scale_color_manual above
  geom_point(data = subset(sig_snps, trait == "BSLMM"), aes(x = win_num, y = -log10(p_wald), color = "BSLMM"), size = 2.7, shape=17) +
  geom_point(data = subset(sig_snps, trait == "LMM"), aes(x = win_num, y = -log10(p_wald), color = "LMM"), size = 1.8) +
  geom_label_repel(data=subset(sig_snps,annotated=="yes"), aes(x=win_num, y=-log10(p_wald), label=gene), size=3,
                   segment.color = "black",segment.size = .2,
                   min.segment.length = unit(0, 'lines'),nudge_y = .03,max.overlaps = 50)+ #this creates little boxed with lines that have the gene names in them, with a repel to keep them from overlapping
  theme_classic() +
  labs(x = "Chromosome", y = "-log10(p-value)") +
  guides(color = guide_legend(override.aes = list(shape =c(17, 19)))) + #changed the legend to have different shapes from the BSLMM and LMM models
  scale_color_manual(name = "Model",
                     values = c("BSLMM" = "#A63C2E",
                                "LMM" = "#C69D44"),
                     breaks = c("BSLMM","LMM"),
                     labels = c("BSLMM","LMM")) #specifying the colors and what the legend will say
option2 #I only ever run this if the c2 (sample) data is being plotted, otherwise R just can't handle it and usually crashes

#quick export for powerpoints and stuff
png("~/Desktop/GCRF/GCRF-Present/files/results/figures/GCRF.feather_PC_labels.png", width=11.7, height=2.5, units="in", res=300)
ggarrange(option2)
dev.off()