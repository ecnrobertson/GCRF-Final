library(data.table)
library(tidyverse)
library(dplyr)

setwd("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/genes/")
lmm_gene_list <- list.files(pattern="lmm_100kb.txt", full.names = TRUE)
gene_names <- lapply(lmm_gene_list, function(file) read.table(file, fill = TRUE, header=T, sep="\t"))

BCRF_genes <- read.table("~/Desktop/GCRF/GCRF-Present/BCRF/BCRF_gene.csv", header=T,sep = ",")
BCRF_genes <- BCRF_genes[,2]

#This loops through the lmm_results and finds overlapping genes! Very neat
for (file in seq_along(gene_names)) {
  x <- gene_names[[file]]
  gene_names_col <- x[, 3]
  for (row in 1:nrow(x)) {
    gene <- gene_names_col[row]
    results <- grep(pattern = gene, x = BCRF_genes)
    if (length(results) == 0) {
      cat("Gene", gene, "not found in BCRF_genes.\n")
      next  # Skip to the next iteration if gene is not found
    }
    print(gene)
  }
}

#Same idea but for the bslmm
bslmm_gene_list <- list.files(pattern="bslmm_100kb.txt", full.names = TRUE)
bslmm_gene_names <- lapply(bslmm_gene_list, function(file) read.table(file, fill = TRUE, header=T, sep="\t"))

for (file in seq_along(bslmm_gene_names)) {
  x <- bslmm_gene_names[[file]]
  gene_names_col <- x[, 3]
  for (row in 1:nrow(x)) {
    gene <- gene_names_col[row]
    results <- grep(pattern = gene, x = BCRF_genes)
    if (length(results) == 0) {
      cat("Gene", gene, "not found in BCRF_genes.\n")
      next  # Skip to the next iteration if gene is not found
    }
    print(gene)
  }
}

master_snps <- read.csv("~/Desktop/GCRF/GCRF-Final/analysis/analysis_output/results_master.csv")
gene_names_col <- master_snps$gene_name
for (row in 1:nrow(master_snps)) {
    gene <- gene_names_col[row]
    results <- grep(pattern = gene, x = BCRF_genes)
    if (length(results) == 0) {
      cat("Gene", gene, "not found in BCRF_genes.\n")
      next  # Skip to the next iteration if gene is not found
    }
    print(gene)
  }
}

grep("ADGB", BCRF_genes)
