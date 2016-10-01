library(shiny)
library(readr)
library(dplyr)
library(ReactomePA)
library(visNetwork)
library(pathview)
library(png)
#library(org.Hs.eg.db)
probe_list <- read_csv("PROBELIST.csv")
gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(na.omit(gene_names))
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names
#if error, check for nas


#check out seurat package
