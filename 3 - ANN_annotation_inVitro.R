library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
library(biomaRt)
library(reshape2)

ANNident <- read.csv("predconsensus.txt", header = F)
SeuratEmatoEndo10_13 <- readRDS("SeuratObject.rds")
#assign ANN identit to a new metadata 
SeuratEmatoEndo10_13$ANN <- ANNident$V1
#store clustering identities into a new metadata so it can be recalled later
SeuratEmatoEndo10_13$cluster_idents <- Idents(SeuratEmatoEndo10_13)
#assign ANN identities to the idents
Idents(SeuratEmatoEndo10_13) <- SeuratEmatoEndo10_13$ANN
DimPlot(SeuratEmatoEndo10_13)
