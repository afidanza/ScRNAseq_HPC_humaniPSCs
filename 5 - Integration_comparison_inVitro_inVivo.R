library(ggplot2)
library(Seurat)
library(dplyr)
library(plyr)
library(Matrix)
library(tibble)
library(biomaRt)
library(reshape2)

seurat.list <- c(InVitroHSCMPP, InVivoHSCMPP)
rm (InVitroHSCMPP, InVivoHSCMPP)
integration.features <- SelectIntegrationFeatures(object.list = seurat.list, 
                                                  nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list,
                                  anchor.features = integration.features,
                                  verbose = FALSE)
HSC.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                      normalization.method = "SCT", 
                                      anchor.features = integration.features, 
                                      verbose = FALSE)
Integrated_HSCMPP <- IntegrateData(anchorset = HSC.anchors, 
                                   normalization.method = "SCT", 
                                   verbose = FALSE)

rm(seurat.list)
rm(HSC.anchors)

########## DIMENSION REDUCTION
Integrated_HSCMPP <- RunPCA(Integrated_HSCMPP, verbose = FALSE)
Integrated_HSCMPP <- RunTSNE(Integrated_HSCMPP, dims = 1:15)
Integrated_HSCMPP <- RunUMAP(Integrated_HSCMPP, dims = 1:15)
Integrated_HSCMPP <- FindNeighbors(object = Integrated_HSCMPP, 
                                   dims = 1:15, verbose = FALSE)
Integrated_HSCMPP <- FindClusters(object = Integrated_HSCMPP, 
                                  resolution = 0.1, verbose = FALSE)

###### DEG analysis

#store cluster identity in another metadata
Integrated_HSCMPP$CellClusters <- Idents(Integrated_HSCMPP)
#rename origin identities
Integrated_HSCMPP$orig.ident<- revalue(Integrated_HSCMPP$orig.ident, c ("10x_2019"= "in vitro", "Embryo-liver"= "in vivo"))
#assign origin of cells to identities
Idents(Integrated_HSCMPP) <- Integrated_HSCMPP$orig.ident
#assign clusters to identities
Idents(Integrated_HSCMPP) <- Integrated_HSCMPP$CellClusters
#calculate DEGs
DEGinvitro<- FindMarkers(Integrated_HSCMPP, 
                         ident.1 = "10x_2019", ident.2 = "Embryo-liver",
                         min.pct = 0.25, logfc.threshold = 0.25)

#order by pval
DEGinvitro <- DEGinvitro[order(DEGinvitro$p_val_adj),]
DEGinvitro <- DEGinvitro[DEGinvitro$p_val_adj < 0.05,]

###filter TFs
genelist <- DEGinvitro
names(genelist)[1] <- "gene"
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
go_id_TF = c('GO:0003700', 'GO:0000130', 'GO:0001071')
genelist$gene <- rownames(genelist)
filteredgenes <- getBM(attributes = c('hgnc_symbol'),
                       filters = c('go', 'hgnc_symbol'), 
                       values = list(go_id_TF, genelist$gene), 
                       mart = ensembl, uniqueRows = T)
TF.ids <- genelist$gene %in% filteredgenes$hgnc_symbol
genelist.TF <- genelist[TF.ids,]
genelist.TF <- genelist.TF[order(genelist.TF$avg_logFC, decreasing = TRUE),]

