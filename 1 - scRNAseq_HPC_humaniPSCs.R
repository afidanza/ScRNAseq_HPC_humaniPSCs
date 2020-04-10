library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
library(biomaRt)

######################################################
# CREATE THE SEURAT OBJECT FROM THE SEQUENCING MATRIX#
######################################################

# Change to directory containing the 10x data files (genes.tsv, barcodes.tsv, matrix.mtx)
setwd("~/scRNAseq")
data10x <- Read10X(".")
seuratobject <- CreateSeuratObject(raw.data = data10x, 
                                   project = "scRNAseq_diff_hIPSC",
                                   min.cells = 3, # Minimum 3 cells/gene
                                   min.genes = 200) # Minimum 200 genes/cell
rm (data10x)

##########################
#QUALITY CONTROL ANALYSIS#
##########################

# Names of mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratobject@data), value = TRUE)
# Proportion of mitochondrial genes
percent.mito <- Matrix::colSums(seuratobject@raw.data[mito.genes, ])/Matrix::colSums(seuratobject@raw.data)
# Add the proportion to seurat metadata
seuratobject <- AddMetaData(object = seuratobject, metadata = percent.mito, col.name = "percent.mito")

# Diagnostic plots
VlnPlot(object = seuratobject, features.plot = c("percent.mito", "nUMI", "nGene"), nCol = 3, point.size.use = 0.01)

par(mfrow = c(1, 3))
GenePlot(object = seuratobject, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.1)
GenePlot(object = seuratobject, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.1)
GenePlot(object = seuratobject, gene1 = "nGene", gene2 = "percent.mito", cex.use = 0.1)

# Filtering limits shown in red
par(mfrow = c(1, 1))
GenePlot(object = seuratobject, ylim=c(0, 0.4), gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.1)
abline(v=5000, h=0.005, col= 'red', lwd = 2)
abline(v=50000, h=0.05, col= 'green', lwd = 2)

# Filter cells according to limits above for % mitochondrial genes and number of UMI
seuratobject <- FilterCells(object = seuratobject, 
                            subset.names = c("percent.mito","nUMI"),
                            low.thresholds =  c(0.005, 5000), 
                            high.thresholds = c(0.05, 50000))

# Normalise data
seuratobject <- NormalizeData(object = seuratobject, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
# Find most variable genes
seuratobject <- FindVariableGenes(object = seuratobject, 
                                  mean.function = ExpMean, 
                                  dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, 
                                  x.high.cutoff = 3, 
                                  y.cutoff = 0.5)

###############################################################
# Scaling the data and removing unwanted sources of variation #
###############################################################

# Find batch number for each cell based on the 10X library tag, which is character 18 on the tag 
batch <- as.integer (substring(seuratobject@cell.names, 18))
# Add batch information to metadata
seuratobject@meta.data$batch = batch
# Regress out any possible UMI and batch effects
seuratobject <- ScaleData(object = seuratobject,
                          vars.to.regress = c("nUMI", "batch"))

###############################################
# PERFORM LINEAR DIMENSION REDUCTION WITH PCA #
###############################################

# Run the PCA analysis, using the most variable genes
seuratobject <- RunPCA(object = seuratobject, 
                       pc.genes = seuratobject@var.genes)

# Elbow plot to choose the number of PCs to use
PCElbowPlot(object = seuratobject)

#######################
# CLUSTERING and tSNE #
#######################

# Cluster the cells
seuratobject <- FindClusters(object = seuratobject, 
                             reduction.type = "pca", 
                             dims.use = 1:15, 
                             resolution = 0.6, 
                             print.output = 0, 
                             save.SNN = TRUE, 
                             force.recalc = TRUE)

# Run the tSNE analysis and plot its result
seuratobject <- RunTSNE(object = seuratobject, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = seuratobject, pt.size = 0.5)
# Identify the markers for each cluster
seuratobject.markers <- FindAllMarkers(object = seuratobject, 
                                       only.pos = TRUE, 
                                       min.pct = 0.25)

############################################
# DETERMINATION OF MEMBRANE AND TF MARKERS #
############################################

for (clusternumber in unique(seuratobject@ident))
  {
  clustermarkers <- seuratobject.markers [seuratobject.markers$cluster == clusternumber, 
                                  c("gene","avg_logFC", "p_val_adj")]
  write.csv (clustermarkers, 
             file = paste0 ("markers/cluster_markers_", clusternumber, ".csv"), 
             row.names = FALSE)
  }

######## Filtering TF and membrane markers ############
for (clusternumber in unique(seuratobject@ident)) 
  {
  genelist <- read.csv(file = paste0("markers/cluster_markers_", clusternumber, ".csv"),
                       stringsAsFactors = FALSE)

  # We are using the Homo Sapiens mart
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  # GO IDS we are interested in
  go_id_TF = c('GO:0003700', 'GO:0000130', 'GO:0001071')
  go_id_membr = c('GO:0005886')
  #go_id_receptor = 'GO:0098802'
  filteredgenes <- getBM(attributes = c('hgnc_symbol'),
                         filters = c('go', 'hgnc_symbol'),
                         values = list(go_id_TF, genelist$gene),
                         mart = ensembl, uniqueRows = T)

  # We can create a subset of values of the original gene list
  # Find the lines in the original list of genes which contain the result of our query
  # %in% returns for each position in genelist$gene
  # TRUE of FALSE if they are or not present in filterdgenes$hgnc-symbol
  TF.ids <- genelist$gene %in% filteredgenes$hgnc_symbol
  # Just get those lines
  genelist.TF <- genelist[TF.ids,]

  # reassign to genelist.TF the values ordered by logFC
  genelist.TF <- genelist.TF[order(genelist.TF$avg_logFC, decreasing = TRUE),]

  # filter for log FC >=1
  # genelist.TF <- genelist.TF [(genelist.TF$avg_logFC >= 1),]

  write.csv(genelist.TF, file = paste0("TF_cluster_", clusternumber, ".csv"), row.names = FALSE)

  barplot ( rev(genelist.TF$avg_logFC),
            names.arg = rev(genelist.TF$gene),
            las = 1,
            horiz = TRUE,
            xlim = c(0,1.5),
            col = "black",
            main = paste0 ("TF Cluster ", clusternumber),
            cex.names = 1,
            cex.axes = 1)

  filteredgenes <- getBM(attributes = c('hgnc_symbol'),
                         filters = c('go', 'hgnc_symbol'),
                         values = list(go_id_membr, genelist$gene),
                         mart = ensembl, uniqueRows = T)

  membr.ids <- genelist$gene %in% filteredgenes$hgnc_symbol
  # Just get those lines
  genelist.membr <- genelist[membr.ids,]

  genelist.membr <- genelist.membr[order(genelist.membr$avg_logFC, decreasing = TRUE),]

  write.csv ( genelist.membr, file = paste0("Membrane_cluster_", clusternumber, ".csv"), row.names = FALSE)

  barplot ( rev(genelist.membr$avg_logFC),
            names.arg = rev(genelist.membr$gene),
            las = 1,
            horiz = TRUE,
            xlim = c(0,4),
            col = "black",
            main = paste0("Membrane Marker Cluster ", clusternumber),
            cex.names = 0.3,
            cex.axes = 1)
  }

# Removal of doublet/mixed cluster (number 8)
seuratobject <- SubsetData(seuratobject, ident.remove = 8)
# Remove markers of cluster 8
seuratobject.markers <- seuratobject.markers[seuratobject.markers$cluster != 8,]

# Save object to file for later use
# save(seuratobject, file = "seurat.RData")

#################################
# TRAJECTORY ANALYSIS - MONOCLE #
#################################

monocleobject <- importCDS(seuratobject, import_all=T)

# Calculate size factors and dispersions
monocleobject <- estimateSizeFactors(monocleobject)
monocleobject <- estimateDispersions(monocleobject)
monocleobject <- detectGenes(monocleobject, min_expr = 0.1)
# Filters for genes that are expressed in at least 50 cells,
expressed_genes <- row.names(subset(fData(monocleobject), num_cells_expressed >= 50))
  
# Assign tSNE coordinates from Seurat object
monocleobject@reducedDimA <- t(seuratobject@dr$tsne@cell.embeddings)
# Copy phenodata from that imported from Seurat
monocleobject@phenoData@data$Cluster <- pData(monocleobject)$res.0.6

plot_cell_clusters(monocleobject, 1, 2, color_by = "Cluster", cell_size = 3)

# Get the top 20 Seurat markers and use them as ordering genes for Monocle
monocleobject_ordering_genes <- seuratobject.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  
monocleobject <- setOrderingFilter(monocleobject,
                                    ordering_genes = monocleobject_ordering_genes)
  
# Plot the elbow plot of the PC contribution to the variance
plot_pc_variance_explained(monocleobject, return_all = FALSE)
  
# Run dimensionality reduction
monocleobject <- reduceDimension (monocleobject,
                                   method = 'DDRTree',
                                   max_components = 6)
  
# Order the cells
monocleobject <- orderCells(monocleobject, num_paths = 3)
  
plot_cell_trajectory(monocleobject, color_by = "Cluster", x=1, y=2)

# Save object to file for later use
# save(monocleobject, file = "monocle.RData")
  
##########################
# DIFFUSION MAP ANALYSIS #
##########################

seuratobject <- RunDiffusion(seuratobject)
# Plot results
DMPlot(seuratobject)
DimPlot(seuratobject, reduction.use = "dm", do.return = T)