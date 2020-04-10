library(Seurat)
library(ggplot2)

# Read 10X data
data10x <- Read10X("outs/raw_feature_bc_matrix/")
# Create seurat object from the gene expression data only
seuratObj <- CreateSeuratObject(data10x$`Gene Expression`,
                                min.cells = 3,
                                min.features = 200,
                                project = "10x 2019")
# Create the ADT assay and add to Seurat object
adt <- CreateAssayObject(data10x$`Antibody Capture`)
adt <- SubsetData(adt, cells = Cells(seuratObj))
seuratObj[["ADT"]] <- adt
rm(data10x)

# Add metadata
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj[["libraryID"]] <- substr(colnames(seuratObj), 18, 18)
seuratObj[["Day"]] <- ifelse(seuratObj[["libraryID"]] %in% c("1", "2", "3"), 
                             "D10", "D13")

# QC plots
VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01) 

ggplot(seuratObj@meta.data, aes(nFeature_RNA, percent.mt)) +
  geom_point(aes(col = Day), size = 0.5) +
  geom_hline(yintercept = c(4,20)) +
  geom_vline(xintercept = c(2000, 7500)) 

boxplot(percent.mt ~ Day, seuratObj@meta.data)

table(seuratObj@meta.data$percent.mt > 5, seuratObj@meta.data$Day)

FeatureScatter(seuratObj, "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  geom_hline(yintercept = c(4,20)) +
  geom_vline(xintercept = c(2000, 7500)) 

# Filtering
seuratObj <- subset(seuratObj, 
                    subset = percent.mt > 4 & percent.mt < 20 & nFeature_RNA > 2000 & 7500)
# SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
table(seuratObj@meta.data$Day)
seuratObj <- SCTransform(seuratObj, vars.to.regress = "percent.mt")

seuratObj <- RunPCA(object = seuratObj)
ElbowPlot(seuratObj, ndims = 50)
seuratObj <- RunUMAP(object = seuratObj, dims = 1:15)
seuratObj <- FindNeighbors(object = seuratObj, dims = 1:15, verbose = FALSE)
seuratObj <- FindClusters(object = seuratObj, verbose = FALSE)
DimPlot(object = seuratObj, label = TRUE)