library(Seurat)
library(ggplot2)

# In vivo data was downloaded from 
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/
# The E-MTAB-7407.sdrf.txt file contains information about the samples
samples <- read.table("E-MTAB-7407.sdrf.txt", sep = "\t", head = TRUE, 
                      stringsAsFactors = FALSE)

# Filter our samples
# Some data is not available as processed, we remove that
available.dirs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
samples <- samples[which(samples$Source.Name %in% available.dirs),]
nrow(samples)

proj.name <- "Embryo-liver"
samples <- samples[samples$Factor.Value.organism.part. == "liver",]

# Read annotated identities
idents.list <- sapply(1:nrow(samples), function(i)
{
  s <- samples$Source.Name[i]
  idpath <- (paste0(s, "/", s, ".csv"))
  
  print(paste("Reading", idpath))
  ids <- read.csv(idpath, stringsAsFactors = F)
  ids$Barcodes <- paste0(ids$Barcodes, "_", i)
  ids$Cell.Labels <- as.factor(ids$Cell.Labels)
  
  ids
}, simplify = FALSE)

# Loop through the samples and create seurat objects
i <- 1
seur.obj.list <- apply(samples, 1, function(s)
  {
  print(paste("Processing sample", i, " - ", s["Source.Name"]))
  path <- (paste0(s["Source.Name"], "/GRCh38"))
  data10X <- Read10X(path)
  # Convert into a Seurat object + initial filtering
  seuratObject <- CreateSeuratObject(data10X,
                                     min.cells = 3, # min 3 cells/gene
                                     min.features = 200, # min 200 genes/cell
                                     project = proj.name)
  print(paste("Read", length(Cells(seuratObject)), "cells"))
  # Save metadata
  seuratObject <- AddMetaData(seuratObject, s["Characteristics.sex."], "Sex")
  wk <- as.integer(substr(s["Characteristics.age."], 1, 2))
  seuratObject <- AddMetaData(seuratObject, wk, "Weeks.gestation")
  seuratObject <- AddMetaData(seuratObject, s["Characteristics.clinical.information."], "Trimester")
  seuratObject <- AddMetaData(seuratObject, s["Factor.Value.facs.sorting."], "FACS.sorting")
  # Assign cell labels
  Idents(seuratObject) <- idents.list[[i]]$Cell.Labels

  i <<- i + 1

  seuratObject
  })

# Save as RDS object to avoid to recalculate later on
saveRDS(seur.obj.list, "ListSeuratObjLiver.rds")

# Merge the objects in the list
seuratObject <- seur.obj.list[[1]]
# Append _1 to cell names
seuratObject <- RenameCells(seuratObject, new.names = paste0(Cells(seuratObject), "_1"))

# We loop through the samples and merge them one by one
# merge could also be run directly on the list once, but alot of memory is required
for (i in 1:length(seur.obj.list))
  {
  print(paste("Adding sample", i, "-", length(Cells(seur.obj.list[[2]])), "cells"))
  # Add samples one by one and delete them from the list to save memory
  seur.obj.list[[2]] <- RenameCells(seur.obj.list[[2]], new.names = paste0(Cells(seur.obj.list[[2]]), "_", i))

  seuratObject <- merge(seuratObject, seur.obj.list[[2]])
  seur.obj.list[[2]] <- NULL
  }

print(paste("Merging complete.", length(Cells(seuratObject)), "processed."))

# store mitochondrial percentage in object meta data
seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
batch <- substr(Cells(seuratObject), 18, 20)
seuratObject <- AddMetaData(seuratObject, batch, "batch")
# run sctransform
seuratObject <- SCTransform(seuratObject, vars.to.regress = c("batch", "percent.mt"))
seuratObject <- RunPCA(seuratObject)
seuratObject <- RunTSNE(seuratObject, dims = 1:20)
seuratObject <- RunUMAP(seuratObject, dims = 1:20)

DimPlot(seuratObject) + 
  facet_wrap(factor(substr(Cells(seuratObject), 18, 20))) +
  NoLegend()

saveRDS(seuratObject, "SeuratLiverSamples.rds")