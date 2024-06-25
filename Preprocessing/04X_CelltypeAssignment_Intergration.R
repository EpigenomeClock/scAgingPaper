#setwd("~/Google Drive/sc_Aging/")
setwd("C:/Users/admm414r/Google Drive/EBI/sc_Aging/")
source("FabianScripts/functions.R")

library(dplyr)
library(Seurat)
library(scran)
library(stringr)
library(biomaRt)

load("RDataFiles/scMetaData.RDa")
load("RDataFiles/QCpassedCells.RDa")
load("RDataFiles/scMetaData.Rda")
load("RDataFiles/CorrrelationAssignmentResult_v2.RDa")
load("RDataFiles/Reference_Data_fake_SCs_V2.RDa")

dat.counts.blood.QCpassed.cpm <- correctedCPM(dat.counts.blood.QCpassed)
merged.simulated.normalized.counts = ((exp(1)^(merged.simulated.logNormalized.counts))-1)


###create seurat object, consider all genes which are expressed in at least 1% of the cells
snd <- CreateSeuratObject(counts = cbind(dat.counts.blood.QCpassed.cpm,merged.simulated.normalized.counts),
                          project = "single cells",
                          min.cells = 0.01*ncol(dat.counts.blood.QCpassed))

###normalize data
snd <- NormalizeData(snd, normalization.method = "LogNormalize", scale.factor = 10000)

###scale data
all.genes <- rownames(snd)
snd <- ScaleData(snd, features = all.genes)

markerGenes = unique(unlist(markers)) 
markerGenes = markerGenes[which(markerGenes %in% all.genes)]

###run pca on marker genes.
snd <- RunPCA(snd, features =markerGenes)
#snd <- RunPCA(snd, features = VariableFeatures(object = snd))
#plot pca
DimPlot(snd, reduction = "pca",#group.by = "age.weeks",
        dims = c(1,2))

###perform SNN clustering
snd <- FindNeighbors(snd, dims = 1:10)
snd <- FindClusters(snd, resolution = 1.0)

### UMAP
snd <- RunUMAP(snd, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

#filename <- "~/Desktop/UMAP_SeuratClusters.pdf"
#pdf(file = filename, width=6.5, height=5)
DimPlot(snd, reduction = "umap",#,
        label = T) #,
#dev.off()

