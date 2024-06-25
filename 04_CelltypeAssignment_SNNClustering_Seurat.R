#setwd("~/Google Drive/sc_Aging/")
setwd("C:/Users/admm414r/Google Drive/EBI/sc_Aging/")
source("FabianScripts/functions.R")

library(dplyr)
library(Seurat)
library(stringr)
library(biomaRt)

load("RDataFiles/scMetaData.RDa")
load("RDataFiles/QCpassedCells.RDa")
load("RDataFiles/scMetaData.Rda")
load("RDataFiles/CorrrelationAssignmentResult_v2.RDa")
load("RDataFiles/Reference_Data_fake_SCs_V2.RDa")

###create seurat object, consider all genes which are expressed in at least 1% of the cells
snd <- CreateSeuratObject(counts = dat.counts.blood.QCpassed,
                          project = "single cells",
                          min.cells = 0.01*ncol(dat.counts.blood.QCpassed))
#assign mitochondrial genes
mt_genes <- subset(dat.features,Chr == "MT")
mt <- as.vector(mt_genes$Geneid)
mt <- mt[which(mt %in% rownames(snd))]
snd[["percent.mt"]] <- PercentageFeatureSet(snd,features = mt)
VlnPlot(snd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
        ,ncol = 3)


###normalize data
snd <- NormalizeData(snd, normalization.method = "LogNormalize", scale.factor = 10000)

###Add metadata
subs <- dat.info[rownames(dat.info)%in%colnames(snd),]

snd <- AddMetaData(object = snd,
                   metadata = subs$tissue,
                   col.name = "tissue")
snd <- AddMetaData(object = snd,
                   metadata = subs$mouse,
                   col.name = "mouse.number")
snd <- AddMetaData(object = snd,
                   metadata = subs$age_weeks,
                   col.name = "age.weeks")

##also add cell labels according to correlation assignment

snd <- AddMetaData(object = snd,
                   metadata = assigned.celltypes_byCorrelation$V3,
                   col.name = "correlation_lineage")

snd <- AddMetaData(object = snd,
                   metadata = assigned.celltypes_byCorrelation$V2,
                   col.name = "correlation_cell_type")

snd <- AddMetaData(object = snd,
                   metadata = assigned.celltypes_byCorrelation$bestLineage,
                   col.name = "correlation_best_lineage")

snd <- AddMetaData(object = snd,
                   metadata = assigned.celltypes_byCorrelation$bestCellType,
                   col.name = "correlation_best_cell_type")

###scale data
all.genes <- rownames(snd)
snd <- ScaleData(snd, features = all.genes)


markerGenes = unique(unlist(markers))
markerGenes = markerGenes[which(markerGenes %in% all.genes)]

###find top 1000 HVG
snd <- FindVariableFeatures(snd, selection.method = "vst", nfeatures = (length(markerGenes)))


###run pca
snd <- RunPCA(snd, features =union(markerGenes,VariableFeatures(object = snd)))
#snd <- RunPCA(snd, features =(markerGenes))
#snd <- RunPCA(snd, features = VariableFeatures(object = snd))
#plot pca
DimPlot(snd, reduction = "pca",#group.by = "age.weeks",
        dims = c(1,2))

###perform SNN clustering
snd <- FindNeighbors(snd, dims = 1:10)
snd <- FindClusters(snd, resolution = 0.8)

### UMAP
snd <- RunUMAP(snd, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

#filename <- "~/Desktop/UMAP_SeuratClusters.pdf"
#pdf(file = filename, width=6.5, height=5)
DimPlot(snd, reduction = "umap",#,
        label = T) #,
#dev.off()

#filename <- "~/Desktop/PCA_SeuratClusters.pdf"
#pdf(file = filename, width=6.5, height=5)
DimPlot(snd, reduction = "pca",#group.by = "age.weeks",
        dims = c(1,2))
#dev.off()


######################Identify Marker Genes per Cluster#######################

#snd.markers <- FindAllMarkers(snd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# ###assign gene names according to ENSEMBL IDs
# snd.markers <- assign.GeneNames(EnsemblIDs = snd.markers$gene,
#                                 target.DF = snd.markers)

###subset top 3 marker genes for every cluster according to adjusted p value
#top3 <- snd.markers %>% group_by(cluster) %>% top_n(n = 3, wt = -(p_val_adj))



######Color Code UMAP by Corrlation Assignment 
#filename <- "~/Desktop/UMAP_CorrAss.pdf"
#pdf(file = filename, width=6.5, height=5)
DimPlot(snd, reduction = "umap",group.by = "correlation_lineage",
        label = T) #,

DimPlot(snd, reduction = "umap",group.by = "correlation_best_lineage",
        label = T) #,


#dev.off()

View(cbind(snd$seurat_clusters, snd$correlation_best_lineage))
all(rownames(snd$seurat_clusters) == rownames(snd$correlation_cell_type))
#snd$nCount_RNA[1]
#snd$nFeature_RNA

##Step one lineage by correlation. 
table(paste(snd$seurat_clusters,snd$correlation_best_lineage,sep="-"))
lineageCleaned = snd$correlation_lineage
lineageCleaned[snd$seurat_clusters==0] = "B cell Lineage"
lineageCleaned[snd$seurat_clusters==1] = "B cell Lineage"
lineageCleaned[snd$seurat_clusters==2] = "T cell Lineage"
lineageCleaned[snd$seurat_clusters==3] = "T cell Lineage"
lineageCleaned[snd$seurat_clusters==4] = "T cell Lineage"
lineageCleaned[snd$seurat_clusters==5] = "B cell Lineage"
lineageCleaned[snd$seurat_clusters==6] = "Mixed"
lineageCleaned[snd$seurat_clusters==7] = "Mixed2"
lineageCleaned[snd$seurat_clusters==8] = "NK cell Lineage"


snd <- AddMetaData(object = snd,
                   metadata = lineageCleaned,
                   col.name = "Combined_cell_type_classification")

snd <- SetIdent(snd, value = "Combined_cell_type_classification")

snd.markers <- FindAllMarkers(snd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DimPlot(snd, reduction = "umap",group.by = "Combined_cell_type_classification",
        label = T) #,


##Step 2. Check cluster 6 & cluster 7, what are they?
snd_Mixed1 = snd[,which(snd$seurat_clusters==6)]
snd_Mixed1 <- FindVariableFeatures(snd_Mixed1, selection.method = "vst", nfeatures = (1000))
snd.markers[which(snd.markers$cluster=="Mixed"),]

for(n in names(markers)){
  print(paste(n,length(which(markers[[n]] %in% snd.markers$gene[which(snd.markers$cluster=="Mixed")]))))
}
###
#Mixed seem to be Bcells still.
lineageCleaned[snd$seurat_clusters==6] = "B cell Lineage"

###
snd_Mixed2 = snd[,which(snd$seurat_clusters==7)]
snd_Mixed2 <- FindVariableFeatures(snd_Mixed1, selection.method = "vst", nfeatures = (1000))
#snd.markers[which(snd.markers$cluster=="Mixed2"),]

for(n in names(markers)){
  print(paste(n,length(which(markers[[n]] %in% snd.markers$gene[which(snd.markers$cluster=="Mixed2")]))))
}

#Mixed2 seem to be a mixed bag of the other cell types mostly: Neutrophils & Macrophages.

# markerGenesMixedcells = unique(unlist(c(markers["Neutrophil Lineage"],markers["Macrophage Lineage"],markers["Eosinophil Lineage"],markers["Erythrocyte Lineage"],markers["Basophil Lineage"],markers["Dendritic Cell Lineage"])))
# markerGenesMixedcells = markerGenesMixedcells[which(markerGenesMixedcells %in% rownames(snd_Mixed2))]
# 
# ###find top 1000 HVG
# snd_Mixed2 <- FindVariableFeatures(snd_Mixed2, selection.method = "vst", nfeatures = (length(markerGenesMixedcells)))
# 
# ###run pca
# snd_Mixed2 <- RunPCA(snd_Mixed2, features =union(markerGenes,VariableFeatures(object = snd_Mixed2)))
# #snd_Mixed2 <- RunPCA(snd_Mixed2, features =(markerGenes))
# #snd_Mixed2 <- RunPCA(snd_Mixed2, features = VariableFeatures(object = snd_Mixed2))
# 
# ###perform SNN clustering
# snd_Mixed2 <- FindNeighbors(snd_Mixed2, dims = 1:10)
# snd_Mixed2 <- FindClusters(snd_Mixed2, resolution = 0.8)
# 
# ### UMAP
# snd_Mixed2 <- RunUMAP(snd_Mixed2, dims = 1:10)
# 
# # note that you can set `label = TRUE` or use the LabelClusters function to help label
# # individual clusters
# 
# #filename <- "~/Desktop/UMAP_SeuratClusters.pdf"
# #pdf(file = filename, width=6.5, height=5)
# DimPlot(snd_Mixed2, reduction = "umap",#,
#         label = T) #,
# table(paste(snd_Mixed2$seurat_clusters,snd_Mixed2$correlation_lineage,sep="-"))

snd <- AddMetaData(object = snd,
                   metadata = lineageCleaned,
                   col.name = "Combined_cell_type_classification")

snd <- SetIdent(snd, value = "Combined_cell_type_classification")
DimPlot(snd, reduction = "umap",group.by = "Combined_cell_type_classification",
        label = T) #,

write.table(cbind((snd$age.weeks),(snd$mouse.number),(snd$Combined_cell_type_classification)),"./RDataFiles/broadClassification.txt",sep="\t",quote=F)

##Step 3. T-cells.
snd_Tcells = snd[,which(snd$Combined_cell_type_classification=="T cell Lineage")]

assigned.celltypes_byCorrelation_Tcells = assigned.celltypes_byCorrelation[which(rownames(assigned.celltypes_byCorrelation) %in% colnames(snd_Tcells)), ]
assigned.celltypes_byCorrelation_Tcells = assigned.celltypes_byCorrelation_Tcells[,c(20,21,34,35,37,38,40:42,44)]
Best_tCell = NULL
for(i in 1:nrow(assigned.celltypes_byCorrelation_Tcells)){
   Best_tCell = c(Best_tCell,colnames(assigned.celltypes_byCorrelation_Tcells)[which(assigned.celltypes_byCorrelation_Tcells[i,] ==max(assigned.celltypes_byCorrelation_Tcells[i,]))])       
}
assigned.celltypes_byCorrelation_Tcells = cbind(assigned.celltypes_byCorrelation_Tcells, Best_tCell)
all(colnames(snd_Tcells) == rownames(assigned.celltypes_byCorrelation_Tcells))

snd_Tcells <- AddMetaData(object = snd_Tcells,
                   metadata = assigned.celltypes_byCorrelation_Tcells$Best_tCell,
                   col.name = "correlation_best_cell_type_T")


###scale data
snd_Tcells <- ScaleData(snd_Tcells, features = all.genes)
topMarkersT = readRDS(file = "./RDataFiles/Tmarkers.Rds")

markerGenesTcells = topMarkersT$gene[which(unique(topMarkersT$gene) %in% rownames(snd_Tcells))]

snd_Tcells <- FindVariableFeatures(snd_Tcells, selection.method = "vst", nfeatures = length(markerGenesTcells))

###run pca
snd_Tcells <- RunPCA(snd_Tcells, features =union(markerGenesTcells,VariableFeatures(object = snd_Tcells)))
#snd_Tcells <- RunPCA(snd_Tcells, features =markerGenesTcells)
#snd_Tcells <- RunPCA(snd_Tcells, features =VariableFeatures(object = snd_Tcells))
#plot pca
DimPlot(snd_Tcells, reduction = "pca",#group.by = "age.weeks",
        dims = c(1,2))

maxDim = 5
snd_Tcells <- FindNeighbors(snd_Tcells, dims = 1:maxDim)
snd_Tcells <- FindClusters(snd_Tcells, resolution = 0.5)
snd_Tcells.markers <- FindAllMarkers(snd_Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
count(snd_Tcells)
snd_Tcells <- RunUMAP(snd_Tcells, dims = 1:maxDim)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

#filename <- "~/Desktop/UMAP_SeuratClusters.pdf"
#pdf(file = filename, width=6.5, height=5)
DimPlot(snd_Tcells, reduction = "umap",#,
        label = T) #,

##Select from correlations only the T-cells
DimPlot(snd_Tcells, reduction = "umap",group.by = "correlation_best_cell_type_T",
        label = T) #,

table(paste(snd_Tcells$seurat_clusters,snd_Tcells$correlation_best_cell_type_T,sep="-"))

tMarkers = FindAllMarkers(snd_Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
for(c in unique(tMarkers$cluster)){
        print(c)
        for(n in unique(topMarkersT$cluster)){
          print(paste("   ",n,length(which(topMarkersT$gene[which(topMarkersT$cluster==n)] %in% tMarkers$gene[which(tMarkers$cluster==c)]))))
        }
}

lineageCleaned_T = as.character(snd_Tcells$seurat_clusters)
lineageCleaned_T[snd_Tcells$seurat_clusters==0] = "Naive T cell"
lineageCleaned_T[snd_Tcells$seurat_clusters==1] = "Mem T cell"
lineageCleaned_T[snd_Tcells$seurat_clusters==2] = "Eff T cell"
lineageCleaned_T[snd_Tcells$seurat_clusters==3] = "Mem T cell"

snd_Tcells <- AddMetaData(object = snd_Tcells,
                   metadata = lineageCleaned_T,
                   col.name = "correlation_Tcell_lineage")

DimPlot(snd_Tcells, reduction = "umap",group.by = "correlation_Tcell_lineage",
        label = T)

tCellExp = as.data.frame(GetAssayData(snd_Tcells, slot="data"))
cd4.8.Exp = rbind(tCellExp[which(rownames(tCellExp)=="ENSMUSG00000023274"),],tCellExp[which(rownames(tCellExp)=="ENSMUSG00000053977"),])

cdStatus = rep("CD8-CD4-",ncol(cd4.8.Exp))
cdStatus[which(cd4.8.Exp[1,]>cd4.8.Exp[2,])] = "CD4+"
cdStatus[which(cd4.8.Exp[2,]>cd4.8.Exp[1,])] = "cD8+"

cbind(cdStatus,lineageCleaned_T)

write.table(cbind(colnames(snd_Tcells),lineageCleaned_T,cdStatus),"./RDataFiles/Tcell_Classification.txt",sep="\t",quote=F,row.names=F)


# dat.info$nCount_RNA <- snd$nCount_RNA[match(rownames(dat.info),colnames(snd))]
# dat.info$nFeature_RNA <- snd$nFeature_RNA[match(rownames(dat.info),colnames(snd))]
# dat.info$Seurat_cluster_resolution0.5 <- snd$seurat_clusters[match(rownames(dat.info),colnames(snd))]
# 
# save(dat.ages,dat.info,dat.features,file = "RDataFiles/scMetaData_v2.Rda")
