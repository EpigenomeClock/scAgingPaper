#setwd("~/Google Drive/sc_Aging/")
setwd("G:/My Drive/EBI/sc_Aging/")

source("./Fabian/FabianScripts/functions.R")

######Assign single cells to a cell type based on the reference data 

library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(scater)
library(scran)

load("RDataFiles/scMetaData.RDa")
#load("RDataFiles/Reference_Data_new.RDa")
load("RDataFiles/Reference_Data_V2.RDa")
load("RDataFiles/Reference_Data_fake_SCs_V2.RDa")
load("RDataFiles/QCpassedCells.RDa")


##############Step 1: Extract 50 marker genes for every cell lineage of the reference cells############
########################################################################################################
snd <- CreateSeuratObject(counts = Reference_counts_logCPM, project = "mouse_reference")

###add information on cell lineages to the reference samples and set the cell lineages as actual indent
snd <- AddMetaData(object = snd,
                   metadata = Reference_connectionDF$cell_lineage,
                   col.name = "cell_lineage")
###manually set idents
snd <- SetIdent(snd, value = "cell_lineage")

###Calculate all marker genes, distinguishing the different cell lineages
snd.markers <- FindAllMarkers(snd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#extract top 50 marker genes for every cell lineage, based on avg fold change
top50 <- snd.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top50 <- assign.GeneNames(EnsemblIDs = top50$gene,target.DF = top50)









##############Step 2: Correlate the log-normalized gene expression of the single cells with the ###########
#reference samples. Assign celltype label based on highest correlation with the reference samples.#########
#Important: For the celltype assignment, the specific subpopulation information of the T Cells is considered
###########################################################################################################

###First, realebel the T Cells of the reference samples to their subtype (naive CD4, naive CD8,...)
tcell_sub <- subset(Reference_connectionDF,cell_lineage == "T Cell Lineage")
Reference_connectionDF <- Reference_connectionDF[!((Reference_connectionDF$cell_lineage == "T Cell Lineage")),]

tcell_sub$cell_lineage <- tcell_sub$cell_type_short
total.cd4cd8 <- c("CD4T","CD8T")
tcell_sub <- tcell_sub[!(tcell_sub$cell_lineage %in% total.cd4cd8),]
Reference_connectionDF <- rbind(Reference_connectionDF,tcell_sub)

###Log-Normalize QC filtered scSamples

counts.corrCPM <- correctedCPM(dat.counts.blood.QCpassed)
log_normCPM.counts <- log_xp1(counts.corrCPM)

###keep only genes which are among the top50 marker genes per lineage
sc.Data <- log_normCPM.counts[which(rownames(log_normCPM.counts)%in%top50$gene),]

###subset the reference cells accordingly, to contain the same set of genes
ref <- Reference_counts_logCPM[which(rownames(Reference_counts_logCPM)%in%rownames(sc.Data)),which(colnames(Reference_counts_logCPM)%in%rownames(Reference_connectionDF))]
summary(rownames(sc.Data)==rownames(ref))

###########Perform correlation for every sample of the single cells with the reference samples###########
###Calculate correlation Matrix
cor.mat <- data.frame(matrix(nrow = ncol(sc.Data),ncol = ncol(ref)))
rownames(cor.mat) <- colnames(sc.Data)
colnames(cor.mat) <- colnames(ref)

for(i in seq(1,ncol(sc.Data))){
  curr <- sc.Data[,i]
  for(k in seq(1,ncol(ref))){
    y <- ref[,k]
    cor.coeff <- cor(curr,y,method = "pearson")
    cor.mat[i,k] <- cor.coeff
  }
}

##########Assign every single cell sample to the reference sample with the highes correlation value
#

correlation.mat <- cor.mat
max.cor.df <- as.data.frame(matrix(nrow=nrow(correlation.mat)))
max.cor.df[,1] <- NULL
rownames(max.cor.df) <- rownames(correlation.mat)

N <- 1
for (i in seq(1,nrow(correlation.mat))){
  ordered_cor <- order(correlation.mat[i,], decreasing = T)
  
  #assign top correlation cell lineage
  top.idx <- ordered_cor[N]#[1:N]
  max.cor.df[i,1] <- as.numeric(correlation.mat[i,top.idx])
  max.cor.df[i,2] <- colnames(correlation.mat)[top.idx]
  max.cor.df[i,3] <- Reference_connectionDF$cell_lineage[match(max.cor.df[i,2],rownames(Reference_connectionDF))]
  
  #assign all cell lineages to all samples according to reference sample
  top.cell.lineages <- (correlation.mat)[i,ordered_cor]
  top.cell.lineages[2,] <- Reference_connectionDF$cell_lineage[match(colnames(top.cell.lineages),rownames(Reference_connectionDF))]
  top.cell.lineages <- t(top.cell.lineages)
  
  #check which of the reference cells is the first one with a different cell_lineage than the one with the highest correlation
  top_lineage <- top.cell.lineages[1,2]
  z <- 1
  while (top.cell.lineages[z,2] == top_lineage){
    z <- z+1
  }
  diff.lineage.ct <- rownames(top.cell.lineages)[z]
  diff.lineage.corr <- top.cell.lineages[z,1]
  diff.lineage.lineage <- top.cell.lineages[z,2]
  
  max.cor.df[i,(4:6)] <- c(diff.lineage.lineage,diff.lineage.ct,diff.lineage.corr)
  colnames(max.cor.df) <- c("max_corr","celltype","cell_lineage","first_diff.cell_lineage","first_diff.cell_ct","first_diff.cell.corr")
}
max.cor.df$corr_difference <- max.cor.df$max_corr - as.numeric(max.cor.df$first_diff.cell.corr)
max.cor.df$corr_difference_ratio <- max.cor.df$max_corr / as.numeric(max.cor.df$first_diff.cell.corr)

table(max.cor.df$cell_lineage)
assigned.celltypes_byCorrelation <- max.cor.df

save(assigned.celltypes_byCorrelation,file = "RDataFiles/CorrrelationAssignmentResult.RDa")




