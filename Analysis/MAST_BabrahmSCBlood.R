setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/")
library(MAST)
library(scater)
library(SummarizedExperiment)
library(data.table)

load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")

##Main Settings
aggregateTcells = T
pseudoBulking = T
estimateHbLevels = T
plotHbInfo = F

##Testing settings
flag_n_gene=TRUE
flag_correct_Hb=T
flag_correct_HbF=F
flag_subCellType=T
flag_subCellType_randomEff=F
flag_mouse=T
flag_mouse_randomEff=F

##First do single cell QC + HB estimate.

se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(dat.counts.blood.QCpassed)), colData=colDataCells)
se0 <- addPerCellQC(se0)
se0 <- addPerFeatureQC(se0)
se0 <- scater::logNormCounts(se0)

if(estimateHbLevels){
  hbGenes = c("ENSMUSG00000073940", "ENSMUSG00000052305", "ENSMUSG00000069919", "ENSMUSG00000069917")
  
  ##Add code to flag HB cells / and remove them.
  sce0.hb = as(se0[which(rownames(se0) %in% hbGenes),], "SingleCellExperiment")
  
  if(plotHbInfo){
    vioplot::vioplot(logcounts(sce0.hb)[1,])
    vioplot::vioplot(logcounts(sce0.hb)[2,])
    vioplot::vioplot(logcounts(sce0.hb)[3,])
    vioplot::vioplot(logcounts(sce0.hb)[4,])
    
    hist(colSums(logcounts(sce0.hb)))
    hist(colSums(logcounts(sce0.hb[c(1,4),])))
    hist(colSums(logcounts(sce0.hb[c(2:3),])))
    
    hist((logcounts(sce0.hb)[1,]))
    hist((logcounts(sce0.hb)[2,]))
    hist((logcounts(sce0.hb)[3,]))
    hist((logcounts(sce0.hb)[4,]))
  }
  
  g1Annot = rep(0,ncol(sce0.hb))
  g1Annot[which(logcounts(sce0.hb)[1,]>0)]=1
  g1Annot[which(logcounts(sce0.hb)[1,]>10)]=2
  
  g2Annot = rep(0,ncol(sce0.hb))
  g2Annot[which(logcounts(sce0.hb)[2,]>0)]=1
  g2Annot[which(logcounts(sce0.hb)[2,]>12)]=2
  
  g3Annot = rep(0,ncol(sce0.hb))
  g3Annot[which(logcounts(sce0.hb)[3,]>0)]=1
  g3Annot[which(logcounts(sce0.hb)[3,]>12)]=2
  
  g4Annot = rep(0,ncol(sce0.hb))
  g4Annot[which(logcounts(sce0.hb)[4,]>0)]=1
  g4Annot[which(logcounts(sce0.hb)[4,]>10)]=2
  
  annotHbGenes = cbind(g1Annot,g2Annot,g3Annot,g4Annot)
  colnames(annotHbGenes) = rownames(sce0.hb)
  annotHbGenes = cbind(annotHbGenes,rowMeans(annotHbGenes),ceiling(rowMeans(annotHbGenes)))
  rm(g1Annot,g2Annot,g3Annot,g4Annot)
  rownames(annotHbGenes) = colnames(se0)
  
  if(plotHbInfo){
    View(annotHbGenes)
    hist(annotHbGenes[,6])
  }
  
  colData(se0)["hbProportion"] = annotHbGenes[,6]
  colData(se0)["hbProportionFine"] = annotHbGenes[,5]
  
  colDataCells["hbProportion"] = annotHbGenes[,6]
  colDataCells["hbProportionFine"] = annotHbGenes[,5]
  
  se0$hbProportion
  se0$hbProportionFine
  
}

if(pseudoBulking){
  ##Add code to pseudo bulk.
  
  if(aggregateTcells){
    ##Calcualte T subtype fractions.
    tcellFractions = list()
    for(m in unique(colDataCells$mouse)){
      infoTable = table(colDataCells$annotation[intersect(which(colDataCells$mouse==m),grep(colDataCells$annotation,pattern = "T"))])
      tcellFractions[[m]] = infoTable / sum(infoTable)
    }
    ##Adding empty columns.
    columnsToAdd = unique(names(unlist(tcellFractions)))
    colDataCells[columnsToAdd] = NA
    for(m in unique(colDataCells$mouse)){
      for(cId in columnsToAdd){
        if(length(which(names(tcellFractions[[m]])==cId))==0){
          colDataCells[which(colDataCells$mouse==m),which(colnames(colDataCells) == cId)] = 0
        } else {
          colDataCells[which(colDataCells$mouse==m),which(colnames(colDataCells) == cId)] = tcellFractions[[m]][which(names(tcellFractions[[m]])==cId)]
        }
      }
    }
    colDataCells$annotation[grep(colDataCells$annotation,pattern = "T")] = "TCell"
  }
  
  ##Add average HB proportion.
  if(estimateHbLevels){
    for(m in unique(colDataCells$mouse)){
      colDataCells$hbProportion[which(colDataCells$mouse==m)] = mean(colDataCells$hbProportion[which(colDataCells$mouse==m)])
      colDataCells$hbProportionFine[which(colDataCells$mouse==m)] = mean(colDataCells$hbProportionFine[which(colDataCells$mouse==m)])
    }
  }
  
  colDataCells["pbCt"] = paste0(colDataCells$mouse,"_",colDataCells$annotation)
  pbcolDataCells = unique(colDataCells[,c(7,8,9,10,14,16:ncol(colDataCells))])
  
  pb.dat.counts.blood.QCpassed = as.data.frame(t(dat.counts.blood.QCpassed[which(rowSums(dat.counts.blood.QCpassed)>0),]))
  pb.dat.counts.blood.QCpassed["pbCt"] = colDataCells["pbCt"]
  
  pb.dat.counts.blood.QCpassed = data.table(pb.dat.counts.blood.QCpassed)
  pb.dat.counts.blood.QCpassed = pb.dat.counts.blood.QCpassed[, lapply(.SD, sum, na.rm=TRUE), by=pbCt]
  
  pb.dat.counts.blood.QCpassed = as.data.frame(pb.dat.counts.blood.QCpassed)
  rownames(pb.dat.counts.blood.QCpassed) = pb.dat.counts.blood.QCpassed[,1]
  pb.dat.counts.blood.QCpassed = pb.dat.counts.blood.QCpassed[,2:ncol(pb.dat.counts.blood.QCpassed)]
  
  se0 <- SummarizedExperiment(assays=SimpleList(counts=t(as.matrix(pb.dat.counts.blood.QCpassed))), colData=pbcolDataCells)
  rm(pb.dat.counts.blood.QCpassed,pbcolDataCells)
  
  ##Need to do QC again!
  se0 <- addPerCellQC(se0)
  se0 <- addPerFeatureQC(se0)
  se0 <- scater::logNormCounts(se0)
}

rm(dat.counts.blood.QCpassed)

####Filter expressed in 25% of the cells (cell type specific), protein coding only

gtf_df <- as.data.frame(rtracklayer::import('./Mus_musculus.GRCm38.96.gtf.gz'))
gtf_df  = gtf_df[which(gtf_df$gene_biotype=="protein_coding"),]

##Select B-cells (to single cell assay)
sceB = as(se0[,which(se0$annotation=="BCell")], "SingleCellExperiment")
rowData(sceB) = NULL

sceB <- addPerFeatureQC(sceB)
sceB = sceB[which(rowData(sceB)$detected>25),]
sceB = sceB[which(rownames(sceB) %in% gtf_df$gene_id),]

##Do two filter steps in parallel.

par(mfrow=c(1,2))
dec.sca <- scran::modelGeneVarByPoisson(sceB)
plot(dec.sca$mean, dec.sca$total, pch=16,cex=.6, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec.sca)$trend(x), col="dodgerblue", add=TRUE)

chosen <- scran::getTopHVGs(dec.sca, var.threshold = 1)

plot(dec.sca[chosen,]$mean, dec.sca[chosen,]$total, pch=16, cex=.6, xlab="Mean of log-expression",
     ylab="Variance of log-expression",ylim=c(0, max(dec.sca[chosen,]$total)))
curve(metadata(dec.sca)$trend(x), col="dodgerblue", add=TRUE)

par(mfrow=c(1,1))

##Check if in at least 2 ages there are at least 5 cells.
testableGenes = list()
for(ag in unique(sceB$age_weeks.group)){
  if(!pseudoBulking){
    testableGenes[[ag]] = names(which(rowSums(logcounts(sceB)[,which(sceB$age_weeks.group==ag)]>0)>4))
  } else {
    testableGenes[[ag]] = names(which(rowSums(logcounts(sceB)[,which(sceB$age_weeks.group==ag)]>0)>1))
  }
}
testableGenes = names(which(table(unlist(testableGenes))>1))

##Intersect two lists
testableGenes = intersect(chosen,testableGenes)
sceB = sceB[which(rownames(sceB) %in% testableGenes),]

scaB = SceToSingleCellAssay(sceB, class = "SingleCellAssay")

colData(scaB)$detected= scale(colData(scaB)$detected) # detected == n_gene (CDR)
colData(scaB)$mouse= as.factor(colData(scaB)$mouse)
colData(scaB)$annotation= as.factor(colData(scaB)$annotation)
scaB = scaB[rowSums(assay(scaB)) != 0, ]


##Select T-cells (to single cell assay)
sceT = as(se0[,grep(se0$annotation,pattern = "T")], "SingleCellExperiment")
rowData(sceT) = NULL

sceT <- addPerFeatureQC(sceT)
sceT = sceT[which(rowData(sceT)$detected>25),]
sceT = sceT[which(rownames(sceT) %in% gtf_df$gene_id),]

par(mfrow=c(1,2))
dec.sca <- scran::modelGeneVarByPoisson(sceT)
plot(dec.sca$mean, dec.sca$total, pch=16,cex=.6, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec.sca)$trend(x), col="dodgerblue", add=TRUE)

chosen <- scran::getTopHVGs(dec.sca, var.threshold = 1)

plot(dec.sca[chosen,]$mean, dec.sca[chosen,]$total, pch=16, cex=.6, xlab="Mean of log-expression",
     ylab="Variance of log-expression",ylim=c(0, max(dec.sca[chosen,]$total)))
curve(metadata(dec.sca)$trend(x), col="dodgerblue", add=TRUE)
par(mfrow=c(1,1))

##Check if in at least 2 ages there are at least 5 cells.
testableGenes = list()
for(ag in unique(sceT$age_weeks.group)){
  if(!pseudoBulking){
    testableGenes[[ag]] = names(which(rowSums(logcounts(sceT)[,which(sceT$age_weeks.group==ag)]>0)>4))
  } else {
    testableGenes[[ag]] = names(which(rowSums(logcounts(sceT)[,which(sceT$age_weeks.group==ag)]>0)>1))  
  }
}
testableGenes = names(which(table(unlist(testableGenes))>1))
##Intersect two lists
testableGenes = intersect(chosen,testableGenes)
sceT = sceT[which(rownames(sceT) %in% testableGenes),]

scaT = SceToSingleCellAssay(sceT, class = "SingleCellAssay")

colData(scaT)$detected = scale(colData(scaT)$detected) # detected == n_gene (CDR)
colData(scaT)$mouse= as.factor(colData(scaT)$mouse)
colData(scaT)$annotation= as.factor(colData(scaT)$annotation)
scaT = scaT[rowSums(assay(scaT)) != 0, ]

rm(sceT,sceB)


print(paste0('flag_correct_Hb: ', flag_correct_Hb))
print(paste0('flag_n_gene: ', flag_n_gene))
print(paste0('flag_subCellType: ', flag_subCellType))
print(paste0('flag_subCellType_randomEff: ', flag_subCellType_randomEff))
print(paste0('flag_mouse: ', flag_mouse))
print(paste0('flag_mouse_randomEff: ', flag_mouse_randomEff))

# DGE testing 
covariate = ''
if (flag_n_gene==TRUE){
  covariate = paste0(covariate, " + detected")
}
if (flag_correct_Hb==TRUE){
  covariate = paste0(covariate, " + hbProportion")
}
if (flag_mouse==TRUE & flag_mouse_randomEff!=TRUE){
  if(!pseudoBulking){
    covariate = paste0(covariate, " + mouse")  
  }
}
if (flag_mouse_randomEff==TRUE){
  covariate = paste0(covariate, " + (1|mouse)")
}

print(covariate)

##Test B cells (Twice)
if(!flag_subCellType_randomEff & !flag_mouse_randomEff){
  zlmCond_B <- zlm(formula = as.formula(paste0("~age_weeks.group", covariate)), sca=scaB, parallel = TRUE)
  if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
    if(!pseudoBulking){
      covariate = paste0(covariate, " + annotation")
    } else {
      covariate = paste0(covariate, " + EffCD4T + MemCD8T + MemCD4T + NveCD4T + NveCd8T + RegT")
    }
  }
  if (flag_subCellType_randomEff==TRUE){
    covariate = paste0(covariate, " + (1|annotation)")
  }
  zlmCond_T <- zlm(formula = as.formula(paste0("~age_weeks.group", covariate)), sca=scaT, parallel = TRUE)
  
} else {
  zlmCond_B <- zlm(formula = as.formula(paste0("~age_weeks.group", covariate)), method='glmer', sca=scaB, parallel = TRUE, ebayes=FALSE)
  if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
    covariate = paste0(covariate, " + annotation")
  }
  if (flag_subCellType_randomEff==TRUE){
    covariate = paste0(covariate, " + (1|annotation)")
  }
  zlmCond_T <- zlm(formula = as.formula(paste0("~age_weeks.group", covariate)), method='glmer', sca=scaT, parallel = TRUE, ebayes=FALSE)
}
term = "age_weeks.group"
summaryCond_B <- summary(zlmCond_B, doLRT="age_weeks.group")
# Summarize results
summaryDt <- summaryCond_B$datatable
dt1 = summaryDt[contrast==term & component=="H", .(primerid, `Pr(>Chisq)`)]
dt2 = summaryDt[contrast==term & component=="logFC", .(primerid, coef, z)]
dt3 = summaryDt[contrast==term & component=="D", .(primerid, coef, z,`Pr(>Chisq)`)]
dt4 = summaryDt[contrast==term & component=="C", .(primerid, coef, z,`Pr(>Chisq)`)]
de_res = merge(dt2, dt1, by="primerid")
de_res = merge(de_res, dt3, by="primerid")
de_res = merge(de_res, dt4, by="primerid")

colnames(de_res) <- c("gene", "age.logFC", 'age.logFC_z', "age.H_p", "age.D.coef", 'age.D_z', "age.D_p", "age.C.coef", 'age.C_z', "age.C_p")
##Add expression level and expressed
tempStats = perFeatureQCMetrics(scaB)
de_res = cbind(de_res,as.matrix(tempStats[match(de_res$gene,rownames(tempStats)),]))

##Remove effects that are nominal significant but in opposite direction.
toDrop = c(intersect(which(de_res$age.C_z<0 & de_res$age.D_z>0),intersect(which(de_res$age.C_p<0.05),which(de_res$age.D_p<0.05))),
           intersect(which(de_res$age.C_z>0 & de_res$age.D_z<0),intersect(which(de_res$age.C_p<0.05),which(de_res$age.D_p<0.05))))
de_res$age.H_p[toDrop] = NA
de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")

if(pseudoBulking){
  write.csv(de_res[order(de_res$age.H_fdr),],paste0('./Bcells_pb_Aging_Aug2022.csv'), row.names=FALSE)
} else {
  write.csv(de_res[order(de_res$age.H_fdr),],paste0('./Bcells_sc_Aging_Aug2022.csv'), row.names=FALSE)
}

summaryCond_T <- summary(zlmCond_T, doLRT="age_weeks.group")
# Summarize results
summaryDt <- summaryCond_T$datatable
dt1 = summaryDt[contrast==term & component=="H", .(primerid, `Pr(>Chisq)`)]
dt2 = summaryDt[contrast==term & component=="logFC", .(primerid, coef, z)]
dt3 = summaryDt[contrast==term & component=="D", .(primerid, coef, z,`Pr(>Chisq)`)]
dt4 = summaryDt[contrast==term & component=="C", .(primerid, coef, z,`Pr(>Chisq)`)]
de_res = merge(dt2, dt1, by="primerid")
de_res = merge(de_res, dt3, by="primerid")
de_res = merge(de_res, dt4, by="primerid")

colnames(de_res) <- c("gene", "age.logFC", 'age.logFC_z', "age.H_p", "age.D.coef", 'age.D_z', "age.D_p", "age.C.coef", 'age.C_z', "age.C_p")
##Add expression level and expressed
tempStats = perFeatureQCMetrics(scaT)
de_res = cbind(de_res,as.matrix(tempStats[match(de_res$gene,rownames(tempStats)),]))

##Remove effects that are nominal significant but in opposite direction.
toDrop = c(intersect(which(de_res$age.C_z<0 & de_res$age.D_z>0),intersect(which(de_res$age.C_p<0.05),which(de_res$age.D_p<0.05))),
           intersect(which(de_res$age.C_z>0 & de_res$age.D_z<0),intersect(which(de_res$age.C_p<0.05),which(de_res$age.D_p<0.05))))
de_res$age.H_p[toDrop] = NA
de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
if(pseudoBulking){
  write.csv(de_res[order(de_res$age.H_fdr),],paste0('./Tcells_pb_Aging_Bonder_subcellTypeCorrected_Aug2022.csv'), row.names=FALSE)
} else {
  write.csv(de_res[order(de_res$age.H_fdr),],paste0('./Tcells_sc_Aging_Bonder_subcellTypeCorrected_Aug2022.csv'), row.names=FALSE)
}


