setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/")
library(MAST)
library(scater)
library(SummarizedExperiment)

load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")


##Methylation info.
agePredictionScOrg <- read.delim("./Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p1.extended.txt",as.is=T)
agePredictionScOrg2 <- read.delim("./Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p2.extended.txt",as.is=T) 

#Drop non-passing QC cells.
qcInfoP1 = read.delim("./Trapp_Re-engineerd/passQc_Cells.txt",header = F)
qcInfoP2 = read.delim("./Trapp_Re-engineerd/passQc_Cells_p2.txt",header = F)

agePredictionScOrg = agePredictionScOrg[grep(paste(qcInfoP1$V1,collapse = "|"),rownames(agePredictionScOrg)),]
agePredictionScOrg2 = agePredictionScOrg2[grep(paste(qcInfoP2$V1,collapse = "|"),rownames(agePredictionScOrg2)),]

##Make into short IDs to do this analysis easier.
for(rN in 1:nrow(agePredictionScOrg)){
  tmp = strsplit(rownames(agePredictionScOrg)[rN],"_")
  tmp = paste0(tmp[[1]][4],"_",tmp[[1]][5])
  tmp = gsub(pattern = "lood",tmp,replacement = "")
  rownames(agePredictionScOrg)[rN] = tmp
}

for(rN in 1:nrow(agePredictionScOrg2)){
  tmp = strsplit(rownames(agePredictionScOrg2)[rN],"_")
  tmp = paste0(tmp[[1]][2],"_",tmp[[1]][3])
  tmp = gsub(pattern = "lood",tmp,replacement = "")
  #print(tmp)
  tmp = gsub(pattern = "B3",tmp,replacement = "B03")
  rownames(agePredictionScOrg2)[rN] = tmp
}

agePredictionScOrg = rbind(agePredictionScOrg,agePredictionScOrg2)

agePredictionScOrg["predictedAgeCor"] = agePredictionScOrg$predictedAge - (agePredictionScOrg$medianRandomAge - agePredictionScOrg$actualAge)

colDataCells = cbind(colDataCells,agePredictionScOrg[match(colDataCells$babraham_short, rownames(agePredictionScOrg)),c(4:10)])

aggregateTcells = F
estimateHbLevels = T
removeHbCells = F
## replaceHighFdrAgeValues = F ##To be implemented

se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(dat.counts.blood.QCpassed)), colData=colDataCells)

rm(dat.counts.blood.QCpassed)
se0 <- addPerCellQC(se0)
se0 <- addPerFeatureQC(se0)
se0 <- scater::logNormCounts(se0)

plotInfo = F

if(estimateHbLevels){
  hbGenes = c("ENSMUSG00000073940", "ENSMUSG00000052305", "ENSMUSG00000069919", "ENSMUSG00000069917")
  
  ##Add code to flag HB cells / and remove them.
  sce0.hb = as(se0[which(rownames(se0) %in% hbGenes),], "SingleCellExperiment")
  
  if(plotInfo){
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
  
  if(plotInfo){
    View(annotHbGenes)
    hist(annotHbGenes[,6])
  }
  
  colData(se0)["hbProportion"] = annotHbGenes[,6]
  colData(se0)["hbProportionFine"] = annotHbGenes[,5]
  
  #se0$hbProportion
  #se0$hbProportionFine
}
##

if(removeHbCells){
  se0 = se0[,which(se0$hbProportion==0)]
}

##Filter expressed in 25% of the cells (cell type specific), protein coding only

gtf_df <- as.data.frame(rtracklayer::import('./Mus_musculus.GRCm38.96.gtf.gz'))
gtf_df  = gtf_df[which(gtf_df$gene_biotype=="protein_coding"),]

##Select B-cells (to single cell assay)
sceB = as(se0[,which(se0$annotation=="BCell" & !is.na(se0$predictedAge))], "SingleCellExperiment")
rowData(sceB) = NULL

sceB <- addPerFeatureQC(sceB)
sceB = sceB[which(rowData(sceB)$detected>=25),]
sceB = sceB[which(rownames(sceB) %in% gtf_df$gene_id),]

##Check if in at least 2 ages there are at least 5 cells.
testableGenes = list()
for(ag in unique(sceB$age_weeks.group)){
  testableGenes[[ag]] = names(which(rowSums(logcounts(sceB)[,which(sceB$age_weeks.group==ag)]>0)>4))
}
testableGenes = names(which(table(unlist(testableGenes))>1))
sceB = sceB[which(rownames(sceB) %in% testableGenes),]

scaB = SceToSingleCellAssay(sceB, class = "SingleCellAssay")

colData(scaB)$detected= scale(colData(scaB)$detected) # detected == n_gene (CDR)
colData(scaB)$mouse= as.factor(colData(scaB)$mouse)
colData(scaB)$annotation= as.factor(colData(scaB)$annotation)
scaB = scaB[rowSums(assay(scaB)) != 0, ]

##Select T-cells (to single cell assay)
sceT = as(se0[,intersect(grep(se0$annotation,pattern = "T"),which(!is.na(se0$predictedAge)))], "SingleCellExperiment")
rowData(sceT) = NULL

sceT <- addPerFeatureQC(sceT)
sceT = sceT[which(rowData(sceT)$detected>=25),]
sceT = sceT[which(rownames(sceT) %in% gtf_df$gene_id),]

##Check if in at least 2 ages there are at least 5 cells.
testableGenes = list()
for(ag in unique(sceT$age_weeks.group)){
  testableGenes[[ag]] = names(which(rowSums(logcounts(sceT)[,which(sceT$age_weeks.group==ag)]>0)>4))
}
testableGenes = names(which(table(unlist(testableGenes))>1))
sceT = sceT[which(rownames(sceT) %in% testableGenes),]

scaT = SceToSingleCellAssay(sceT, class = "SingleCellAssay")

colData(scaT)$detected = scale(colData(scaT)$detected) # detected == n_gene (CDR)
colData(scaT)$mouse= as.factor(colData(scaT)$mouse)
colData(scaT)$annotation= as.factor(colData(scaT)$annotation)
scaT = scaT[rowSums(assay(scaT)) != 0, ]

rm(sceT,sceB)


##Testing flags
flag_n_gene=TRUE
flag_correct_Hb=T
flag_correct_chronAge = F
flag_correct_HbF=F
flag_subCellType=T
flag_subCellType_randomEff=F
flag_mouse=T
flag_mouse_randomEff=F

print(paste0('Remove HB cells: ', removeHbCells))
print(paste0('flag_correct_Hb: ', flag_correct_Hb))
print(paste0('flag_correct_chronAge: ', flag_correct_chronAge))
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
if (flag_correct_chronAge==TRUE){
  covariate = paste0(covariate, " + age_weeks.group")
}

if (flag_mouse==TRUE & flag_mouse_randomEff!=TRUE){
  covariate = paste0(covariate, " + mouse")
}
if (flag_mouse_randomEff==TRUE){
  covariate = paste0(covariate, " + (1|mouse)")
}

print(covariate)

##Test B cells (Twice)
if(!flag_subCellType_randomEff & !flag_mouse_randomEff){
  
  zlmCond_B <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), sca=scaB, parallel = TRUE)
  
  if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
    covariate = paste0(covariate, " + annotation")
  }
  zlmCond_T <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), sca=scaT, parallel = TRUE)
  
} else{
  zlmCond_B <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), method='glmer', sca=scaB, parallel = TRUE, ebayes=FALSE)
  if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
    covariate = paste0(covariate, " + annotation")
  }
  if (flag_subCellType_randomEff==TRUE){
    covariate = paste0(covariate, " + (1|annotation)")
  }
  zlmCond_T <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), method='glmer', sca=scaT, parallel = TRUE, ebayes=FALSE)
}

term = "predictedAge"

summaryCond_B <- summary(zlmCond_B, doLRT=term)
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
write.csv(de_res,paste0('./Bcells_epi_Age_Bonder.csv'), row.names=FALSE)

summaryCond_T <- summary(zlmCond_T, doLRT="predictedAge")
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
write.csv(de_res,paste0('./Tcells_epi_Aging_Bonder_subcellTypeCorrected.csv'), row.names=FALSE)

##Per age.
if(!flag_correct_chronAge){
  for(ageG in unique(se0$age_weeks.group) ){
    
    ##Select B-cells (to single cell assay)
    sceB = as(se0[,which(se0$annotation=="BCell" & !is.na(se0$predictedAge) & se0$age_weeks.group==ageG)], "SingleCellExperiment")
    rowData(sceB) = NULL
    sceB <- addPerFeatureQC(sceB)
    
    ##Check if the gene is present in at least 25% of the cells.
    sceB = sceB[which(rowData(sceB)$detected>=25),]
    
    scaB = SceToSingleCellAssay(sceB, class = "SingleCellAssay")
    
    colData(scaB)$detected= scale(colData(scaB)$detected) # detected == n_gene (CDR)
    colData(scaB)$mouse= as.factor(colData(scaB)$mouse)
    colData(scaB)$annotation= as.factor(colData(scaB)$annotation)
    scaB = scaB[rowSums(assay(scaB)) != 0, ]
    
    
    ##Select T-cells (to single cell assay)
    sceT = as(se0[,intersect(grep(se0$annotation,pattern = "T"),which(!is.na(se0$predictedAge) & se0$age_weeks.group==ageG))], "SingleCellExperiment")
    rowData(sceT) = NULL
    sceT <- addPerFeatureQC(sceT)
    
    ##Check if the gene is present in at least 25% of the cells.
    sceT = sceT[which(rowData(sceT)$detected>=25),]
    
    scaT = SceToSingleCellAssay(sceT, class = "SingleCellAssay")
    
    colData(scaT)$detected = scale(colData(scaT)$detected) # detected == n_gene (CDR)
    colData(scaT)$mouse= as.factor(colData(scaT)$mouse)
    colData(scaT)$annotation= as.factor(colData(scaT)$annotation)
    scaT = scaT[rowSums(assay(scaT)) != 0, ]
    
    rm(sceT,sceB)
    
    
    ##Testing.
    flag_n_gene=TRUE
    flag_correct_Hb=T
    
    flag_correct_HbF=F
    flag_subCellType=T
    flag_subCellType_randomEff=F
    flag_mouse=F
    flag_mouse_randomEff=F
    
    print(paste0('Remove HB cells: ', removeHbCells))
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
    
    if (flag_mouse==TRUE & flag_mouse_randomEff!=TRUE){
      covariate = paste0(covariate, " + mouse")
    }
    if (flag_mouse_randomEff==TRUE){
      covariate = paste0(covariate, " + (1|mouse)")
    }
    
    print(covariate)
    
    ##Test B cells (Twice)
    if(!flag_subCellType_randomEff & !flag_mouse_randomEff){
      
      zlmCond_B <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), sca=scaB, parallel = TRUE)
      
      if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
        covariate = paste0(covariate, " + annotation")
      }
      if (flag_subCellType_randomEff==TRUE){
        covariate = paste0(covariate, " + (1|annotation)")
      }
      zlmCond_T <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), sca=scaT, parallel = TRUE)
      
    } else{
      zlmCond_B <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), method='glmer', sca=scaB, parallel = TRUE, ebayes=FALSE)
      if (flag_subCellType==TRUE & flag_subCellType_randomEff!=TRUE){
        covariate = paste0(covariate, " + annotation")
      }
      if (flag_subCellType_randomEff==TRUE){
        covariate = paste0(covariate, " + (1|annotation)")
      }
      zlmCond_T <- zlm(formula = as.formula(paste0("~predictedAge", covariate)), method='glmer', sca=scaT, parallel = TRUE, ebayes=FALSE)
    }
    
    summaryCond_B <- summary(zlmCond_B, doLRT="predictedAge")
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
    write.csv(de_res,paste0('./Bcells_epi_Age_Bonder_',ageG,'.csv'), row.names=FALSE)
    
    summaryCond_T <- summary(zlmCond_T, doLRT="predictedAge")
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
    write.csv(de_res,paste0('./Tcells_epi_Aging_Bonder_subcellTypeCorrected_',ageG,'.csv'), row.names=FALSE)
    
  }
}
