setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/")
library(MAST)
library(scater)
library(SummarizedExperiment)
library(data.table)

load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")
rownames(colDataCells) = colDataCells$X

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
rm(agePredictionScOrg2)

#Drop predictions with less then 40 sites.
agePredictionScOrg = agePredictionScOrg[which(agePredictionScOrg$sitesUsed>39),]

agePredictionScOrg["predictedAgeCor"] = agePredictionScOrg$predictedAge - (agePredictionScOrg$medianRandomAge - agePredictionScOrg$actualAge)

colDataCells = cbind(colDataCells,agePredictionScOrg[match(colDataCells$babraham_short, rownames(agePredictionScOrg)),c(4:10)])

##Main Settings
aggregateTcells = T
pseudoBulking = F
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
#gtf_df <-as.data.frame(rtracklayer::import("Mus_musculus.GRCm38.96.gtf.gz"))
gtf<-read.delim("./Mus_musculus.GRCm38.96.gtf.gz",as.is=T,skip=5,header=F)
gtf <- gtf[which(gtf$V3=="gene"),]
gtf = cbind(gtf,matrix(unlist(strsplit(gtf$V9,split="; ")),ncol = 5,byrow = T))
gtf = gtf[,-9]
##Fix new columns.
gtf[,9] = gsub(gtf[,9],pattern = "gene_id ",replacement = "")
gtf[,10] = gsub(gtf[,10],pattern = "gene_version ",replacement = "")
gtf[,11] = gsub(gtf[,11],pattern = "gene_name ",replacement = "")
gtf[,12] = gsub(gtf[,12],pattern = "gene_source ",replacement = "")
gtf[,13] = gsub(gtf[,13],pattern = "gene_biotype ",replacement = "")
gtf[,13] = gsub(gtf[,13],pattern = ";",replacement = "")
colnames(gtf) = c("chr","source","feature_type","start","end",".","strand","","gene_id","version","name","source2","gene_biotype")
gtf_df <- as.data.frame(gtf)
rm(gtf)
gtf_df  = gtf_df[which(gtf_df$gene_biotype=="protein_coding"),]

##All cells
sce = as(se0[,which(!is.na(se0$annotation) & !is.na(se0$predictedAge))], "SingleCellExperiment")
rowData(sce) = NULL

sce <- addPerFeatureQC(sce)
sce = sce[which(rowData(sce)$detected>25),]
sce = sce[which(rownames(sce) %in% gtf_df$gene_id),]

scaB = SceToSingleCellAssay(sce, class = "SingleCellAssay")

colData(scaB)$detected= colData(scaB)$detected # detected == n_gene (CDR)
colData(scaB)$mouse= as.factor(colData(scaB)$mouse)
colData(scaB)$annotation= as.factor(colData(scaB)$annotation)
sca = scaB[rowSums(assay(scaB)) != 0, ]


##Select B-cells (to single cell assay)
sceB = as(se0[,which(se0$annotation=="BCell")], "SingleCellExperiment")
rowData(sceB) = NULL

sceB <- addPerFeatureQC(sceB)
sceB = sceB[which(rowData(sceB)$detected>25),]
sceB = sceB[which(rownames(sceB) %in% gtf_df$gene_id),]

scaB = SceToSingleCellAssay(sceB, class = "SingleCellAssay")

colData(scaB)$detected= colData(scaB)$detected # detected == n_gene (CDR)
colData(scaB)$mouse= as.factor(colData(scaB)$mouse)
colData(scaB)$annotation= as.factor(colData(scaB)$annotation)
scaB = scaB[rowSums(assay(scaB)) != 0, ]


##Select T-cells (to single cell assay)
sceT = as(se0[,grep(se0$annotation,pattern = "T")], "SingleCellExperiment")
rowData(sceT) = NULL

sceT <- addPerFeatureQC(sceT)
sceT = sceT[which(rowData(sceT)$detected>25),]
sceT = sceT[which(rownames(sceT) %in% gtf_df$gene_id),]

scaT = SceToSingleCellAssay(sceT, class = "SingleCellAssay")

colData(scaT)$detected = (colData(scaT)$detected) # detected == n_gene (CDR)
colData(scaT)$mouse= as.factor(colData(scaT)$mouse)
colData(scaT)$annotation= as.factor(colData(scaT)$annotation)
scaT = scaT[rowSums(assay(scaT)) != 0, ]

##sceT: cd4 & cd8.
##cd8
sceT8 = as(se0[,intersect(grep(se0$annotation,pattern = "T"),grep(se0$annotation,pattern = "8"))], "SingleCellExperiment")
rowData(sceT8) = NULL

sceT8 <- addPerFeatureQC(sceT8)
sceT8 = sceT8[which(rowData(sceT8)$detected>25),]
sceT8 = sceT8[which(rownames(sceT8) %in% gtf_df$gene_id),]

scaT8 = SceToSingleCellAssay(sceT8, class = "SingleCellAssay")

colData(scaT8)$detected = (colData(scaT8)$detected) # detected == n_gene (CDR)
colData(scaT8)$mouse= as.factor(colData(scaT8)$mouse)
colData(scaT8)$annotation= as.factor(colData(scaT8)$annotation)
scaT8 = scaT8[rowSums(assay(scaT8)) != 0, ]

##cd4
sceT4 = as(se0[,intersect(grep(se0$annotation,pattern = "T"),grep(se0$annotation,pattern = "8",invert = T))], "SingleCellExperiment")
rowData(sceT4) = NULL

sceT4 <- addPerFeatureQC(sceT4)
sceT4 = sceT4[which(rowData(sceT4)$detected>25),]
sceT4 = sceT4[which(rownames(sceT4) %in% gtf_df$gene_id),]


scaT4 = SceToSingleCellAssay(sceT4, class = "SingleCellAssay")

colData(scaT4)$detected = (colData(scaT4)$detected) # detected == n_gene (CDR)
colData(scaT4)$mouse= as.factor(colData(scaT4)$mouse)
colData(scaT4)$annotation= as.factor(colData(scaT4)$annotation)
scaT4 = scaT4[rowSums(assay(scaT4)) != 0, ]

rm(sceT,sceB,sceT8,sceT4,sce)


summary(lm(colData(scaB)$nFeature_RNA ~ colData(scaB)$nCount_RNA + colData(scaB)$age_weeks.group + colData(scaB)$predictedAge  + colData(scaB)$hbProportion + as.factor(colData(scaB)$mouse)))
resB = residuals(lm(colData(scaB)$nFeature_RNA ~ colData(scaB)$nCount_RNA + colData(scaB)$hbProportion))
boxplot(colData(scaB)$nFeature_RNA ~ colData(scaB)$predictedAge)
boxplot(resB ~ colData(scaB)$predictedAge)

summary(lm(colData(scaT)$nFeature_RNA ~ colData(scaT)$nCount_RNA + colData(scaT)$predictedAge+ colData(scaT)$hbProportion + as.factor(colData(scaT)$mouse)))
resT = residuals(lm(colData(scaT)$nFeature_RNA ~ colData(scaT)$nCount_RNA + colData(scaT)$hbProportion))
boxplot(colData(scaT)$nFeature_RNA ~ colData(scaT)$predictedAge)
boxplot(resT ~ colData(scaT)$predictedAge)

summary(lm(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$nCount_RNA + colData(scaT4)$predictedAge+ colData(scaT4)$hbProportion + as.factor(colData(scaT4)$mouse)))
resT4 = residuals(lm(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$nCount_RNA + colData(scaT4)$hbProportion))
boxplot(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$predictedAge)
boxplot(resT4 ~ colData(scaT4)$predictedAge)

summary(lm(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$nCount_RNA + colData(scaT8)$predictedAge+ colData(scaT8)$hbProportion + as.factor(colData(scaT8)$mouse)))
resT8 = residuals(lm(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$nCount_RNA + colData(scaT8)$hbProportion))
boxplot(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$predictedAge)
boxplot(resT8 ~ colData(scaT8)$predictedAge)


summary(lm(colData(sca)$nFeature_RNA ~ colData(sca)$nCount_RNA + colData(sca)$predictedAge+ colData(sca)$hbProportion + as.factor(colData(sca)$mouse)+ as.factor(colData(sca)$annotation)))
res = residuals(lm(colData(sca)$nFeature_RNA ~ colData(sca)$nCount_RNA + colData(sca)$hbProportion+ as.factor(colData(sca)$annotation)))
boxplot(colData(sca)$nFeature_RNA ~ colData(sca)$predictedAge)
boxplot(res ~ colData(sca)$predictedAge)


##In 101, 77, 36, 10
summary(lm(colData(scaB)$nFeature_RNA[colData(scaB)$age_weeks.group==101] ~ colData(scaB)$nCount_RNA[colData(scaB)$age_weeks.group==101]  + colData(scaB)$predictedAge[colData(scaB)$age_weeks.group==101]  + colData(scaB)$hbProportion[colData(scaB)$age_weeks.group==101]))
summary(lm(colData(scaB)$nFeature_RNA[colData(scaB)$age_weeks.group==77] ~ colData(scaB)$nCount_RNA[colData(scaB)$age_weeks.group==77]  + colData(scaB)$predictedAge[colData(scaB)$age_weeks.group==77]  + colData(scaB)$hbProportion[colData(scaB)$age_weeks.group==77]))
summary(lm(colData(scaB)$nFeature_RNA[colData(scaB)$age_weeks.group==36] ~ colData(scaB)$nCount_RNA[colData(scaB)$age_weeks.group==36]  + colData(scaB)$predictedAge[colData(scaB)$age_weeks.group==36]  + colData(scaB)$hbProportion[colData(scaB)$age_weeks.group==36]))
summary(lm(colData(scaB)$nFeature_RNA[colData(scaB)$age_weeks.group==10] ~ colData(scaB)$nCount_RNA[colData(scaB)$age_weeks.group==10]  + colData(scaB)$predictedAge[colData(scaB)$age_weeks.group==10]  + colData(scaB)$hbProportion[colData(scaB)$age_weeks.group==10]))

