setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/DATA/RNA")

##
file="mt"
#file="all"

# cellSNP.tag.AD.mtx: a file in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles.
# cellSNP.tag.DP.mtx: a file in “Matrix Market exchange formats”, containing the sum of allele depths of the reference and alternative alleles (REF + ALT).
# cellSNP.tag.OTH.mtx: a file in “Matrix Market exchange formats”, containing the sum of allele depths of all the alleles other than REF and ALT.

samples = list.files("./",include.dirs = T,pattern = "EpiAge")

##
cells = NULL
nVar = NULL
nAltVars = NULL
altACount = NULL
allACount = NULL
for(s in samples){
  print(s)
  ad <- read.delim(paste0("./",s,"/",file,"/cellSNP.tag.AD.mtx"),skip=2,header=F)[,c(1,3)]
  oth <- read.delim(paste0("./",s,"/",file,"/cellSNP.tag.OTH.mtx"),skip=2,header=F)[,c(1,3)]
  
  noAlt = F
  if(nrow(ad)!=1){
    ##This means there are alternative alleles.
    ad = ad[-1,]
  } else {
    noAlt=T
  }
  
  if(nrow(oth)!=1){
    #This means there are both alternative and other alternative alleles.
    oth = oth[-1,]
    ad[match(oth[,1],ad[,1]),2] = ad[match(oth[,1],ad[,1]),2]+ oth[,2]
  }
  rm(oth)
  
  dp <- read.delim(paste0("./",s,"/",file,"/cellSNP.tag.DP.mtx"),skip=2,header=F)[,c(1,3)]
  if(nrow(dp)!=1 & noAlt){
    #Skip cell
    next
  }
  cells = c(cells,s)
  allACount = c(allACount,sum(dp[,2]))
  nVar = c(nVar,nrow(dp))
  if(!noAlt){
    altACount = c(altACount,sum(ad[,2]))  
    nAltVars = c(nAltVars,nrow(ad))
  } else {
    altACount = c(altACount,0)
    nAltVars = c(nAltVars,0)
  }
  
  
}


cellInfo = as.data.frame(cbind(cells,nVar,allACount,nAltVars,altACount))

load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
expStats=(cbind(colSums(dat.counts.blood.QCpassed),colSums(dat.counts.blood.QCpassed>0)))
colnames(expStats)=c("sequencedReads","expGenes")
rm(dat.counts.blood.QCpassed)

##annotationInformation.
##Test if epi-age is related to nAlt alleles over nRef alleles, corrected for total nVar
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")
rownames(colDataCells) = colDataCells$X
##Methylation info.
agePredictionScOrg <- read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p1.extended.txt",as.is=T)
agePredictionScOrg2 <- read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p2.extended.txt",as.is=T) 

#Drop non-passing QC cells.
qcInfoP1 = read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd/passQc_Cells.txt",header = F)
qcInfoP2 = read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd/passQc_Cells_p2.txt",header = F)
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

#Drop predictions with less then 40 sites. [Filter]
agePredictionScOrg = agePredictionScOrg[which(agePredictionScOrg$sitesUsed>39),]

agePredictionScOrg["predictedAgeCor"] = agePredictionScOrg$predictedAge - (agePredictionScOrg$medianRandomAge - agePredictionScOrg$actualAge)
colDataCells = cbind(colDataCells,agePredictionScOrg[match(colDataCells$babraham_short, rownames(agePredictionScOrg)),c(4:10)])
colDataCells$sanger_id_rna = gsub(x=colDataCells$sanger_id_rna,replacement = "",pattern = ".cram")

colDataCells = cbind(colDataCells,expStats[match(rownames(colDataCells),rownames(expStats)),])

colDataCellsE = cbind(colDataCells,cellInfo[match(colDataCells$sanger_id_rna,cellInfo$cells),])
colDataCellsE = colDataCellsE[which(!is.na(colDataCellsE$predictedAge)),]

library(lme4)
##Add extra info (Partially specfic for the MTs)
load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
expStats=(cbind(colSums(dat.counts.blood.QCpassed),colSums(dat.counts.blood.QCpassed>0)))
colnames(expStats)=c("sequencedReads","expGenes")

##Calculate expression of MT genes.
if(file=="mt"){
  library(SummarizedExperiment)
  library(scater)
  se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(dat.counts.blood.QCpassed)), colData=colDataCells)
  se0 <- addPerCellQC(se0)
  se0 <- addPerFeatureQC(se0)
  se0 <- scater::logNormCounts(se0)
  
  gtf_df <- as.data.frame(rtracklayer::import('../../Mus_musculus.GRCm38.96.gtf.gz'))
  gtf_df  = gtf_df[which(gtf_df$gene_biotype=="protein_coding"),]
  se0.mt= se0[which(rownames(se0)%in% gtf_df$gene_id[which(gtf_df$seqnames=="MT")]),]
  avgMt = cbind(colMeans(as.matrix(se0.mt@assays@data$logcounts)),colSums(as.matrix(se0.mt@assays@data$counts)>0),colSums(as.matrix(se0.mt@assays@data$counts)))
  colnames(avgMt)=c("avgMtExp","expMtGenes","sequencedMtReads")
}

rm(dat.counts.blood.QCpassed)
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")
rownames(colDataCells) = colDataCells$X

colDataCells["annotationFine"] = colDataCells$annotation
colDataCells$annotation[which(colDataCells$annotationFine %in% c("NveCd8T","MemCD8T"))] = "CD8T"
colDataCells$annotation[which(colDataCells$annotationFine %in% c("EffCD4T","RegT","NveCD4T","MemCD4T"))] = "CD4T"

##Methylation info.
agePredictionScOrg <- read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p1.extended.txt",as.is=T)
agePredictionScOrg2 <- read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd_V5/Blood/Github_Code/preditionScBlood.p2.extended.txt",as.is=T) 


#Drop non-passing QC cells.
qcInfoP1 = read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd/passQc_Cells.txt",header = F)
qcInfoP2 = read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trapp_Re-engineerd/passQc_Cells_p2.txt",header = F)

agePredictionScOrg = agePredictionScOrg[grep(paste(qcInfoP1$V1,collapse = "|"),rownames(agePredictionScOrg)),]
agePredictionScOrg2 = agePredictionScOrg2[grep(paste(qcInfoP2$V1,collapse = "|"),rownames(agePredictionScOrg2)),]


agePredictionSc <- read.delim("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Trap_scAge/NewBloodModel_RealBloodCells_TrapPredicition.txt",as.is=T,sep=" ")
agePredictionSc = agePredictionSc[which(rownames(agePredictionSc) %in% c(rownames(agePredictionScOrg),rownames(agePredictionScOrg2))),]
agePredictionScOrg = agePredictionScOrg[which(rownames(agePredictionScOrg) %in% rownames(agePredictionSc)),]
agePredictionScOrg2 = agePredictionScOrg2[which(rownames(agePredictionScOrg2) %in% rownames(agePredictionSc)),]

all(rownames(agePredictionSc) == c(rownames(agePredictionScOrg),rownames(agePredictionScOrg2)))
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
rownames(agePredictionSc) = rownames(agePredictionScOrg)
rm(agePredictionScOrg,agePredictionScOrg2)

colDataCells = cbind(colDataCells,agePredictionSc[match(colDataCells$babraham_short, rownames(agePredictionSc)),c(3:6)])
colDataCells = cbind(colDataCells,expStats[match(rownames(colDataCells),rownames(expStats)),])
colDataCells$sanger_id_rna = gsub(x=colDataCells$sanger_id_rna,replacement = "",pattern = ".cram")

colDataCellsE = cbind(colDataCells,cellInfo[match(colDataCells$sanger_id_rna,cellInfo$cells),])
if(file=="mt"){
  colDataCellsE = cbind(colDataCellsE,avgMt[match(rownames(colDataCellsE),rownames(avgMt)),])
}
colDataCellsE = colDataCellsE[which(!is.na(colDataCellsE$predictedAge)),]

###

##Quick test
summary(lm(colDataCellsE$altACount~colDataCellsE$predictedAge))
summary(lm(colDataCellsE$nAltVars~colDataCellsE$age_weeks))

varT = t(rbind(as.numeric(colDataCellsE$nAltVars), (as.numeric(colDataCellsE$nVar)-as.numeric(colDataCellsE$nAltVars))))
#summary(glmer(varT ~ colDataCellsE$predictedAge , family = binomial))

summary(glm(varT ~ as.numeric(colDataCellsE$sequencedReads)+as.numeric(colDataCellsE$expGenes)+colDataCellsE$predictedAge+as.numeric(colDataCellsE$age_weeks.group)+as.factor(colDataCellsE$mouse), family="binomial"))
summary(glm(varT ~ as.numeric(colDataCellsE$sequencedReads)+as.numeric(colDataCellsE$expGenes)+colDataCellsE$predictedAge+as.numeric(colDataCellsE$age_weeks)+as.factor(colDataCellsE$annotation), family="binomial"))
summary(glm(varT ~ as.numeric(colDataCellsE$sequencedReads)+as.numeric(colDataCellsE$expGenes)+colDataCellsE$predictedAge+as.numeric(colDataCellsE$age_weeks)+as.factor(colDataCellsE$mouse)+as.factor(colDataCellsE$annotation), family="binomial"))

corrected = residuals(glm(varT ~ as.numeric(colDataCellsE$sequencedReads)+as.numeric(colDataCellsE$expGenes)+as.factor(colDataCellsE$mouse)+as.factor(colDataCellsE$annotation), family="binomial"))

corrected = residuals(glm(varT ~ as.numeric(colDataCellsE$sequencedReads)+as.numeric(colDataCellsE$expGenes)+as.factor(colDataCellsE$annotation), family="binomial"))


cor.test(corrected,as.numeric(colDataCellsE$age_weeks)[-which(is.na(colDataCellsE$annotation))])
cor.test(corrected,as.numeric(colDataCellsE$predictedAge)[-which(is.na(colDataCellsE$annotation))])

##What about within an age?

summary(glm(varT[which(colDataCellsE$age_weeks.group==10),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==10)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==36),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==36)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==77),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==77)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==101),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==101)], family="binomial"))


summary(glm(varT[which(colDataCellsE$age_weeks.group==10),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==10)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==36),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==36)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==77),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==77)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==101),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$sequencedReads)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$expGenes)[which(colDataCellsE$age_weeks.group==101)], family="binomial"))


##Random effect on mouse.
colnames(varT) = c("AltAlleles","refAlleles")
colDataCellsE = cbind(colDataCellsE, varT)
colnames(colDataCellsE)[which(colnames(colDataCellsE)=="age_Week.group")]
#gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
#summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge + age_weeks.group + annotation +expGenes + (1|mouse), data=colDataCellsE , family = binomial))

## summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge + (1|age_weeks.group) + annotation + expGenes + (1|mouse), data=colDataCellsE , family = binomial))

summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge + age_weeks.group + annotation + expGenes + (1|mouse), data=colDataCellsE , family = binomial))


##Final results all
summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge1 + age_weeks.group + scale(nCount_RNA) + (1|mouse) + (1|annotation), data=colDataCellsE , family = binomial))

summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge1 + scale(nCount_RNA) + (1|mouse) + (1|annotation), data=colDataCellsE , family = binomial))

summary(glmer(cbind(AltAlleles,refAlleles) ~ age_weeks + scale(nCount_RNA) + (1|mouse) + (1|annotation), data=colDataCellsE , family = binomial))


##For MT we need to correct for the expression level of MT genes.
summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge + age_weeks.group + annotation + avgMtExp + expMtGenes + (1|mouse), data=colDataCellsE , family = binomial))
#summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge + age_weeks.group + annotation + expGenes + avgMtExp + expMtGenes + (1|mouse), data=colDataCellsE , family = binomial))

summary(glm(varT[which(colDataCellsE$age_weeks.group==10),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==10)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==36),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==36)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==77),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==77)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==101),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$mouse)[which(colDataCellsE$age_weeks.group==101)], family="binomial"))

summary(glm(varT[which(colDataCellsE$age_weeks.group==10),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==10)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==10)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==10)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==36),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==36)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==36)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==36)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==77),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==77)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==77)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==77)], family="binomial"))
summary(glm(varT[which(colDataCellsE$age_weeks.group==101),] ~ colDataCellsE$predictedAge[which(colDataCellsE$age_weeks.group==101)]+as.factor(colDataCellsE$annotation)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$avgMtExp)[which(colDataCellsE$age_weeks.group==101)]+as.numeric(colDataCellsE$expMtGenes)[which(colDataCellsE$age_weeks.group==101)], family="binomial"))


##Final results MT results
summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge1 + age_weeks.group + scale(sequencedMtReads) + (1|mouse) + (1|annotation), data=colDataCellsE , family = binomial))
summary(glmer(cbind(AltAlleles,refAlleles) ~ predictedAge1 + scale(sequencedMtReads) + (1|mouse)  + (1|annotation), data=colDataCellsE , family = binomial))
summary(glmer(cbind(AltAlleles,refAlleles) ~ age_weeks + scale(sequencedMtReads) + (1|mouse)  + (1|annotation), data=colDataCellsE , family = binomial))
