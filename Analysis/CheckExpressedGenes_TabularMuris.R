setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Tabular_muris_Senis/")

library(MAST)
library(scater)
library(SummarizedExperiment)
library(data.table)

##Do all tissues

ds <-  c("./facs_RDS/facs.normalized.Aorta.rds","./facs_RDS/facs.normalized.BAT.rds","./facs_RDS/facs.normalized.Bladder.rds","./facs_RDS/facs.normalized.Brain_Myeloid.rds","./facs_RDS/facs.normalized.Diaphragm.rds","./facs_RDS/facs.normalized.GAT.rds","./facs_RDS/facs.normalized.Heart.rds","./facs_RDS/facs.normalized.Kidney.rds",
         "./facs_RDS/facs.normalized.Large_Intestine.rds","./facs_RDS/facs.normalized.Limb_Muscle.rds","./facs_RDS/facs.normalized.Liver.rds","./facs_RDS/facs.normalized.Lung.rds","./facs_RDS/facs.normalized.Mammary_Gland.rds","./facs_RDS/facs.normalized.Marrow.rds","./facs_RDS/facs.normalized.MAT.rds","./facs_RDS/facs.normalized.Pancreas.rds",
         "./facs_RDS/facs.normalized.SCAT.rds","./facs_RDS/facs.normalized.Skin.rds","./facs_RDS/facs.normalized.Spleen.rds","./facs_RDS/facs.normalized.Thymus.rds","./facs_RDS/facs.normalized.Trachea.rds")

#ds <-  c("./facs_RDS/facs.normalized.Kidney.rds","./facs_RDS/facs.normalized.Large_Intestine.rds","./facs_RDS/facs.normalized.Liver.rds","./facs_RDS/facs.normalized.Lung.rds","./facs_RDS/facs.normalized.Mammary_Gland.rds""./facs_RDS/facs.normalized.Spleen.rds","./facs_RDS/facs.normalized.Skin.rds","./facs_RDS/facs.normalized.Thymus.rds")


gtf<-read.delim("../Mus_musculus.GRCm38.96.gtf.gz",as.is=T,skip=5,header=F)
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


#Spleen, Thymus
for(d in ds){
  
  data = readRDS(d)
  data = data[which(rownames(data) %in% gtf_df$name),]
  colData(data)["n_protCodingGenes"] = colSums(assay(data)>0)
  
  data <- SceToSingleCellAssay(data, class = "SingleCellAssay")
  
  print(gsub("./facs_RDS/facs.normalized.","",d))
  #boxplot(colData(data)$n_genes ~ colData(data)$age_num,main=gsub("./facs_RDS/facs.normalized.","",d))
  boxplot(colData(data)$n_protCodingGenes ~ colData(data)$age_num,main=gsub("./facs_RDS/facs.normalized.","",d))
  
  if(length(levels(as.factor(colData(data)$sex)))!=1){
    lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$age_num + as.factor(colData(data)$sex)+ as.factor(colData(data)$cell_ontology_class)))
    print(coefficients(lmout)[3,])
    
    lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))
    
    boxplot(lmout ~ colData(data)$age_num,main=gsub("./facs_RDS/facs.normalized.","corrected ",d))
  } else {
    lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$age_num + as.factor(colData(data)$cell_ontology_class)))
    print(coefficients(lmout)[3,])
    
    lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))
    
    boxplot(lmout ~ colData(data)$age_num,main=gsub("./facs_RDS/facs.normalized.","corrected ",d))
  }
}


##binary (between 77 weeks and 104 weeks (~29 months)) [no data at 29, testing 24]
cutOff = 24
#cutOff = 29
##At 24, Thymus,spleen, Liver sig in right direction.
for(d in ds){
  
  data = readRDS(d)
  print(d)
  # print(table(colData(data)$age_num))
  
  data = data[which(rownames(data) %in% gtf_df$name),]
  colData(data)["n_protCodingGenes"] = colSums(assay(data)>0)
  colData(data)["binaryAge"] = colData(data)$age_num
  colData(data)$binaryAge[which(colData(data)$age_num<cutOff)] <- "Y"
  colData(data)$binaryAge[which(colData(data)$age_num>=cutOff)] <- "O"

  print(length(unique(colData(data)$binaryAge)))
  print(table(unique(colData(data)$age_num)))

  if(length(unique(colData(data)$binaryAge))==2){
    data <- SceToSingleCellAssay(data, class = "SingleCellAssay")

    print(gsub("./facs_RDS/facs.normalized.","",d))
    #boxplot(colData(data)$n_genes ~ colData(data)$age_num,main=gsub("./facs_RDS/facs.normalized.","",d))

    boxplot(colData(data)$n_protCodingGenes ~ colData(data)$binaryAge,main=gsub("./facs_RDS/facs.normalized.","",d))

    if(length(levels(as.factor(colData(data)$sex)))!=1){
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$binaryAge + as.factor(colData(data)$sex)+ as.factor(colData(data)$cell_ontology_class)))
      print(lmout)

      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))

      boxplot(lmout ~ colData(data)$binaryAge,main=gsub("./facs_RDS/facs.normalized.","corrected ",d))
    } else {
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$binaryAge + as.factor(colData(data)$cell_ontology_class)))
      print(lmout)

      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))

      boxplot(lmout ~ colData(data)$binaryAge,main=gsub("./facs_RDS/facs.normalized.","corrected ",d))
    }
  }
}

# 
# Aorta <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Aorta.rds"), class = "SingleCellAssay")
# BAT <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.BAT.rds"), class = "SingleCellAssay")
# Bladder <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Bladder.rds"), class = "SingleCellAssay")
# Brain_Myeloid <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Brain_Myeloid.rds"), class = "SingleCellAssay")
# Diaphragm <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Diaphragm.rds"), class = "SingleCellAssay")
# GAT <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.GAT.rds"), class = "SingleCellAssay")
# Heart <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Heart.rds"), class = "SingleCellAssay")
# Kidney <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Kidney.rds"), class = "SingleCellAssay")
# Large_Intestine <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Large_Intestine.rds"), class = "SingleCellAssay")
# Limb_Muscle <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Limb_Muscle.rds"), class = "SingleCellAssay")
# Liver <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Liver.rds"), class = "SingleCellAssay")
# Lung <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Lung.rds"), class = "SingleCellAssay")
# Mammary_Gland <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Mammary_Gland.rds"), class = "SingleCellAssay")
# Marrow <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Marrow.rds"), class = "SingleCellAssay")
# MAT <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.MAT.rds"), class = "SingleCellAssay")
# Pancreas <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Pancreas.rds"), class = "SingleCellAssay")
# SCAT <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.SCAT.rds"), class = "SingleCellAssay")
# Skin <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Skin.rds"), class = "SingleCellAssay")
# Spleen <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Spleen.rds"), class = "SingleCellAssay")
# Thymus <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Thymus.rds"), class = "SingleCellAssay")
# Trachea <- SceToSingleCellAssay(readRDS("./facs_RDS/facs.normalized.Trachea.rds"), class = "SingleCellAssay")
# 
# 
# cbind()
# 
# #colData(Aorta)$detected = (colData(Aorta)$detected) # detected == n_gene (CDR)
# # colData(Aorta)$mouse= as.factor(colData(Aorta)$mouse)
# # colData(Aorta)$annotation= as.factor(colData(Aorta)$annotation)
# Aorta = Aorta[rowSums(assay(Aorta)) != 0, ]
# 
# 
# summary(lm(colData(Aorta)$n_genes ~ colData(Aorta)$n_counts + colData(Aorta)$age_num + as.factor(colData(Aorta)$mouse.id) + as.factor(colData(Aorta)$cell_ontology_class)))
# boxplot(colData(Aorta)$n_genes ~ colData(Aorta)$age_num)
# # resB = residuals(lm(colData(scaB)$nFeature_RNA ~ colData(scaB)$nCount_RNA + colData(scaB)$hbProportion))
# # boxplot(resB ~ colData(scaB)$age_weeks.group)
# # 
# 
# summary(lm(colData(BAT)$n_genes ~ colData(BAT)$n_counts + colData(BAT)$age_num + as.factor(colData(BAT)$mouse.id) + as.factor(colData(BAT)$cell_ontology_class)))
# boxplot(colData(BAT)$n_genes ~ colData(BAT)$age_num)
# # resB = residuals(lm(colData(scaB)$nFeature_RNA ~ colData(scaB)$nCount_RNA + colData(scaB)$hbProportion))
# # boxplot(resB ~ colData(scaB)$age_weeks.group)
# 
# boxplot(colData(BAT)$n_genes[which(BAT$cell_ontology_class=="T cell")] ~ colData(BAT)$age_num[which(BAT$cell_ontology_class=="T cell")])
# boxplot(colData(BAT)$n_genes[which(BAT$cell_ontology_class=="B cell")] ~ colData(BAT)$age_num[which(BAT$cell_ontology_class=="B cell")])
# 
# # summary(lm(colData(scaT)$nFeature_RNA ~ colData(scaT)$nCount_RNA + colData(scaT)$age_weeks+ colData(scaT)$hbProportion + as.factor(colData(scaT)$mouse)))
# # resT = residuals(lm(colData(scaT)$nFeature_RNA ~ colData(scaT)$nCount_RNA + colData(scaT)$hbProportion))
# # boxplot(colData(scaT)$nFeature_RNA ~ colData(scaT)$age_weeks.group)
# # boxplot(resT ~ colData(scaT)$age_weeks.group)
# # 
# # summary(lm(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$nCount_RNA + colData(scaT4)$age_weeks+ colData(scaT4)$hbProportion + as.factor(colData(scaT4)$mouse)))
# # resT4 = residuals(lm(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$nCount_RNA + colData(scaT4)$hbProportion))
# # boxplot(colData(scaT4)$nFeature_RNA ~ colData(scaT4)$age_weeks.group)
# # boxplot(resT4 ~ colData(scaT4)$age_weeks.group)
# # 
# # summary(lm(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$nCount_RNA + colData(scaT8)$age_weeks+ colData(scaT8)$hbProportion + as.factor(colData(scaT8)$mouse)))
# # resT8 = residuals(lm(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$nCount_RNA + colData(scaT8)$hbProportion))
# # boxplot(colData(scaT8)$nFeature_RNA ~ colData(scaT8)$age_weeks.group)
# # boxplot(resT8 ~ colData(scaT8)$age_weeks.group)
# # 
# # 
# # summary(lm(colData(sca)$nFeature_RNA ~ colData(sca)$nCount_RNA + colData(sca)$age_weeks+ colData(sca)$hbProportion + as.factor(colData(sca)$mouse)+ as.factor(colData(sca)$annotation)))
# # res = residuals(lm(colData(sca)$nFeature_RNA ~ colData(sca)$nCount_RNA + colData(sca)$hbProportion+ as.factor(colData(sca)$annotation)))
# # boxplot(colData(sca)$nFeature_RNA ~ colData(sca)$age_weeks.group)
# # boxplot(res ~ colData(sca)$age_weeks.group)
