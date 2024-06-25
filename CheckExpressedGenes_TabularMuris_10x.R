setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/Tabular_muris_Senis/")

library(MAST)
library(scater)
library(SummarizedExperiment)
library(data.table)

##Do all tissues

ds <-  list.files("./droplet_RDS/",pattern = ".rds",full.names = T)

ds = ds[c(2,4,5,7,8,9,11,12)]

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

for(d in ds){
  data = readRDS(d)
  
  if(length(unique(colData(data)$age_num))>1){
    
    data <- SceToSingleCellAssay(data[which(rownames(data) %in% gtf_df$name),], class = "SingleCellAssay")
    colData(data)["n_protCodingGenes"] = colSums(assay(data)>0)
    
    print(gsub("./droplet_RDS/droplet.normalized.","",d))
    boxplot(colData(data)$n_protCodingGenes ~ colData(data)$age_num,main=gsub("./droplet_RDS/droplet.normalized.","",d))
    
    if(length(levels(as.factor(colData(data)$sex)))!=1){
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$age_num + as.factor(colData(data)$sex) + as.factor(colData(data)$cell_ontology_class)))
      print(coefficients(lmout)[3,])
      
      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$sex) + as.factor(colData(data)$cell_ontology_class)))
      
      boxplot(lmout ~ colData(data)$age_num,main=gsub("./droplet_RDS/droplet.normalized.","corrected ",d))
    } else {
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$age_num + as.factor(colData(data)$cell_ontology_class)))
      print(coefficients(lmout)[3,])
      
      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))
      
      boxplot(lmout ~ colData(data)$age_num,main=gsub("./droplet_RDS/droplet.normalized.","corrected ",d))
    }
  }
  
}


##Binary:
#cutOff = 20
cutOff = 24
#cutOff = 29
drop=21


for(d in ds){
  data = readRDS(d)
  print(d)
  print(table(colData(data)$age_num))
  colData(data)["binaryAge"] = colData(data)$age_num

  colData(data)$binaryAge[which(colData(data)$age_num<cutOff)] <- "Y"
  colData(data)$binaryAge[which(colData(data)$age_num>=cutOff)] <- "O"
  colData(data)$binaryAge[which(colData(data)$age_num %in% drop)] <- NA

  if(length(unique(na.omit(colData(data)$binaryAge)))>1){

    data <- SceToSingleCellAssay(data[which(rownames(data) %in% gtf_df$name),], class = "SingleCellAssay")
    colData(data)["n_protCodingGenes"] = colSums(assay(data)>0)

    print(gsub("./droplet_RDS/droplet.normalized.","",d))
    boxplot(colData(data)$n_protCodingGenes ~ colData(data)$age_num,main=gsub("./droplet_RDS/droplet.normalized.","",d))

    if(length(levels(as.factor(colData(data)$sex)))!=1){
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$binaryAge + as.factor(colData(data)$sex) + as.factor(colData(data)$cell_ontology_class)))
      print((lmout))

      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$sex) + as.factor(colData(data)$cell_ontology_class)))

      boxplot(lmout ~ colData(data)$binaryAge,main=gsub("./droplet_RDS/droplet.normalized.","corrected ",d))
    } else {
      lmout = summary(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + colData(data)$binaryAge + as.factor(colData(data)$cell_ontology_class)))
      print((lmout))

      lmout = residuals(lm(colData(data)$n_protCodingGenes ~ colData(data)$n_counts + as.factor(colData(data)$cell_ontology_class)))

      boxplot(lmout ~ colData(data)$binaryAge,main=gsub("./droplet_RDS/droplet.normalized.","corrected ",d))
    }
  }
  
}
