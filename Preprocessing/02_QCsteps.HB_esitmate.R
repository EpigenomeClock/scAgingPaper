setwd("C:/Users/admm414r/Google Drive/EBI/sc_Aging/")
setwd("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/")

load("./RDataFiles/QCpassedCells.RDa")

# dat.counts.blood.QCpassed.2 = dat.counts.blood.QCpassed[which(rowSums(dat.counts.blood.QCpassed)!=0),]
# write.table(dat.counts.blood.QCpassed.2,sep="\t",quote=F, file = "./RDataFiles/sc_MT_Babraham_Blood_aging_QCpassedCells.txt")

#"ENSMUSG00000073940" #Hbb-bt
#"ENSMUSG00000052305" # Hbb-bs
#"ENSMUSG00000069919" #Hba-a1
#"ENSMUSG00000069917" #Hba-a2


hbGenes = c("ENSMUSG00000073940", "ENSMUSG00000052305", "ENSMUSG00000069919", "ENSMUSG00000069917")

sce <- SingleCellExperiment::SingleCellExperiment(list(counts=dat.counts.blood.QCpassed))

sce <- scater::logNormCounts(sce)

sce.QCpassed.hb = sce[which(rownames(sce) %in% hbGenes),]
colSums(SingleCellExperiment::logcounts(sce.QCpassed.hb))
vioplot::vioplot(SingleCellExperiment::logcounts(sce.QCpassed.hb)[1,])
vioplot::vioplot(SingleCellExperiment::logcounts(sce.QCpassed.hb)[2,])
vioplot::vioplot(SingleCellExperiment::logcounts(sce.QCpassed.hb)[3,])
vioplot::vioplot(SingleCellExperiment::logcounts(sce.QCpassed.hb)[4,])

hist(colSums(SingleCellExperiment::logcounts(sce.QCpassed.hb)))
hist(colSums(SingleCellExperiment::logcounts(sce.QCpassed.hb)[c(1,4),]))
hist(colSums(SingleCellExperiment::logcounts(sce.QCpassed.hb)[c(2:3),]))

hist(SingleCellExperiment::logcounts(sce.QCpassed.hb)[1,])
g1Annot = rep(0,ncol(sce.QCpassed.hb))
g1Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[1,]>0)]=1
g1Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[1,]>10)]=2

hist(SingleCellExperiment::logcounts(sce.QCpassed.hb)[2,])
g2Annot = rep(0,ncol(sce.QCpassed.hb))
g2Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[2,]>0)]=1
g2Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[2,]>12)]=2

hist(SingleCellExperiment::logcounts(sce.QCpassed.hb)[3,])
g3Annot = rep(0,ncol(sce.QCpassed.hb))
g3Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[3,]>0)]=1
g3Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[3,]>12)]=2

hist(SingleCellExperiment::logcounts(sce.QCpassed.hb)[4,])
g4Annot = rep(0,ncol(sce.QCpassed.hb))
g4Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[4,]>0)]=1
g4Annot[which(SingleCellExperiment::logcounts(sce.QCpassed.hb)[4,]>10)]=2

annotHbGenes = cbind(g1Annot,g2Annot,g3Annot,g4Annot)
hist(rowSums(annotHbGenes))
hist(ceiling(rowMeans(annotHbGenes)))
