setwd("D:/OnlineFolders/BitSync/CurrentWork/sc_epigen_clock/")

library(speckle)
library(limma)
library(ggplot2)
library(SummarizedExperiment)
library(scater)

##################
## Example data ##
##################

# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
cell_info <- speckle_example_data()
head(cell_info)

# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = cell_info$clusters, sample = cell_info$samples, 
          group = cell_info$group)

# Plot cell type proportions
plotCellTypeProps(clusters=cell_info$clusters, sample=cell_info$samples)


##################################
##  Real data group based test  ##
##################################

load("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/RDataFiles/QCpassedCells.RDa")
colDataCells = read.delim("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/MethylationAging_Shared/Expression/DifferentialExpression/annotation_meta_data_inc_meth_V2.tab")

#Add fine cell type annoataion.
colDataCells["fineAnnotation"] = colDataCells$annotation
colDataCells$annotation[grep("T",colDataCells$annotation)] = "TCell"

se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(dat.counts.blood.QCpassed)), colData=colDataCells)
se0 <- addPerCellQC(se0)
se0 <- addPerFeatureQC(se0)
se0 <- scater::logNormCounts(se0)

# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = se0$annotation, sample = se0$mouse, group = se0$age_weeks.group)

# Plot cell type proportions
plotCellTypeProps(clusters = se0$annotation, sample = se0$mouse)

######################
## Complex modeling ##
######################

prop.logit <- getTransformedProps(clusters = se0$annotation, sample=se0$mouse, transform = "logit")

ageNumeric = colDataCells$age_weeks.group[match(colnames(prop.logit$Counts),colDataCells$mouse)]
ageFactor = as.factor(ageNumeric)
ageGrouped = ageNumeric
ageGrouped[which(ageGrouped==101)] = 77
ageGrouped = as.factor(ageGrouped)

design.anova <- model.matrix(~0+ageFactor)
design.anova

des.dose <- model.matrix(~ageNumeric)
des.dose

##Same as above.
propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3,4), robust=TRUE, trend = FALSE, sort=TRUE)
write.table(propeller.anova(prop.logit,design = design.anova, coef=c(1,2,3,4), robust=TRUE, trend = FALSE, sort=TRUE), file = "./cellProportionEffects.factor.txt",quote = F,sep = "\t")

propeller.anova(prop.logit,design = des.dose, coef=c(1,2), robust=TRUE, trend = FALSE, sort=TRUE)
##Alt for the line above here.
# fit <- lmFit(prop.logit$TransformedProps,des.dose)
# fit <- eBayes(fit, robust=TRUE)
# topTable(fit,coef=2)
write.table(propeller.anova(prop.logit,design = des.dose, coef=c(1,2), robust=TRUE, trend = FALSE, sort=TRUE), file = "./cellProportionEffects.continous.txt",quote = F,sep = "\t")
# Plot cell type proportions
plotCellTypeProps(clusters = se0$annotation, sample = se0$mouse)

