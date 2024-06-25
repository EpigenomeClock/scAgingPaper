###Load the count matrix and the feature information from the compressed folder into R

#setwd("~/Google Drive/sc_Aging/")
setwd("D:/OnlineFolders/Google_Drive/EBI/sc_Aging/")

library(stringr)



#############Prepare Count Matrix and Feature information matrix##############
###Load Count Matrix and Feature information file
dat.counts <- as.data.frame(read.table("AgingSC.102019.featureCounts.genes.counts.unique.tsv_counts.tsv.gz",sep = '\t', header = TRUE))
dat.features <- as.data.frame(read.table("AgingSC.102019.featureCounts.genes.counts.unique.tsv_features.tsv.gz",sep = '\t', header = TRUE))

###Assign Rownames to gene names
rownames(dat.counts) <- dat.counts$X
dat.counts$X <- NULL

###Reduce column names of count matrix to sample numbers
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\d*")
} 
for (i in seq(1,ncol(dat.counts))){
  colnames(dat.counts)[i] <- numextract(colnames(dat.counts)[i])
  
}

###Remove rows which are no gene count information but summary statistics
dat.counts.all <- dat.counts
dat.counts <- dat.counts.all[(substr(rownames(dat.counts.all),1,4) == "ENSM"),]
summary_statistics <- dat.counts.all[(substr(rownames(dat.counts.all),1,4) != "ENSM"),]



###Import information of birth date and death date of mice and calculate age
dat.ages <- read.table("AgeInfoData_sorted1.txt",sep = "\t",header = TRUE)
col <- ncol(dat.ages)+1
for (i in seq(1,nrow(dat.ages))){
  DOB <- dat.ages[i,which(colnames(dat.ages)=="DOB")]
  DOD <- dat.ages[i,which(colnames(dat.ages)=="DOD")]
  lifetime_weeks <- round(as.double(difftime(strptime(DOD, format = "%d.%m.%Y"),
                                             strptime(DOB, format = "%d.%m.%Y"),units="weeks")),digits = 1)
  dat.ages[i,col] <- lifetime_weeks
}
colnames(dat.ages)[ncol(dat.ages)] <- "age_weeks"


###Import Information of Measurment and link between sampleID and mouse
dat.info <- as.data.frame(read.table("2019-05-29_filenames_aging_rna-batch01.tsv"))
dat.info <- dat.info[-1,]

#Export information of mouse from Babraham_Id
for(i in (seq(1,nrow(dat.info)))){
  sample <- as.numeric(substr(dat.info[i,2],start = 2,stop = 3))
  dat.info[i,4] <- sample
}

#Export sampleID from SangerID
for (i in seq(1,nrow(dat.info))){
  rownames(dat.info)[i] <- numextract(dat.info[i,1])
  
}
colnames(dat.info) <- c("sanger_id","babraham_id","tissue","mouse")

#Add information about age of mice
dat.info$age_weeks <- dat.ages$age_weeks[match(dat.info$mouse,dat.ages$Mouse)]

#assign age groups
unique(dat.info$age_weeks)
dat.info$age_weeks.group[dat.info$age_weeks == 9.6] <- 10
dat.info$age_weeks.group[dat.info$age_weeks == 9.7] <- 10
dat.info$age_weeks.group[dat.info$age_weeks == 36.1] <- 36
dat.info$age_weeks.group[dat.info$age_weeks == 77.0] <- 77
dat.info$age_weeks.group[dat.info$age_weeks == 77.4] <- 77
dat.info$age_weeks.group[dat.info$age_weeks == 99.1] <- 101
dat.info$age_weeks.group[dat.info$age_weeks == 101.1] <- 101


############Save CountMatrix, MetaData and CountStatistics######
#define savefolder

save(dat.counts,file = "RDataFiles/scCountMatrix.RDa")
save(dat.ages,dat.info,dat.features,file = "RDataFiles/scMetaData.Rda")
save(summary_statistic,file = "RDataFiles/scCountStatistics.RDa")






