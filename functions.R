###functions file


correctedCPM <- function(dataframe){
  ###Calculates corrected CPM value
  
  #compute sum factors
  mat <- as.matrix(dataframe)
  
  sce <- SingleCellExperiment(list(counts=mat))
  sce <- computeSumFactors(sce)
  #summary(sizeFactors(sce))
  
  ####convert and get factors
  converted <- convertTo(sce, type="edgeR")
  df1 <- converted$samples
  norm.factors <- df1$norm.factors
  
  #cpm normalization
  cpm <- dataframe
  for (i in seq(1,ncol(dataframe))){
    cs <- colSums(dataframe)[i]
    cpm[,i] <- dataframe[,i] / cs
  }
  cpm <- cpm*1E6
  
  ####corr CPM with sumfactors
  DT_corr.cpm <- as.data.frame(matrix(ncol = ncol(dataframe),nrow = nrow(dataframe)))
  
  rownames(DT_corr.cpm) <- rownames(dataframe)
  colnames(DT_corr.cpm) <- colnames(dataframe)
  for (i in seq(1,ncol(cpm))){
    DT_corr.cpm[,i] <- cpm[,i] * norm.factors[i]
  }
  
  return(DT_corr.cpm)
  #return(cpm)
}


log_xp1 <- function(x){
  value <- log(x+1)
  return(value)
}



assign.GeneNames <- function(EnsemblIDs,target.DF){
  ####get gene_ids from ensembl ids
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- (EnsemblIDs)
  gene_ids = getBM(
    values = genes,
    filters = c("ensembl_gene_id"),
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    mart = mouse
  )
  idx <- sapply(EnsemblIDs, function(x) {
    which(gene_ids$ensembl_gene_id == x)
  })
  
  #only keep ensemblIDs which were assigned correctly
  #target.DF <- target.DF[which(EnsemblIDs%in%gene_ids$ensembl_gene_id),]
  #target.DF <- as.data.frame(target.DF)
  
  target.DF$gene_id <- gene_ids$external_gene_name[match(EnsemblIDs,gene_ids$ensembl_gene_id)]
  target.DF$gene_desc <- gene_ids$description[match(EnsemblIDs,gene_ids$ensembl_gene_id)]
  
  return(target.DF)
  
}



weighted.means_celltype_mouse <- function(mouse.list,current.celltype.ids,count.matrix,all_genes_rawCount.matrix,remove.zeroes = FALSE,suppress.print = FALSE){
  ###given a set of sample ids (eg. for a specific celltype), calculate the weighted mean expression per gene per mouse
  ###weights represent the total number of counts per sample
  
  weighted.means.mouse <- list()
  variances.mouse <- list()
  readdepth.sum.mouse <- list()
  
  ###iterate over list of mice and extract sample ids by mice
  for(m in seq(1,length(mouse.list))){
    curr.mouse.name <- names(mouse.list)[m]
    curr.mouse.ids <- mouse.list[[m]]
    
    ###overlapping ids from current celltype and current mouse
    curr.ids_mouse_ct <- current.celltype.ids[current.celltype.ids%in%curr.mouse.ids]
    curr.mouse.readdepth <- sum(colSums(as.data.frame(all_genes_rawCount.matrix[,colnames(all_genes_rawCount.matrix)%in%curr.ids_mouse_ct])))
    curr.weights <- dat.info$nCount_RNA[rownames(dat.info)%in%curr.ids_mouse_ct]
    
    ###get expression data for current subset of sample ids
    curr.data <- as.data.frame(count.matrix[,colnames(count.matrix)%in%curr.ids_mouse_ct])
    rownames(curr.data) <- rownames(count.matrix)
    
    ###calculate weighted means per gene per mouse and variance per gene per mouse, extract readdepth per sample
    curr.weighted.means <- apply(curr.data, 1, function(x) weighted.mean(x, curr.weights))
    curr.variance <- apply(curr.data, 1,var)
    #curr.mouse.readdepth <- sum(dat.info$nCount_RNA[match(colnames(curr.data),rownames(dat.info))])
    
    ###store weighted means and variance 
    weighted.means.mouse[[m]] <- curr.weighted.means
    variances.mouse[[m]] <- curr.variance
    readdepth.sum.mouse[[m]] <- curr.mouse.readdepth
    
    names(weighted.means.mouse)[m] <- paste("Mouse",curr.mouse.name,sep = " ")
    names(variances.mouse)[m] <- paste("Mouse",curr.mouse.name,sep = " ")
    names(readdepth.sum.mouse)[m] <- paste("Mouse",curr.mouse.name,sep = " ")
    
    if(suppress.print == FALSE){
      print(paste("Mouse Nr.",curr.mouse.name))
    }
  }
  
  ########Create matrix with weighted means and variance per mouse#########
  weighted.means.df <- data.frame(matrix(nrow = length(weighted.means.mouse$`Mouse 23`),ncol = length(weighted.means.mouse)))
  rownames(weighted.means.df) <- names(weighted.means.mouse$`Mouse 23`)
  colnames(weighted.means.df) <- names(weighted.means.mouse)
  variance.df <- weighted.means.df
  readdepth.df <- data.frame(matrix(nrow = 1,ncol = length(readdepth.sum.mouse)),row.names = "readdepth")
  colnames(readdepth.df) <- names(readdepth.sum.mouse)
  
  for(l in seq(1,length(weighted.means.mouse))){
    curr.means.mouse <- weighted.means.mouse[[l]]
    curr.var.mouse <- variances.mouse[[l]]
    
    weighted.means.df[,l] <- curr.means.mouse[match(rownames(weighted.means.df),names(curr.means.mouse))]
    variance.df[,l] <- curr.var.mouse[match(rownames(variance.df),names(curr.var.mouse))]
    
  }
  for(r in seq(1,length(readdepth.sum.mouse))){
    readdepth.df[1,r] <- readdepth.sum.mouse[[r]]
  }
  
  ###remove all mice in which the current celltype is not present at all
  weighted.means.df <- weighted.means.df[,colSums(is.na(weighted.means.df)) == 0]
  variance.df <- as.data.frame(variance.df[,colSums(is.na(variance.df)) == 0])
  
  result.list <- list(weighted.means.df,variance.df,readdepth.df)
  names(result.list) <- c("weighted.means.df","variance.df","readdepth.df")
  ###check if at least 1 mouse per timepoint is present
  present.mouse.timeponts <- names(sorted.mouse_timepoint.pair)[paste("Mouse",(sorted.mouse_timepoint.pair))%in%colnames(weighted.means.df)]
  if(length(unique(present.mouse.timeponts))!=length(unique(names(sorted.mouse_timepoint.pair)))){
    return("Not at least 1 mouse per timepoint")
  }
  else{
    return(result.list)
  }
}


sorted.Genes_by_avg_CPM <- function(weighted.means.df,sort.by.var = NULL,var.sort.cpm.threshold = NULL){
  ###based on weighted mean expression per mouse, this function calculates the mean expression of every
  ###gene over every mouse and outputs a sorted list from highhest to lowest, based on CPM
  
  
  ###calculate mean and variance expression per gene between mice
  means.per.gene <- apply(weighted.means.df, 1, mean)
  names(means.per.gene) <- rownames(weighted.means.df)
  
  ###also calculate variances between weighted mean gene expression between mice
  var.per.gene <- apply(weighted.means.df, 1, var)
  names(var.per.gene) <- rownames(weighted.means.df)
  
  ###sort genes according to mean expression over all mice and analyze top 1000 genes per celltype
  sorted_means.per.gene <- means.per.gene[order(means.per.gene,decreasing = T)]
  
  ##sort by variance, but consider only genese which are expressed in average >5 CPM
  if(!(is.null(sort.by.var))){
    print("loop1")
    if(!(is.null(var.sort.cpm.threshold))){
      print("loop2")
      avgCPM.thresh <- var.sort.cpm.threshold
    }
    else{avgCPM.thresh <- 5}
    
    var.per.gene <- var.per.gene[means.per.gene > log_xp1(avgCPM.thresh)]
    sorted_var.per.gene <- var.per.gene[order(var.per.gene,decreasing = T)]
    
    return(sorted_var.per.gene)
  }
  else{
    return(sorted_means.per.gene)
  }
  
  
}

linear.model_per.gene <- function(gene.means.mouse,model.for.readdepth = NULL){
  ###builds linear model for every gene of the respective cell ids ,based on the weighted mean expression per mouse
  ###outputs a list in which the lmer element is stored, named according to the gene
  
  gene.means.mouse <- gene.means.mouse
  z <- 1
  lm.list <- list()
  i <- 1
  for(i in seq(1,nrow(gene.means.mouse$weighted.means.df))){
    print(paste("Gene",i))
    
    curr.dat <- as.data.frame(t(gene.means.mouse$weighted.means.df[i,]))
    colnames(curr.dat) <- "expression"
    curr.dat$mouse <- rownames(curr.dat)
    curr.dat$age.group <- dat.info$age_weeks.group[match(curr.dat$mouse,paste("Mouse",dat.info$mouse))]
    curr.dat$readdepth <- as.numeric(gene.means.mouse$readdepth.df[1,])[match(curr.dat$mouse,colnames(gene.means.mouse$readdepth.df))]
    
    if(length(unique(curr.dat$age.group)) < 4){
      print("Not data for mice of all age groups")
      
      lm.list[[z]] <- "No data for mice of all ages"
      names(lm.list)[z] <- rownames(count.data)[i]
      z <- z+1
      next
    }
    
    if(!(is.null(model.for.readdepth))){
      lm <- lm(expression ~ readdepth, data = curr.dat)
    } else{
      lm <- lm(expression ~ age.group, data = curr.dat)
    }
    
    
    lm.list[[z]] <- lm
    names(lm.list)[z] <- rownames(gene.means.mouse$weighted.means.df)[i]
    z <- z+1
  }
  
  return(lm.list)
  
}  


determine.significant.genes_LinearModel <- function(linear.model.list,FDR.threshold,test.N.genes,sorted.genes.vec){
  ###performs multiple testing correction for the p-values of the linear model for every gene
  ###outputs a vector containing all gene names which have a FDR below the FDR threhsold 
  #model.list <- lm.list
  model.list <- linear.model.list
  sorted.genes <- sorted.genes.vec
  ###order genes in model.list, according to the topN genes, for which the model is not singular
  g <- names(model.list)
  
  model.list_sorted <- model.list[order(match(g,names(sorted.genes)))]
  
  
  ####Extract p-values from all models which are not singular
  topNgenes_model.list <- model.list_sorted[1:test.N.genes]
  p_vals <- c()
  for(i in seq(1,length(topNgenes_model.list))){
    print(i)
    if(names(topNgenes_model.list)[i] != ""){
      p_vals[i] <- anova(topNgenes_model.list[[i]])$"Pr(>F)"[1]   
      names(p_vals)[i] <- names(topNgenes_model.list)[i]
    }
  }
  qobj <- qvalue(p = as.numeric(p_vals))
  #hist(qobj$qvalues,breaks = 100)
  sign.genes <- (names(p_vals)[qobj$qvalues < FDR.threshold])
  
  return(sign.genes)
}



rel.freq <- function(df_col){
  table <- as.data.frame(table(df_col))
  table$relFreq <- (table$Freq / sum(table$Freq))*100
  table$totalSum[1] <- sum(table$Freq)
  return(table)
}



Cluster.Abundance.celltypes <- function(correlation_assignment.results,seurat.file_withClusters){
  ###Takes the result of a seurat clustering as input as well as the cell type labels, assigned by correlation assignment
  ###calculates the occupancy of every celltype (by correlation assignment) in the seurat clustering
  
  ###seurat cluster result
  
  seurat.clusters.df <- as.data.frame(seurat.file_withClusters$seurat_clusters)
  colnames(seurat.clusters.df) <- "seurat_cluster"
  
  ###correlation assignment result
  assigned_lineages.dataframe <- correlation_assignment.results
  
  
  
  
  clusters <- unique(seurat.clusters.df$seurat_cluster)
  clusters <- clusters[!is.na(clusters)]
  clusters.cellnames.list <- list()
  for(i in seq(1,length(clusters))){
    subcluster.cells <- rownames(subset(seurat.clusters.df,seurat_cluster == clusters[i]))
    
    clusters.cellnames.list[[i]] <- subcluster.cells
    names(clusters.cellnames.list)[i] <- paste("Cluster",clusters[i])
    
  }
  remove(subcluster.cells)
  
  
  ###Calculate cluster_assignment.matrix
  cluster_celltype_assignmend.matrix <- as.data.frame(matrix(nrow=length(clusters.cellnames.list)))
  #colnames(cluster_celltype_assignmend.matrix) <- c("top_lineage","pct_cells.top_lineage","second.top.lineage","pct_cells.second","n_cells")
  
  
  ###get majority of assigned celltypes by reference samples for cells of each cluster
  for(i in seq(1,length(clusters.cellnames.list))){
    cluster.cellnames <- clusters.cellnames.list[[i]]
    reference_cells <- assigned_lineages.dataframe[which(rownames(assigned_lineages.dataframe)%in%cluster.cellnames),]
    lineage.table <- as.data.frame(rel.freq(reference_cells$cell_lineage))
    
    ordered_lineage.table <- lineage.table[order(lineage.table$relFreq,decreasing = T),]
    
    cluster_celltype_assignmend.matrix$top_lineage[i] <- as.character(ordered_lineage.table$df_col[1])
    cluster_celltype_assignmend.matrix$pct_cells.top_lineage[i] <- as.numeric(ordered_lineage.table$relFreq[1])
    cluster_celltype_assignmend.matrix$second.top.lineage[i] <- as.character(ordered_lineage.table$df_col[2])
    cluster_celltype_assignmend.matrix$pct_cells.second[i] <- as.numeric(ordered_lineage.table$relFreq[2])
    cluster_celltype_assignmend.matrix$third.top.lineage[i] <- as.character(ordered_lineage.table$df_col[3])
    cluster_celltype_assignmend.matrix$pct_cells.third[i] <- as.numeric(ordered_lineage.table$relFreq[3])
    
    cluster_celltype_assignmend.matrix$n_cells[i] <- as.numeric(ordered_lineage.table$totalSum[1])
    
    rownames(cluster_celltype_assignmend.matrix)[i] <- names(clusters.cellnames.list)[i]
    #cluster_celltype_assignmend.matrix[,1] <- NULL
  }
  
  return(cluster_celltype_assignmend.matrix)
  
  
}
