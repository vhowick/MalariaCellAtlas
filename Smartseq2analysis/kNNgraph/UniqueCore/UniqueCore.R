library(scater)
library(plotly)
library(pheatmap)
library(Rtsne)
library(viridis)
library(VennDiagram)
setwd("/Users/vh3/Documents/MCA/ANALYSIS_3")
rawmca <- readRDS("MCAqcTMMSLS34_20181026.rds")

pheno <- as.data.frame(colData(rawmca))
library(plyr)
test <- ddply(pheno,.(ShortenedLifeStage4),function(x) x[sample(nrow(x),60),])

cells <- test$sample_id

datalist = list()

for (i in unique(pheno$ShortenedLifeStage4)) {
  stage <- rawmca[, rawmca$ShortenedLifeStage4==i]
  substage <- stage[, stage$sample_id %in% cells]
  counts <- as.data.frame(counts(substage))
  Ooc <- apply(counts, MARGIN = 1, function(x) sum(x > 0))
  Ooc <- as.data.frame(Ooc)
  cores <- droplevels(subset(Ooc, Ooc >= 30))
  length(cores$Ooc)
  Ooccore <- as.data.frame(rownames(cores))
  colnames(Ooccore) <- "gene"
  Ooccore$stage <- rep(i, length(Ooccore$gene))
  #dat <- data.frame(i, j, length(overlap), HVGnotpir, notHVGnotpir, notHVGyespir, test$p.value)
  #colnames(dat) <- c("Clust", "Core", "overlap", "clustnotcore", "neither", "notclustyescore", "pval")
  datalist[[i]] <- Ooccore
  
}
big_data = do.call(rbind, datalist)


datalist2 = list()
for (i in unique(big_data$stage)) {
  getit <- datalist[[i]]
  genes <- getit$gene
  allbut <- droplevels(subset(big_data, stage != i))
  allbutgenes <- allbut$gene
  unique <- as.data.frame(setdiff(genes, allbutgenes))
  colnames(unique) <- "unique"
  unique$stage <- rep(i, length(unique$unique))
  datalist2[[i]] <- unique
}

coreunique <- do.call(rbind, datalist2)
#write.csv(big_data, "CorebyStage_20190111.csv")
#write.csv(coreunique, "CoreUniquebyStage_20190111.csv")




clusters <- read.csv("Table_S2_ClusterAssignments_20180814.csv", header=TRUE)



datalist = list()

for (i in unique(clusters$Cluster_Name)) {
  clust <- clusters[clusters$Cluster_Name==i, ]
  for (j in colnames(coreunique$stage)) {
    stagecore <- coreunique[[j]]
    overlap <- intersect(clust$gene_id, stagecore$unique)
    HVGnotpir <- length(clust$gene_id) - length(overlap)
    notHVGnotpir <- length(clusters$gene_id) - length(stagecore$unique) - length(clust$gene_id)
    notHVGyespir <- length(clust$gene_id) - length(overlap)
    
    B <- matrix(c(HVGnotpir, notHVGnotpir, length(overlap), notHVGyespir), nrow=2, ncol=2)
    
    test <- chisq.test(B)
    dat <- data.frame(i, j, length(overlap), HVGnotpir, notHVGnotpir, notHVGyespir, test$p.value)
    colnames(dat) <- c("Clust", "Core", "overlap", "clustnotcore", "neither", "notclustyescore", "pval")
    datalist[[j]] <- dat
    
  }
}

datalist = list()

for (i in unique(clusters$Cluster_Name)) {
  clust <- clusters[clusters$Cluster_Name==i, ]
  for (j in unique(coreunique$stage)) {
    stagecore <- coreunique[coreunique$stage==j,]
    overlap <- intersect(clust$gene_id, stagecore$unique)
    HVGnotpir <- length(clust$gene_id) - length(overlap)
    notHVGnotpir <- length(clusters$gene_id) - length(stagecore$unique) - length(clust$gene_id)
    notHVGyespir <- length(clust$gene_id) - length(overlap)
    
    B <- matrix(c(HVGnotpir, notHVGnotpir, length(overlap), notHVGyespir), nrow=2, ncol=2)
    
    test <- fisher.test(B, alternative = "less")
    dat <- data.frame(i, j, length(overlap), HVGnotpir, notHVGnotpir, notHVGyespir, test$p.value)
    colnames(dat) <- c("Clust", "Core", "overlap", "clustnotcore", "neither", "notclustyescore", "pval")
    dat$adjpval <- p.adjust(dat$pval, method = "bonferroni", n = length(dat$pval))
    datalist[[j]] <- dat
    
  }
  big_data = do.call(rbind, datalist)
  #write.csv(big_data, file=paste0("bigdatachi_", i, ".csv", sep=""))
  write.table(big_data, file="bd2.csv", append=TRUE, sep=",")
}
big_data = do.call(rbind, datalist)



