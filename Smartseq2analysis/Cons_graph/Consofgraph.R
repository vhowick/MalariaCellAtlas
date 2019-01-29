setwd("/Users/vh3/Documents/MCA/ANALYSIS_3/")
library(viridis)
library(ggplot2)
library(plotly)
library(scater)

coord <- read.csv("Clustname_Gene_pheno_manknn5_4_k20_20180814.csv", header=TRUE)
plasmogem <- read.csv("/Users/vh3/Documents/MCA/ANALYSIS_2/Plasmogem.csv", header=TRUE)
orthtable <- read.csv("/Users/vh3/Documents/MCA/10x/orth_table.csv", header=TRUE)
meg <- merge(coord, plasmogem, by.x="gene_id", by.y="current_version_ID", all.x=TRUE, all.y=FALSE)
orthmeg <- merge(coord, orthtable, by.x="gene_id", by.y="new_gene_name_pb", all.x=TRUE, all.y=FALSE)

orthmeg$pforth <- rep("yes", length(orthmeg$gene_id))
orthmeg[which(is.na(orthmeg$Plasmodium.falciparum.3D7.gene.stable.ID)), ]$pforth <- "no"

df <- as.data.frame(table(orthmeg$pforth, orthmeg$cluster))

yn <- c("yes"="grey88", "no"="blue")
ggplot(orthmeg, aes(X.x, Y)) + 
  geom_point(aes(colour=factor(pforth)), size = 0.3) + 
  theme_classic() +
  scale_color_manual(values = yn) +
  labs(x="Dimension 1", y="Dimension 2") +
  theme(legend.position= "none", axis.title=element_text(size=10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), axis.text = element_text(size=10), axis.text.x = element_blank(), axis.text.y = element_blank())

noorth <- droplevels(subset(orthmeg, pforth=="no"))
noorthgenes <- noorth$gene_id

datalist = list()

for (i in unique(orthmeg$Cluster_Name)) {
  clust <- orthmeg[orthmeg$Cluster_Name==i, ]
  
  overlap <- intersect(clust$gene_id, noorthgenes)
  HVGnotpir <- length(clust$gene_id) - length(overlap)
  notHVGnotpir <- 5156 - length(noorthgenes) - length(clust$gene_id)
  notHVGyespir <- length(noorthgenes) - length(overlap)
  
  B <- matrix(c(HVGnotpir, notHVGnotpir, length(overlap), notHVGyespir), nrow=2, ncol=2)
  
  test <- fisher.test(B, alternative = "less")
  dat <- data.frame(i, length(overlap)-1, HVGnotpir, notHVGnotpir, notHVGyespir, test$p.value)
  colnames(dat) <- c("Clust", "overlap", "inclustnotcat", "neither", "notclustyescat", "pval")
  datalist[[i]] <- dat
  
}
big_data = do.call(rbind, datalist)

big_data$adjpval <- p.adjust(big_data$pval, method = "fdr", n = length(big_data$pval))
big_data$bonpval <- p.adjust(big_data$pval, method = "bonferroni", n = length(big_data$pval))

# write.csv(big_data, file = "noorthFETenrichment_20190122.csv")


p <- ggplot(orthmeg, aes(x = as.factor(Cluster_Name), fill = pforth))
p  + geom_bar() + scale_fill_manual(values = yn) + 
  labs(x="Cluster", y="Gene count") +
  theme_classic() + 
  theme(axis.title=element_text(size=10), legend.text = element_text(size = 10), legend.title = element_blank(), axis.text = element_text(size=8))

)