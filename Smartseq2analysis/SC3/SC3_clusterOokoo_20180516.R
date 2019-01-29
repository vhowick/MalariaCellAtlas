setwd("/Users/vh3/Documents/MCA/ANALYSIS_2")
library(scater)
library(SC3)
##Read in data from QCed expression matrices, make SCE object and subset life stages
molecules <- read.table("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_counts.csv", header = TRUE, sep = ",", row.names=1, stringsAsFactors = TRUE)
anno <- read.delim("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_pheno.csv", header = TRUE, sep = ",")


cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "ooSpz" = "lightskyblue", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")


mca <- SingleCellExperiment(assays = list(
  counts = as.matrix(molecules),
  logcounts = log2(as.matrix(molecules) + 1)
), colData = anno)

mca2 <- mca[, (colData(mca)$ShortenedLifeStage2 == "ook") | 
                         (colData(mca)$ShortenedLifeStage2 == "ookoo") |
                         colData(mca)$ShortenedLifeStage2 == "oocyst"]


mc.qc.ookootmm <- scater::normaliseExprs(mca2, method = "TMM")
plotPCA(mc.qc.ookootmm, colour_by="PBANKA_0412900", shape_by="ShortenedLifeStage2")

rowData(mc.qc.ookootmm)$feature_symbol <- rownames(mc.qc.ookootmm)
#Run SC3, note that this will change cluster number name with each run
#pollen <- sc3(mc.qc.ookootmm, ks = 2:7, biology = TRUE)

PCA <- plotPCA(pollen, colour_by="sc3_5_clusters", shape_by="ShortenedLifeStage2")

dat <- PCA$data

ggplot(dat, aes(PC1, PC2)) + geom_point(aes(colour=factor(colour_by), shape = factor(shape_by))) + scale_color_brewer(palette="Set1") + theme_classic()

dat$sample_id <- row.names(dat)
colData(mc.qc.ookootmm) <- merge(colData(mc.qc.ookootmm), dat, by="sample_id", all.x=TRUE, all.y=FALSE)


nam <- as.list(c("colour_by", "sample_id"))
idx <- match(nam, names(dat))
dat <- dat[, idx]

dat$ookoo <- rep("Ookinete", length(dat$colour_by))
dat[which(dat$colour_by==1), ]$ookoo <- "Oocyst"
dat[which(dat$colour_by==2), ]$ookoo <- "Oocyst"
colnames(dat) <- c("SC3_5_clust", "sample_id", "ookoo")
write.csv(dat, file="OokooInfo_5clusts_20180516.csv")

plotPCA(pollen, colour_by="sc3_5_clusters", shape_by="ShortenedLifeStage2", ncomponents=3)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
