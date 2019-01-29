
setwd("/Users/vh3/Documents/MCA/ANALYSIS_3")
library(ggplot2)
library(SC3)
library(scmap)
library(scater)

#mca <- readRDS("mca.qc_20180625.rds")

##Read in data from QCed expression matrices, make SCE object and subset life stages
molecules <- read.table("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_counts.csv", header = TRUE, sep = ",", row.names=1, stringsAsFactors = TRUE)
anno <- read.delim("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_pheno.csv", header = TRUE, sep = ",")


cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "ooSpz" = "lightskyblue", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")


mca <- SingleCellExperiment(assays = list(
  counts = as.matrix(molecules),
  logcounts = log2(as.matrix(molecules) + 1)
), colData = anno)

mca.blood <- mca[, (colData(mca)$ShortenedLifeStage2 == "Trophozoite") | (colData(mca)$ShortenedLifeStage2 == "Schizont") | (colData(mca)$ShortenedLifeStage2 == "Female") | (colData(mca)$ShortenedLifeStage2 == "Male") | (colData(mca)$ShortenedLifeStage2 == "Ring")]

##Normalize with TMM, and viz data
mca.blood <- scater::normaliseExprs(mca.blood, method = "TMM")

pca <- plotPCA(mca.blood, colour_by="ShortenedLifeStage2", ncomp=3)

pcatab <- pca$data

cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "ooSpz" = "lightskyblue", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")
ggplot(pcatab, aes(PC1, PC2)) + geom_point(aes(colour=factor(colour_by))) + scale_color_manual(values = cols) + theme_classic()

#Run SC3, save RDS with this as cluster numbers change if you rerun

mca.blood2 <- sc3(mca.blood, ks = 5:9, biology = TRUE)
#saveRDS(mca.blood2, "SS2bloodstageswithSC3_20181125.rds")
pca <- plotPCA(mca.blood2, colour_by="sc3_8_clusters", shape_by ="ShortenedLifeStage2", size_by="PBANKA_0915000")
pca <- plotPCA(mca.blood2, colour_by="sc3_8_clusters", shape_by="ShortenedLifeStage2")

dat <- pca$data

ggplot(dat, aes(PC1, PC2)) + geom_point(aes(colour=factor(colour_by), shape=factor(shape_by))) + scale_color_viridis(discrete = TRUE, option = "A") +  theme_classic()
ggplot(dat, aes(PC1, PC2)) + geom_point(aes(colour=factor(colour_by), shape=factor(shape_by))) +
scale_color_colormap("SC3 cluster",
                    colormap = colormaps$rainbow, reverse = T, discrete = T) + 
  theme_classic() +
  labs(x="PC1: 30% variance", y="PC2: 16% variance") +
  theme(legend.position="none", axis.title=element_text(size=8), legend.text = element_text(size = 8), legend.title = element_text(size = 8), axis.text = element_text(size=8), axis.text.x = element_blank(), axis.text.y = element_blank())

ggplot(dat, aes(PC1, PC2)) + geom_point(aes(colour=factor(colour_by), shape=factor(shape_by))) +
  scale_color_colormap("SC3 cluster",
                       colormap = colormaps$rainbow, reverse = T, discrete = T) + 
  theme_classic() +
  labs(x="PC1: 30% variance", y="PC2: 16% variance") +
  theme(axis.title=element_text(size=10), legend.text = element_text(size = 8), legend.title = element_text(size = 10), axis.text = element_text(size=10), axis.text.x = element_blank(), axis.text.y = element_blank())

##Make scmap-cell index with these data
library(scmap)
sce <- selectFeatures(mca.blood2, suppress_plot = FALSE, n_features = 500)
table(rowData(sce)$scmap_features)

set.seed(1)
sce <- indexCell(sce)
names(metadata(sce)$scmap_cell_index)

saveRDS(sce, file = "SS2bloodstageIndex_20181126.rds")
