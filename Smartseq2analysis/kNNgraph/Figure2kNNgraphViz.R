library(readxl)
library(ggplot2)
library(scater)
library(pheatmap)

sheets <- excel_sheets("/Users/vh3/Documents/MCA/MCA_File_S1_20190121.xlsx")

test <- as.data.frame(sheets["Gene Info"])
dats <- lapply(sheets, read_excel, path = "/Users/vh3/Documents/MCA/MCA_File_S1_20190121.xlsx")

bigtable <- as.data.frame(dats[[1]])
bigtable$Asex.Relative.Growth.Rate <- as.numeric(as.character(bigtable$Asex.Relative.Growth.Rate))

bigtabsub <- droplevels(subset(bigtable, Asex.Relative.Growth.Rate != "NA"))

ggplot(bigtabsub, aes(knngraph_X, knngraph_Y)) + 
  geom_point(aes(colour=Asex.Relative.Growth.Rate), size = 0.5) + scale_color_viridis(option = "C") +
  #geom_point(data=nodat, colour="grey", alpha=0.2, size = 0.8) +
  labs(x="Dimension 1", y="Dimension 2") +
  theme_classic() + 
  theme(axis.title=element_text(size=8), legend.text = element_text(size = 12), legend.title = element_blank(), axis.text = element_text(size=8), axis.text.x = element_blank(), axis.text.y = element_blank())


clustcol <- c("4" = "mediumseagreen", "5" ="orange", "19" ="gold2" , "1" ="limegreen", "15" ="palegreen", "12" ="cyan3", "8" ="cornflowerblue", "13" ="darkmagenta", "16" ="cadetblue", "7" ="dodgerblue", "20" ="salmon", "11" ="deepskyblue4", "3" ="purple", "14" ="thistle3", "18" ="lightseagreen", "17" ="blue", "9" ="pink", "10" ="red", "2" ="darkgrey", "6" ="darkslategray3")

bigtable$Cluster <- as.factor(bigtable$Cluster)
ggplot(bigtable, aes(knngraph_X, knngraph_Y)) + 
  geom_point(aes(colour=factor(Cluster)), size = 0.5) + 
  theme_classic() +
  scale_color_manual(values = clustcol) +
  labs(x="Dimension 1", y="Dimension 2") +
  theme(legend.position= "none", axis.title=element_text(size=8), legend.text = element_text(size = 8), legend.title = element_text(size = 8), axis.text = element_text(size=8), axis.text.x = element_blank(), axis.text.y = element_blank())

######MAKE HEATMAP FOR FIG 2
molecules <- read.table("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_counts.csv", header = TRUE, sep = ",", row.names=1, stringsAsFactors = TRUE)
anno <- read.delim("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_pheno.csv", header = TRUE, sep = ",")


cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "ooSpz" = "lightskyblue", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")

###Make SCE and normalize in groups of stages
mca.qc <- SingleCellExperiment(assays = list(
  counts = as.matrix(molecules),
  logcounts = log2(as.matrix(molecules) + 1)
), colData = anno)

rowData(mca.qc) <- bigtable


mca.qc.ookoo <- mca.qc[, (colData(mca.qc)$ShortenedLifeStage == "ook") | (colData(mca.qc)$ShortenedLifeStage == "ookoo") | (colData(mca.qc)$ShortenedLifeStage == "oocyst")]
mca.qc.eef <- mca.qc[, (colData(mca.qc)$ShortenedLifeStage == "EEF")]
mca.qc.spz <- mca.qc[, (colData(mca.qc)$ShortenedLifeStage == "bbSpz") | 
                       (colData(mca.qc)$ShortenedLifeStage == "sgSpz") |
                       (colData(mca.qc)$ShortenedLifeStage == "ooSpz")]
mca.qc.idc <- mca.qc[, (colData(mca.qc)$ShortenedLifeStage == "Merozoite") | 
                       (colData(mca.qc)$ShortenedLifeStage == "Ring") |
                       (colData(mca.qc)$ShortenedLifeStage2 == "Schizont") | 
                       (colData(mca.qc)$ShortenedLifeStage2 == "Trophozoite")]
mca.qc.sex <- mca.qc[, (colData(mca.qc)$ShortenedLifeStage2 == "Male") | (colData(mca.qc)$ShortenedLifeStage2 == "Female")]


mca.qc.ookoo.tmm <- scater::normaliseExprs(mca.qc.ookoo, method = "TMM")
#mca.qc.ooc.tmm <- scater::normaliseExprs(mca.qc.ooc, method = "TMM")
mca.qc.eef.tmm <- scater::normaliseExprs(mca.qc.eef, method = "TMM")
mca.qc.spz.tmm <- scater::normaliseExprs(mca.qc.spz, method = "TMM")
mca.qc.idc.tmm <- scater::normaliseExprs(mca.qc.idc, method = "TMM")
mca.qc.sex.tmm <- scater::normaliseExprs(mca.qc.sex, method = "TMM")

mca.qc.tmm <- cbind(mca.qc.ookoo.tmm, mca.qc.eef.tmm)
#mca.qc.tmm <- cbind(mca.qc.tmm, mca.qc.eef.tmm)
mca.qc.tmm <- cbind(mca.qc.tmm, mca.qc.spz.tmm)
mca.qc.tmm <- cbind(mca.qc.tmm, mca.qc.idc.tmm)
mca.qc.tmm <- cbind(mca.qc.tmm, mca.qc.sex.tmm)


ppt <- read.csv("allpptinfo_20180629.csv", header=TRUE) #read in psuedopseudotime to order cells

co <- merge(colData(mca.qc.tmm), ppt, by="sample_id", all.x=TRUE, all.y=FALSE)
subco <- cbind(co$sample_id, co$ppt)
subco <- as.data.frame(subco)
colnames(subco) <- c("sample_id", "ppt")

cellnames <- mca.qc.tmm$sample_id
subco <- subco[match(cellnames, subco$sample_id), ]

mca.qc.tmm$ppt <- as.numeric(as.character(subco$ppt))


ord <- mca.qc.tmm[, order(mca.qc.tmm$ppt)]
#ord <- ord[order(rowData(ord)$Cluster), ]

#get the colData and sort that based on dev
celldata <- colData(ord)
sortcelldata <- celldata[order(as.numeric(celldata$ppt)), ]


ordexp <- logcounts(ord)


#rename the matrix with sorted cellname, also isolate stage for annotation and make rownames the same as colnames of matrix (must be unique)
colnames(ordexp) <- sortcelldata$sample_id

stage <- as.data.frame(sortcelldata$ShortenedLifeStage2)
row.names(stage) <- colnames(ordexp)
colnames(stage) <- "stage"


geneclust <- as.data.frame(rowData(ord)$Cluster)
row.names(geneclust) <- rowData(ord)$feature_symbol
colnames(geneclust) <- "genecluster"



theclust <- cbind(ordexp, geneclust)

theclust <- as.data.frame(theclust)

clustmean <- aggregate(theclust[, 1:1787], by = list(as.factor(as.character(theclust$genecluster))), mean)


clustmean[1:5, 1:5]

subclustmean <- clustmean[, 2:1788]
rownames(subclustmean) <- clustmean$Group.1
subclustmean[1:5, 1:5]

pheatmap(subclustmean, cluster_cols=FALSE, cluster_rows=TRUE ,  annotation_col = stage,
         show_colnames = FALSE, show_rownames = TRUE, color=inferno(10))


ann_colors = list(stage = c(bbSpz = "navy", EEF="darkorange", Merozoite="lightpink", oocyst="steelblue", ook = "turquoise4", Ring="hotpink", sgSpz= "royalblue", Schizont = "violetred", Male="purple", Female="purple4", ookoo = "mediumturquoise", Trophozoite="violet"))

subclustmean2 <- subclustmean
rownames(subclustmean2) <- as.numeric(rownames(subclustmean2))
subclustmean2 <- subclustmean2[ order(as.numeric(rownames(subclustmean2))), ]
pheatmap(subclustmean2, cluster_cols=FALSE, cluster_rows=FALSE ,  annotation = stage,
         show_colnames = FALSE, show_rownames = TRUE, color=inferno(10), annotation_colors = ann_colors)





