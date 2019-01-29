setwd("/Users/vh3/Documents/MCA/ANALYSIS_3")
library(ggthemes)
library(ggbeeswarm)
library(TSCAN)
library(destiny)
library(reshape2)
library(plyr)
library(ggplot2)
library(devtools)
library(monocle)
library(M3Drop)
library(scater)
library(SLICER)
library("lle")

molecules <- read.table("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_counts.csv", header = TRUE, sep = ",", row.names=1, stringsAsFactors = TRUE)
anno <- read.delim("/Users/vh3/Documents/MCA/forgithub/Expression_Matrices/Smartseq2/SS2_pheno.csv", header = TRUE, sep = ",")


cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "ooSpz" = "lightskyblue", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")


mca <- SingleCellExperiment(assays = list(
  counts = as.matrix(molecules),
  logcounts = log2(as.matrix(molecules) + 1)
), colData = anno)

mca.qc.idc <- mca[, (colData(mca)$ShortenedLifeStage2 == "Merozoite") |(colData(mca)$ShortenedLifeStage2 == "Ring") |(colData(mca)$ShortenedLifeStage2 == "Trophozoite") | (colData(mca)$ShortenedLifeStage2 == "Schizont")]



mca.qc.idc <- scater::normaliseExprs(mca.qc.idc, method = "TMM")


mca.qc.idc$ShortenedLifeStage2 <- factor(
  mca.qc.idc$ShortenedLifeStage2,
  levels = c("Merozoite", "Ring", "Trophozoite", "Schizont")
)


cellLabels <- mca.qc.idc$ShortenedLifeStage2
idc <- logcounts(mca.qc.idc)
colnames(idc) <- cellLabels

idc <- logcounts(mca.qc.idc)
colnames(idc) <- mca.qc.idc$sample_id
#ookoo <- ookoo[!duplicated(ookoo), ]
tidc <- t(idc)
tidc <- tidc[!duplicated(tidc), ]
cellLabels <- mca.qc.idc$ShortenedLifeStage2
slicer_genes <- select_genes(tidc)


k <- select_k(t(idc[slicer_genes,]), kmin = 30, kmax=60)

slicer_traj_lle <- lle(t(idc[slicer_genes, ]), m = 2, k)$Y


reducedDim(mca.qc.idc, "LLE") <- slicer_traj_lle
prd <- plotReducedDim(mca.qc.idc, use_dimred = "LLE", colour_by = "ShortenedLifeStage2") +
  xlab("LLE component 1") + ylab("LLE component 2") +
  ggtitle("Locally linear embedding of cells from SLICER") + aes(geom_text(rownames(colData(mca.qc.idc))))
dat <- prd$data
ggplot(dat, aes(X1, X2)) + geom_point() + geom_text(aes(label=row.names(dat)), check_overlap = TRUE, size=3)

slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
start <- ends[2]


pseudotime_order_slicer <- cell_order(slicer_traj_graph, 1)
branches <- assign_branches(slicer_traj_graph, 1)

pseudotime_slicer2 <-
  data.frame(
    Timepoint = cellLabels,
    pseudotime = NA,
    State = branches,
    sample_id = mca.qc.idc$sample_id)

pseudotime_slicer2$pseudotime[pseudotime_order_slicer] <-
  1:length(pseudotime_order_slicer)
mca.qc.idc$pseudotime_slicer2 <- pseudotime_slicer2$pseudotime


ggplot(as.data.frame(colData(mca.qc.idc)), 
       aes(x = pseudotime_slicer2, 
           y = ShortenedLifeStage2, colour = ShortenedLifeStage2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("SLICER pseudotime (cell ordering)") +
  ylab("Timepoint") +
  theme_classic()


plotPCA(mca.qc.idc, colour_by = "pseudotime_slicer2", shape_by="ShortenedLifeStage2", exprs_values = "logcounts", ntop = 500)+ theme(axis.title=element_text(size=12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size=12))

idcslicer <- as.data.frame(cbind(mca.qc.idc$sample_id, mca.qc.idc$pseudotime_slicer2))
colnames(idcslicer) <- c("sample_id", "pseudotime_slicer2")
idcslicer$pseudotime_slicer2 <- as.numeric(idcslicer$pseudotime_slicer2)     

#write.csv(idcslicer, file="update_allIDC_pseudotime_Slicer_20180629.csv")

#saveRDS(mca.qc.idc, file="mca.qc.idc_pt_20180730.rds")
