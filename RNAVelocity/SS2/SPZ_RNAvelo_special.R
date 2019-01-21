### SPZ
require("scater")
require("velocyto.R")
qc_cells2 <- read.delim("SPZgoodcells.csv", sep=",")
ann2 <- readRDS("SPZ.rds")
#ann2 <- toSingleCellExperiment(ann2)
cell_cols <- ann2$col
names(cell_cols) <- rownames(colData(ann2))
ref_cells <- names(cell_cols)
prefix="SPZ"
anti <- readRDS("AntiSense/SPZ_smoother_Antisense.rds")

stuff <- readRDS(paste(prefix, "_RNAVelo_Counts.rds", sep=""))

cell_cols <- ann2$col


rvel.qf <- readRDS(paste(prefix, "_RNAVelo_Estimates_noAnti.rds", sep=""))

pca.out <- readRDS(paste(prefix, "_RNAVelo_PCAobj_noAnti.rds", sep=""))

emb <- pca.out$epc@scores[,1:2]

pcs <- pca.out$epc@scores
rot <- pca.out$epc@loadings

current_pca <- pca.out$epc@scores
projected_pca <- current_pca + pca.out$delta.pcs*0.75

gridn=20
gridxes <- seq(from=min(projected_pca[,1], current_pca[,1]), to=max(projected_pca[,1], current_pca[,1]), length=gridn)
gridyes <- seq(from=min(projected_pca[,2], current_pca[,2]), to=max(projected_pca[,2], current_pca[,2]), length=gridn)

gridpts <- cbind(rep(gridxes, times=gridn), rep(gridyes, each=gridn))
binned <- cbind(cut(current_pca[,1], gridxes), cut(current_pca[,2], gridyes))
binned[is.na(binned)] <- 1



pdf("SPZ_customCoords2_RNAvelo_noAnti.pdf", width=4, height=4)
par(mar=c(4,4,1,1))
plot(current_pca, pch=21, bg=cell_cols, col=cell_cols, cex=1)
box(lwd=2)
#points(rep(gridxes, times=gridn), rep(gridyes, each=gridn), pch=16, cex=0.3)

for(i in min(binned[,1]):max(binned[,1])) {
for(j in min(binned[,1]):max(binned[,1])) {
	thispts <- binned[,1] == i & binned[,2] == j
	if (sum(thispts) < 3) {next;}
	arrows( mean(current_pca[thispts,1]),
		mean(current_pca[thispts,2]),
		mean(projected_pca[thispts,1]),
		mean(projected_pca[thispts,2]), len=0.07, lwd=1.4, col="grey65" )

}
}

dev.off()
