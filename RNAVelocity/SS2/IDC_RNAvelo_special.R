### IDC
require("scater")
require("velocyto.R")
qc_cells2 <- read.delim("allIDCgoodcells.csv", sep=",")
ann2 <- readRDS("allIDC.rds")
#ann2 <- toSingleCellExperiment(ann2)
cell_cols <- ann2$col
names(cell_cols) <- rownames(colData(ann2))
ref_cells <- names(cell_cols)
prefix="IDC"
anti <- readRDS("AntiSense/IDC_smoother_Antisense.rds")

stuff <- readRDS(paste(prefix, "_RNAVelo_Counts.rds", sep=""))
new_cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")
new_cols2 <- c("6" = "#78C679", "2"="#D1EC9F", "0"="#FEB24C", "1"="#F4CF63", "3" = "#FEEEAA", "4"="#85B1D3", "7"= "#9ecae1", "5" = "#C9E8F1", "Male"="#B7B7D8", "Female"="#9C96C6", "unassigned" = "black")

mapping <- read.table("~/Collaborations/MCA/SS2mappedto10xpbfortallulah.csv", header=T, sep=",")
cell_cols <- new_cols[ann_SS$ShortenedLifeStage2]
cell_cols <- new_cols2[as.character(mapping$stage_pred)]

# Exclude AntiSense
anti_c <- colnames(anti$exonic_expr)

velo_g <- rownames(stuff$smat)
anti_g <- rownames(anti$isAnti_by_cell)
is.anti <- anti$isAnti_by_cell[,match(colnames(stuff$emat),anti_c)]
is.anti <- is.anti[match(velo_g, anti_g), ]
is.anti[is.na(is.anti)] <- FALSE
stuff$smat[is.anti] <- 0
stuff$emat[is.anti] <- 0
stuff$nmat[is.anti] <- 0
###

rvel.qf <- readRDS(paste(prefix, "_RNAVelo_Estimates_noAnti.rds", sep=""))

#set.seed(2819)
#rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells, zero.offset=TRUE, smat=stuff$nmat)


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



pdf("IDC_customCoords2_RNAvelo_noAnti.pdf", width=4, height=4)
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
		mean(projected_pca[thispts,2]), len=0.07, lwd=1.4 )

}
}

dev.off()
