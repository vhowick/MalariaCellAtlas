### OOKOO
require("scater")
require("velocyto.R")
qc_cells2 <- read.delim("OOKOOgoodcells.csv", sep=",")
ann2 <- readRDS("OOKOO.rds")
#ann2 <- toSingleCellExperiment(ann2)
cell_cols <- ann2$col
names(cell_cols) <- rownames(colData(ann2))
ref_cells <- names(cell_cols)
image_name="OOKOO_RNAvelo.png"
prefix="OOKOO"
anti <- readRDS("AntiSense/OOKOO_smoother_Antisense.rds")

stuff <- readRDS(paste(prefix, "_RNAVelo_Counts.rds", sep=""))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]

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

set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells, zero.offset=TRUE, smat=stuff$nmat)


pca.out <- readRDS(paste(prefix, "_RNAVelo_PCAobj_noAnti.rds", sep=""))

emb <- pca.out$epc@scores[,1:2]

desired_emb <- read.delim("~/Collaborations/MCA/ookooPCAcoordforta6.csv", sep=",") # this is just from scater plotPCA
# redo plotPCA manually to get loadings/rotation
ntop=500;
exprs_values="logcounts";
method=c("prcomp","irlba");
ncomponents=2;
scores <- rowVars(assays(ann2)[[exprs_values]])
thresh <- sort(scores, decreasing=TRUE)[500]
exprs_to_plot <- assays(ann2)[[exprs_values]][scores >= thresh,]
pca <- prcomp(t(exprs_to_plot), rank. = ncomponents)
pcs <- pca$x
rot <- pca$rotation


tot_exp <- stuff$smat+stuff$emat
pca_current <- as.matrix(rvel.qf$current)[rownames(rvel.qf$current) %in% rownames(exprs_to_plot),]
pca_projected <- as.matrix(rvel.qf$projected)[rownames(rvel.qf$projected) %in% rownames(exprs_to_plot),]

change <- as.matrix(rvel.qf$deltaE)[rownames(rvel.qf$deltaE) %in% rownames(exprs_to_plot),]
change <- change[match(rownames(exprs_to_plot), rownames(change)),]
change[is.na(change)] <- 0
change <- change*5

exprs_projected <- exprs_to_plot+change

projected_pca <- t(exprs_projected) %*% rot;

current_pca <- t(exprs_to_plot) %*% rot

plot(current_pca, pch=16, col=cell_cols)
arrows(current_pca[,1], current_pca[,2], projected_pca[,1], projected_pca[,2], len=0.1, lty=3)
 

gridn=30
gridxes <- seq(from=min(projected_pca[,1], current_pca[,1]), to=max(projected_pca[,1], current_pca[,1]), length=gridn)
gridyes <- seq(from=min(projected_pca[,2], current_pca[,2]), to=max(projected_pca[,2], current_pca[,2]), length=gridn)

gridpts <- cbind(rep(gridxes, times=gridn), rep(gridyes, each=gridn))
binned <- cbind(cut(current_pca[,1], gridxes), cut(current_pca[,2], gridyes))

pdf("OOKOO_customCoords2_RNAvelo_noAnti.pdf", width=4, height=4)
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
