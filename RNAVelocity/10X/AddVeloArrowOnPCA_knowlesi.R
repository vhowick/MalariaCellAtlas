### knowlesi 10X
require("scater")
#dat_10X <- readRDS("knowlesi_10X_velodata.rds")
#velo_10X <- readRDS("knowlesi_10X_veloEst.rds")
dat_10X <- readRDS("new_pk_10X_velodata.rds")
velo_10X <- readRDS("new_pk_10X_veloEst.rds")

tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)
remap_colors <- function(old_colors) {
        for (i in 1:nrow(tab_color_remap)) {
                old_colors[old_colors == tab_color_remap[i,1]] <- tab_color_remap[i,2]
        }
        return(old_colors);
}

count_mat <- readMM("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/raw_gene_bc_matrices/Pknowlesi_at10_20170115/matrix.mtx");
g <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/raw_gene_bc_matrices/Pknowlesi_at10_20170115/genes.tsv")
c <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/raw_gene_bc_matrices/Pknowlesi_at10_20170115/barcodes.tsv")
c <- sub("-1", "", c[,1])
g <- g[,1]

rownames(count_mat) <- g
colnames(count_mat) <- c

count_mat <- count_mat[,colnames(count_mat) %in% colnames(dat_10X$emat)]

mca_pb <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(count_mat))))

mca_pb_tmm <- scater::normaliseExprs(mca_pb, method = "TMM", return_norm_as_exprs = FALSE)
mca_pb_tmm <- mca_pb_tmm[,match(colnames(dat_10X$emat), colnames(mca_pb_tmm))]
plotPCA(mca_pb_tmm, exprs_values = "norm_exprs", ntop = 50)


# redo plotPCA manually to get loadings/rotation
require("matrixStats")
ntop=50;

exprs_values="norm_exprs";
method=c("prcomp","irlba");
ncomponents=2;
scores <- rowVars(as.matrix(assays(mca_pb_tmm)[[exprs_values]]))
thresh <- sort(scores, decreasing=TRUE)[ntop]
exprs_to_plot <- assays(mca_pb_tmm)[[exprs_values]][scores >= thresh,]
exprs_to_plot_unscaled <- exprs_to_plot
exprs_to_plot <- t(apply(exprs_to_plot, 1, scale))
pca <- prcomp(t(exprs_to_plot), rank. = ncomponents)
pcs <- pca$x
rot <- pca$rotation

plot(pcs[,1], -pcs[,2], pch=16, cex=0.7)


tot_exp <- dat_10X$smat+dat_10X$emat
pca_current <- as.matrix(velo_10X$current)[rownames(velo_10X$current) %in% rownames(exprs_to_plot),]
pca_projected <- as.matrix(velo_10X$projected)[rownames(velo_10X$projected) %in% rownames(exprs_to_plot),]

change <- as.matrix(velo_10X$deltaE)[rownames(velo_10X$deltaE) %in% rownames(exprs_to_plot),]
change <- change[match(rownames(exprs_to_plot), rownames(change)),]
change[is.na(change)] <- 0
change <- change*0.75
#change <- change*2.5*1.2

exprs_projected <- exprs_to_plot+change

projected_pca <- t(exprs_projected) %*% rot;

current_pca <- t(exprs_to_plot) %*% rot

current_pca[,2] <- -current_pca[,2]
projected_pca[,2] <- -projected_pca[,2]

plot(current_pca[,1], current_pca[,2], pch=16, col=dat_10X$ann$cell_color)
arrows(current_pca[,1], current_pca[,2], projected_pca[,1], projected_pca[,2], len=0.1, lty=3)

gridn=25
gridxes <- seq(from=min(projected_pca[,1], current_pca[,1]), to=max(projected_pca[,1], current_pca[,1]), length=gridn)
gridyes <- seq(from=min(projected_pca[,2], current_pca[,2]), to=max(projected_pca[,2], current_pca[,2]), length=gridn)

gridpts <- cbind(rep(gridxes, times=gridn), rep(gridyes, each=gridn))
binned <- cbind(cut(current_pca[,1], gridxes), cut(current_pca[,2], gridyes))
binned[is.na(binned[,1]),1] <- 1
binned[is.na(binned[,2]),2] <- 1

#png("Knowlesi_SuperPCAWithArrows.png", width=4, height=4, units="in", res=300)
pdf("Knowlesi_SuperPCAWithArrows.pdf", width=4, height=4)
par(mar=c(1,1,1,1))
plot(current_pca, pch=21, bg=remap_colors(dat_10X$ann$cell_color), col=remap_colors(dat_10X$ann$cell_color), cex=1, xaxt="n", yaxt="n", xlab="", ylab="")
box(lwd=2)
#points(rep(gridxes, times=gridn), rep(gridyes, each=gridn), pch=16, cex=0.3)

for(i in min(binned[,1]):max(binned[,1])) {
for(j in min(binned[,1]):max(binned[,1])) {
	thispts <- binned[,1] == i & binned[,2] == j
	if (sum(thispts) < 5) {next;}
	arrows( mean(current_pca[thispts,1]),
		mean(current_pca[thispts,2]),
		mean(projected_pca[thispts,1]),
		mean(projected_pca[thispts,2]), len=0.07, lwd=1.4 )

}
}

dev.off()
