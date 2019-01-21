# External Data
pca1 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pb_tmm_PCA.csv", ",", header=TRUE)
#pca2 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pf_tmm_PCA.csv", ",", header=TRUE)
#pca2$cell_id <- sub(".1", "", pca2$cell_id)
#pca3 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pk_tmm_PCA.csv", ",", header=TRUE)

labels <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/stage_pred.csv", ",", header=TRUE)
labels <- labels[-grep("^pf", labels[,1]),]
labels <- labels[-grep("^pk", labels[,1]),]

labels <- labels[match(pca1$cell_id, labels[,1]),]
labels <- labels[!is.na(labels[,1]),]

# Read in Data

emat <- readRDS("pbergei_out_S_Matrix.rds")
smat <- readRDS("pbergei_out_U_Matrix.rds")

require("Matrix")

cell_names <- read.delim("pbergei_out_ca_matrix.csv", ",", header=TRUE)
cell_names <- sub("possorted_genome_bam_UVJ4K:", "", cell_names$CellID)
cell_names <- sub("x$", "", cell_names)

gene_ann <- read.delim("pbergei_out_ra_matrix.csv", ",", header=TRUE)
gene_names <- gene_ann$Accession

rownames(emat) <- gene_names
rownames(smat) <- gene_names
colnames(emat) <- cell_names
colnames(smat) <- cell_names


# Consistent cells
emat <- emat[,match(labels[,1], colnames(emat))]
smat <- smat[,match(labels[,1], colnames(smat))]

keep_g <- rowSums(emat > 0) > 20 & rowSums(smat > 0) > 10

emat <- emat[keep_g,] 
smat <- smat[keep_g,] 

# Nicer annoation info
source("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/Colour_Scheme.R")
ann <- labels
rownames(ann) <- labels[,1]
pca1 <- pca1[match(labels[,1], pca1$cell_id),]
ann$pc1 <- pca1[,2]
ann$pc2 <- pca1[,3]
ann$bulk <- pca1[,4]
ann$cell_id <- pca1[,5]

colors <- colors[order(as.numeric(names(colors)))]
names <- names[order(as.numeric(names(names)))]
ann$cell_color <- colors[labels[,2]+1]
ann$cell_name <- names[labels[,2]+1]

# Save Nice input
dat <- list(emat=emat, smat=smat, ann=ann)
saveRDS(dat, file="bergei_10X_velodata.rds")

# Actual Analysis

dat_local <- dat
require("velocyto.R")

### Fitting ###
set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)

saveRDS(rvel.qf, file="bergei_10X_veloEst.rds")

poor_fits <- apply(rvel.qf$deltaE, 1, function(x) {x <- x[x!=0]; a <- table(sign(x)); max(a/sum(a))})
poor_fit_genes <- names(poor_fits)[poor_fits > 0.90]

dat_local$emat <- dat_local$emat[! rownames(dat_local$emat) %in% poor_fit_genes,]
dat_local$smat <- dat_local$smat[! rownames(dat_local$smat) %in% poor_fit_genes,]


set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)


saveRDS(rvel.qf, file="bergei_10X_veloEst.rds")

### PCA ###
my_colors <- dat_local$ann$cell_color
names(my_colors) <- rownames(dat_local$ann)
png("bergei_10X_veloPCA.png", width=7, height=7, units="in", res=300)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1, show.grid.flow=TRUE, grid.n=25, cell.colors=my_colors)
dev.off()
saveRDS(pca.out, file="bergei_10X_veloPCA.rds")


### Top Genes ###
emb <- pca.out$epc@scores[,1:2]

scores <- Matrix::rowMeans(abs(rvel.qf$deltaE))
saveRDS(sort(scores), "bergei_10x_velo_AllRankedGenes.rds")
print("pca")
top_genes_all <- head(sort(scores, decreasing=T), 50)
top_not_ribo <- c("PBANKA_1452300", "PBANKA_1424300", "PBANKA_1365500", "PBANKA_1309500", "PBANKA_1120700", "PBANKA_1116800", "PBANKA_1101100", "PBANKA_0819900", "PBANKA_0703900", "PBANKA_0701100", "PBANKA_0623300", "PBANKA_0604800", "PBANKA_0112600")
top_genes <- names(head(sort(scores, decreasing=T)))

for (g in top_not_ribo) {

        png(paste("bergei_10x", g, "rnavelo_fits.png", sep="_"), width=4*4, height=4, units="in", res=300)
        gene.relative.velocity.estimates(dat_local$emat+0.00000001, dat_local$smat+0.00000001, deltaT=10, kCells = 50,fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=1, show.gene=g, cell.colors=my_colors, cell.emb=emb, old.fit=rvel.qf, zero.offset=TRUE)
        dev.off()

}

