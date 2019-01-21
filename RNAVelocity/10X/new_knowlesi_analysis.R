# External Data
#pca1 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pb_tmm_PCA.csv", ",", header=TRUE)
#pca2 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pf_tmm_PCA.csv", ",", header=TRUE)
#pca2$cell_id <- sub(".1", "", pca2$cell_id)
#pca3 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pk_tmm_PCA.csv", ",", header=TRUE)

#labels <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/stage_pred.csv", ",", header=TRUE)
#labels <- labels[grep("^pf", labels[,1]),]
#labels[,1] <- sub("^pf_", "", labels[,1])

#labels <- labels[match(pca2$cell_id, labels[,1]),]
#labels <- labels[!is.na(labels[,1]),]
#labels$cell_id <- sub(".1", "", labels[,1])

#pca2 <- pca2[match(labels[,1], pca2$cell_id),]
ann <- read.table("~/Collaborations/MCA/For10XVelocity/pk150scmapclusts2method_20190110.csv", sep=",", header=TRUE)

# Read in Data

emat <- readRDS("new_pk_out_S_matrix.rds")
smat <- readRDS("new_pk_out_U_matrix.rds")

require("Matrix")

cell_names <- read.delim("new_pk_out_ca_matrix.csv", ",", header=TRUE)
cell_names <- sub("possorted_genome_bam_6EQP0:", "", cell_names$CellID)
cell_names <- sub("x$", "", cell_names)

gene_ann <- read.delim("new_pk_out_ra_matrix.csv", ",", header=TRUE)
gene_names <- gene_ann$Accession

rownames(emat) <- gene_names
rownames(smat) <- gene_names
colnames(emat) <- cell_names
colnames(smat) <- cell_names


# Consistent cells
emat <- emat[,match(ann[,"sample_id"], colnames(emat))]
smat <- smat[,match(ann[,"sample_id"], colnames(smat))]

keep_g <- Matrix::rowSums(emat > 0) > 20 & Matrix::rowSums(smat > 0) > 10

emat <- emat[keep_g,] 
smat <- smat[keep_g,] 

# Nicer annoation info
source("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/Colour_Scheme.R")
rownames(ann) <- as.character(ann[,"sample_id"])
ann$pc1 <- ann$PC1mean
ann$pc2 <- ann$PC2mean
ann$bulk <- ann$pbbulk
ann$cell_id <- ann$CellName

colors <- colors[order(as.numeric(names(colors)))]
names <- names[order(as.numeric(names(names)))]
ann$cell_color <- colors[as.numeric(as.character(ann$stage_pred))+1]
ann$cell_color[is.na(ann$cell_color)] <- colors["unassigned"]
ann$cell_name <- names[as.numeric(as.character(ann$stage_pred))+1]
ann$cell_name[is.na(ann$cell_name)] <- "unassigned"

# Save Nice input
dat <- list(emat=emat, smat=smat, ann=ann)
saveRDS(dat, file="new_pk_10X_velodata.rds")

# Actual Analysis

dat_local <- dat
require("velocyto.R")

### Fitting ###
set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)

saveRDS(rvel.qf, file="new_pk_10X_veloEst.rds")

poor_fits <- apply(rvel.qf$deltaE, 1, function(x) {x <- x[x!=0]; a <- table(sign(x)); max(a/sum(a))})
poor_fit_genes <- names(poor_fits)[poor_fits > 0.90]

dat_local$emat <- dat_local$emat[! rownames(dat_local$emat) %in% poor_fit_genes,]
dat_local$smat <- dat_local$smat[! rownames(dat_local$smat) %in% poor_fit_genes,]


set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)


saveRDS(rvel.qf, file="new_pk_10X_veloEst.rds")

### PCA ###
my_colors <- dat_local$ann$cell_color
names(my_colors) <- rownames(dat_local$ann)
png("new_knowlesi_10X_veloPCA.png", width=7, height=7, units="in", res=300)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1, show.grid.flow=TRUE, grid.n=25, cell.colors=my_colors)
dev.off()
saveRDS(pca.out, file="new_pk_10X_veloPCA.rds")


