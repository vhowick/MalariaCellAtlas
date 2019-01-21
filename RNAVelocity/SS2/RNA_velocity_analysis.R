# origin: /lustre/scratch118/malaria/team222/tallulah/SSII/newbatches/forTullulah

#args <- commandArgs(trailingOnly=TRUE)

dat1 <- readRDS("run23362_velocity_counts.rds")
dat2 <- readRDS("run24085_velocity_counts.rds")
dat3 <- readRDS("run25078_7_velocity_counts.rds")
dat4 <- readRDS("run25078_8_velocity_counts.rds")
dat5 <- readRDS("run25173_velocity_counts.rds")
dat6 <- readRDS("run20625_velocity_counts.rds")
dat7 <- readRDS("run24968_velocity_counts.rds")
dat8 <- readRDS("run25616_velocity_counts.rds")

all_dats <- list(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)

# Merge data
alldat <- list(emat=vector(), iomat=vector(), smat=vector(), base.df=vector(), exons=vector(), genes=vector(), expr.lstat=vector())


for (d in all_dats) {
	alldat$exons <- d$exons
        alldat$genes <- d$genes
        alldat$emat <- cbind(alldat$emat, d$emat)
        alldat$iomat <- cbind(alldat$iomat, d$iomat)
        alldat$smat <- cbind(alldat$smat, d$smat)

        if (is.null(dim(alldat$base.df))) {
                alldat$base.df <- d$base.df
        } else {
                new_min <- alldat$base.df[, 1] >  d$base.df[, 1]
                if (sum(new_min) > 0) {
                        alldat$base.df[new_min, 1] <- d$base.df[new_min, 1]
                }
                new_max <- alldat$base.df[, 2] <  d$base.df[, 2]
                if (sum(new_max) > 0) {
                        alldat$base.df[new_max, 2] <- d$base.df[new_max, 2]
                }
                new_max <- alldat$base.df[, 3] <  d$base.df[, 3]
                if (sum(new_max) > 0) {
                        alldat$base.df[new_max, 3] <- d$base.df[new_max, 3]
                }
        }
        if (is.null(dim(alldat$expr.lstat))) {
                alldat$expr.lstat <- d$expr.lstat
        } else {
                new_min <- alldat$expr.lstat[, 1] >  d$expr.lstat[, 1]
                if (sum(new_min) > 0) {
                        alldat$expr.lstat[new_min, 1] <- d$expr.lstat[new_min, 1]
                }
                new_max <- alldat$expr.lstat[, 2] <  d$expr.lstat[, 2]
                if (sum(new_max) > 0) {
                        alldat$expr.lstat[new_max, 2] <- d$expr.lstat[new_max, 2]
                }
                new_max <- alldat$expr.lstat[, 3] <  d$expr.lstat[, 3]
                if (sum(new_max) > 0) {
                        alldat$expr.lstat[new_max, 3] <- d$expr.lstat[new_max, 3]
                }
        }
}

alldat$emat <- cbind(dat1$emat, dat2$emat, dat3$emat, dat4$emat, dat5$emat, dat6$emat, dat7$emat, dat8$emat) 
alldat$smat <- cbind(dat1$smat, dat2$smat, dat3$smat, dat4$smat, dat5$smat, dat6$smat, dat7$smat, dat8$smat) 
alldat$iomat <- cbind(dat1$iomat, dat2$iomat, dat3$iomat, dat4$iomat, dat5$iomat, dat6$iomat, dat7$iomat, dat8$iomat) 

print("built whole mat")

do_QC <- function(dat, qc_cells) {
	velo_cells <- colnames(dat$emat)
	velo_cells <- sub(".bam", "", velo_cells)
	velo_cells <- strsplit(velo_cells, "/")
	velo_cells <- lapply(velo_cells, function(a){a[length(a)]})
	velo_cells <- as.vector(unlist(velo_cells))
	velo_cells <- sub("trim", "", velo_cells)
	good_cells <- paste(qc_cells[,2], qc_cells[,3], sep="_")
	matches <- match(good_cells, velo_cells)
	matches <- matches[!is.na(matches)]
	dat$emat <- dat$emat[, matches]
	dat$smat <- dat$smat[, matches]
	dat$iomat <- dat$iomat[, matches]
	# exonic read (spliced) expression matrix
	new_names <- qc_cells[good_cells %in% velo_cells,1]
	emat <- dat$emat;
	colnames(emat) <- new_names
	# intronic read (unspliced) expression matrix
	nmat <- dat$iomat;
	colnames(nmat) <- new_names
	# spanning read (intron+exon) expression matrix
	smat <- dat$smat;
	colnames(smat) <- new_names
	#clusters <- factor(colData(sce)[,"ShortenedLifeStage"])
	#names(clusters) <- colnames(sce)
	gene_filter <- Matrix::rowMeans(emat > 5) > 0.05 & Matrix::rowMeans(smat > 5) > 0.05 & !grepl(":", rownames(emat))
	emat <- emat[gene_filter, ]
	smat <- smat[gene_filter, ] + nmat[gene_filter, ]
	return(list(emat=emat, smat=smat, nmat=nmat));
}
if (args[1] =="spz") {
### SPZ
require("scater")
require("velocyto.R")
qc_cells1 <- read.delim("SPZgoodcells.csv", sep=",")
ann1 <- readRDS("SPZ.rds")
#ann1 <- toSingleCellExperiment(ann1)
stuff <- do_QC(alldat, qc_cells1)
cell_cols <- ann1$col
names(cell_cols) <- rownames(colData(ann1))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
image_name="SPZ_RNAvelo.png"
prefix="SPZ"

} else if (args[1] =="ookoo") {
### OOKOO
require("scater")
require("velocyto.R")
qc_cells2 <- read.delim("OOKOOgoodcells.csv", sep=",")
ann2 <- readRDS("OOKOO.rds")
#ann2 <- toSingleCellExperiment(ann2)
stuff <- do_QC(alldat, qc_cells2)
cell_cols <- ann2$col
names(cell_cols) <- rownames(colData(ann2))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
image_name="OOKOO_RNAvelo.png"
prefix="OOKOO"
} else if (args[1] =="EEF") {
### EEF
require("scater")
require("velocyto.R")
qc_cells3 <- read.delim("EEFgoodcells.csv", sep=",")
ann3 <- readRDS("EEF.rds")
#ann3 <- toSingleCellExperiment(ann3)
stuff <- do_QC(alldat, qc_cells3)
cell_cols <- colData(ann3)$col
names(cell_cols) <- rownames(colData(ann3))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
image_name="EEF_RNAvelo.png"
prefix="EEF"
} else if (args[1]=="idc") {
### IDC
require("scater")
require("velocyto.R")
qc_cells4 <- read.delim("allIDCgoodcells.csv", sep=",")
ann4 <- readRDS("allIDC.rds")
#ann4 <- toSingleCellExperiment(ann4)
stuff <- do_QC(alldat, qc_cells4)
cell_cols <- ann4$col
new_cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")
names(cell_cols) <- rownames(colData(ann4))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
image_name="IDC_RNAvelo.png"
prefix="IDC"
} else if (args[1]=="blood") {
### BLOOD
require("scater")
require("velocyto.R")
qc_cells5 <- read.delim("allBLOODgoodcells.csv", sep=",")
ann5 <- readRDS("allBLOOD.rds")
#ann5 <- toSingleCellExperiment(ann5)
stuff <- do_QC(alldat, qc_cells5)
#stuff <- filter_vals(stuff, threshold=10)
cell_cols <- ann5$col
names(cell_cols) <- rownames(colData(ann5))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)#[cell_cols=="purple" | cell_cols=="purple4"]
image_name="Blood_RNAvelo.png"
prefix="Blood"
} else {
### All
qc_cells1 <- read.delim("SPZgoodcells.csv", sep=",")
qc_cells2 <- read.delim("OOKOOgoodcells.csv", sep=",")
qc_cells3 <- read.delim("EEFgoodcells.csv", sep=",")
qc_cells4 <- read.delim("allIDCgoodcells.csv", sep=",")
qc_cells5 <- read.delim("allBLOODgoodcells.csv", sep=",")
ann1 <- readRDS("SPZ.rds")
ann2 <- readRDS("OOKOO.rds")
ann3 <- readRDS("EEF.rds")
ann4 <- readRDS("allIDC.rds")
ann5 <- readRDS("allBLOOD.rds")
all_cells <- rbind(qc_cells1, qc_cells2, qc_cells3, qc_cells4, qc_cells5)
stuff <- do_QC(alldat, all_cells)
#stuff <- filter_vals(stuff, 10);
cell_cols <- c(ann1$col, ann2$col, ann3$col, ann4$col, ann5$col)
cell_types <- c(ann1$ShortenedLifeStage, ann2$ShortenedLifeStage, ann3$ShortenedLifeStage, ann4$ShortenedLifeStage, ann5$ShortenedLifeStage)
names(cell_cols) <- c(colnames(ann1), colnames(ann2), colnames(ann3), colnames(ann4), colnames(ann5))
names(cell_types) <- c(colnames(ann1), colnames(ann2), colnames(ann3), colnames(ann4), colnames(ann5))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
cell_types <- cell_types[match(colnames(stuff$emat), names(cell_types))]
ref_cells <- names(cell_cols)
#ref_cells <- names(cell_cols)[cell_cols=="purple" | cell_cols=="steelblue" | cell_cols=="navy" | cell_cols=="purple4"]
image_name <- "All_RNAvelo.png"
prefix="All"
###
}

print("cell-type subset set-up")

filter_vals <- function(stuff, threshold) {
	e <- (stuff$emat < threshold & stuff$smat < threshold) & stuff$emat > 0
	s <- (stuff$emat < threshold & stuff$smat < threshold) & stuff$smat > 0

	stuff$emat[e] <- 0
	stuff$emat[s] <- 0
	return(stuff);
}


saveRDS(stuff, paste(prefix, "_RNAVelo_Counts.rds", sep=""))

set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells, zero.offset=TRUE, smat=stuff$nmat)
#pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), show.grid.flow=TRUE, grid.n=30, cell.colors=cell_cols)	

poor_fits <- apply(rvel.qf$deltaE, 1, function(x) {x <- x[x!=0]; a <- table(sign(x)); max(a/sum(a))})
poor_fit_genes <- names(poor_fits)[poor_fits > 0.90]

stuff$emat <- stuff$emat[! rownames(stuff$emat) %in% poor_fit_genes,]
stuff$smat <- stuff$smat[! rownames(stuff$smat) %in% poor_fit_genes,]
stuff$nmat <- stuff$nmat[! rownames(stuff$nmat) %in% poor_fit_genes,]

print("round 1")

set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells, zero.offset=TRUE, smat=stuff$nmat)
saveRDS(rvel.qf, paste(prefix, "_RNAVelo_Estimates.rds", sep=""))

print("round 2")
dim(rvel.qf$projected)

png(image_name, width=6, height=6, units="in", res=300)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), show.grid.flow=TRUE, grid.n=30, cell.colors=cell_cols)	
saveRDS(pca.out, paste(prefix, "_RNAVelo_PCAobj.rds", sep=""))
dev.off()

emb <- pca.out$epc@scores[,1:2]

scores <- Matrix::rowMeans(abs(rvel.qf$deltaE))
saveRDS(sort(scores), paste(prefix, "_RNAVelo_AllRankedGenes.rds", sep=""))
print("pca")
top_genes <- head(sort(scores, decreasing=T))

for (g in names(top_genes)) {

	png(paste(prefix, g, "rnavelo_fits.png", sep="_"), width=4*4, height=4, units="in", res=300)
	gene.relative.velocity.estimates(stuff$emat+0.00000001, stuff$smat+0.00000001, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, show.gene=g, cell.colors=cell_cols, cell.emb=emb, old.fit=rvel.qf, zero.offset=TRUE, smat=stuff$nmat)
	dev.off()

}

print("top_genes")

exit()
## Random Explanation Stuff

# Poor fit genes
require("scater")
require("velocyto.R")
qc_cells1 <- read.delim("SPZgoodcells.csv", sep=",")
ann1 <- readRDS("SPZ.rds")
#ann1 <- toSingleCellExperiment(ann1)
stuff <- do_QC(alldat, qc_cells1)
cell_cols <- ann1$col
names(cell_cols) <- rownames(colData(ann1))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
prefix="PoorFitExamples"
set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells, zero.offset=TRUE, smat=stuff$nmat)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), show.grid.flow=TRUE, grid.n=30, cell.colors=cell_cols)	
emb <- pca.out$epc@scores[,1:2]

poor_fits <- apply(rvel.qf$deltaE, 1, function(x) {x <- x[x!=0]; a <- table(sign(x)); max(a/sum(a))})
poor_fit_genes <- names(poor_fits)[poor_fits > 0.90]

g<-poor_fit_genes[3]
png(paste(prefix, g, "_velofit.png", sep="_"), width=4*4, height=4, units="in", res=300)
gene.relative.velocity.estimates(stuff$emat+0.00000001, stuff$smat+0.00000001, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, show.gene=g, cell.colors=cell_cols, cell.emb=emb, old.fit=rvel.qf, zero.offset=TRUE, smat=stuff$nmat)
dev.off()

# refrence cells
require("scater")
require("velocyto.R")
qc_cells2 <- read.delim("OOKOOgoodcells.csv", sep=",")
ann2 <- readRDS("OOKOO.rds")
#ann2 <- toSingleCellExperiment(ann2)
stuff <- do_QC(alldat, qc_cells2)
cell_cols <- ann2$col
names(cell_cols) <- rownames(colData(ann2))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
ref_cells <- names(cell_cols)
image_name="OOKOO_RNAvelo.png"
prefix="OOKOO"
ref_cells1 <- names(cell_cols)[cell_cols=="steelblue"]
ref_cells2 <- names(cell_cols)[cell_cols=="turquoise4"]
ref_cells3 <- names(cell_cols)[cell_cols=="mediumturquoise"]
prefix="RefCells"


png(paste(prefix, "_blue_velofit.png", sep="_"), width=4*4, height=4, units="in", res=300)
set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells1, zero.offset=TRUE, smat=stuff$nmat)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), show.grid.flow=TRUE, grid.n=30, cell.colors=cell_cols)
dev.off()

png(paste(prefix, "_mediumturquoise_velofit.png", sep="_"), width=4*4, height=4, units="in", res=300)
set.seed(2819)
rvel.qf <- gene.relative.velocity.estimates(stuff$emat, stuff$smat, deltaT=10, kCells = 5,fit.quantile = 0.075, min.nmat.emat.slope=0.1, min.nmat.emat.correlation=0.5, kGenes=1, steady.state.cells=ref_cells2, zero.offset=TRUE, smat=stuff$nmat)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), show.grid.flow=TRUE, grid.n=30, cell.colors=cell_cols)
dev.off()


### Summary Stats

### All
all_cells <- rbind(qc_cells1, qc_cells2, qc_cells3, qc_cells4, qc_cells5)
stuff <- do_QC(alldat, all_cells)
#stuff <- filter_vals(stuff, 10);
cell_cols <- c(ann1$col, ann2$col, ann3$col, ann4$col, ann5$col)
cell_types <- c(ann1$ShortenedLifeStage, ann2$ShortenedLifeStage, ann3$ShortenedLifeStage, ann4$ShortenedLifeStage, ann5$ShortenedLifeStage)
names(cell_cols) <- c(colnames(ann1), colnames(ann2), colnames(ann3), colnames(ann4), colnames(ann5))
names(cell_types) <- c(colnames(ann1), colnames(ann2), colnames(ann3), colnames(ann4), colnames(ann5))
cell_cols <- cell_cols[match(colnames(stuff$emat), names(cell_cols))]
cell_types <- cell_types[match(colnames(stuff$emat), names(cell_types))]
ref_cells <- names(cell_cols)
#ref_cells <- names(cell_cols)[cell_cols=="purple" | cell_cols=="steelblue" | cell_cols=="navy" | cell_cols=="purple4"]
image_name <- "All_RNAvelo.png"
prefix="All"
###

require("CellTypeProfiles")
e_by_type <- my_row_mean_aggregate(stuff$emat, cell_types)
s_by_type <- my_row_mean_aggregate(stuff$smat, cell_types)
i_by_type <- my_row_mean_aggregate(stuff$nmat, cell_types)

genes1 <- sort(rowMeans(e_by_type), decreasing=TRUE)[1:1000]


require("gplots")
genes <- genes1

png("Summary_ecounts_by_type.png", width=6, height=6, units="in", res=300)
toplot <- e_by_type[match(names(genes), rownames(e_by_type)), ]
thing <- heatmap.2(log2(toplot+1), col=colorRampPalette(c("blue", "white", "red"))(10), trace="none")
dev.off()

png("Summary_scounts_by_type.png", width=6, height=6, units="in", res=300)
toplot <- s_by_type[match(names(genes), rownames(s_by_type)), ]
heatmap.2(log2(toplot+1), col=colorRampPalette(c("blue", "white", "red"))(10), trace="none", breaks=thing$breaks, Rowv=thing$rowDendrogram, Colv=thing$colDendrogram)
dev.off()

png("Summary_icounts_by_type.png", width=6, height=6, units="in", res=300)
toplot <- i_by_type[match(names(genes), rownames(i_by_type)), ]
heatmap.2(log2(toplot+1), col=colorRampPalette(c("blue", "white", "red"))(10), trace="none", breaks=thing$breaks, Rowv=thing$rowDendrogram, Colv=thing$colDendrogram)
dev.off()

#ratio <- (log2(s_by_type+1.1))/(log2(e_by_type+1.1))
#toplot<-ratio[match(names(genes), rownames(ratio)), ]
#heatmap.2(log2(toplot+1), col=colorRampPalette(c("blue", "white", "red"))(10), trace="none", Rowv=thing$rowDendrogram, Colv=thing$colDendrogram)




# CellType Legend 
source("~/R-Scripts/Blank_plot.R")
t <- table(cell_types, cell_cols)
out <- which(t>0, arr.ind=TRUE)
out[,1] <- rownames(t)[as.numeric(out[,1])]
out[,2] <- colnames(t)[as.numeric(out[,2])]
png("Legend.png", width=2, height=4, units="in", res=300)
blank_plot()
legend("left", out[,1], fill=out[,2])
dev.off()















	alt_FS <- Matrix::rowSums(abs(rvel.qf$deltaE))
	all_FS <- names(alt_FS[alt_FS > quantile(alt_FS, probs=0.90)])
	all_FS <- unique(c(as.character(e_fs$Gene), as.character(s_fs$Gene)))

do_FS <- function(rvel.qf, features) {
	keep <- rownames(rvel.qf$conv.emat.norm) %in% features
	rvel.qf$conv.emat.norm <- rvel.qf$conv.emat.norm[keep,]
	rvel.qf$conv.nmat.norm <- rvel.qf$conv.nmat.norm[keep,]
	keep <- rownames(rvel.qf$projected) %in% features
	rvel.qf$projected <- rvel.qf$projected[keep,]
	rvel.qf$current <- rvel.qf$current[keep,]
	rvel.qf$ko <- rvel.qf$ko[rownames(rvel.qf$ko) %in% features,]
	rvel.qf$mval <- rvel.qf$mval[rownames(rvel.qf$mval) %in% features,]
	rvel.qf$gamma <- rvel.qf$gamma[names(rvel.qf$gamma) %in% features]
	rvel.qf$deltaE <- rvel.qf$deltaE[rownames(rvel.qf$deltaE) %in% features,]
	return(rvel.qf);
}

	forArthur <- list(Count_Mats = dat, Vel_Est=rvel.qf, PCA=stuff1)
	saveRDS(forArthur, file="SSII_Velocity_Out.rds")
