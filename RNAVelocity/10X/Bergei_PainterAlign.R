# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

## 10X
dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

up_mat <- velo_10X$deltaE
up_mat[up_mat < 0] <- 0;

# Read in published

paint_trans <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T, stringsAsFactors=FALSE)
paint_stab <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_stabilitzation_vs_time_FromSuppl1.csv", sep=",", header=T, stringsAsFactors=FALSE)
paint_stab[,51] <- as.numeric(paint_stab[,51])
identical(paint_stab[,1], paint_trans[,1])
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% paint_trans$Gene.ID,]


# Match genes
require("scater")
obj <- SingleCellExperiment(assays=list(emat=dat_10X$emat, smat=dat_10X$smat), colData=dat_10X$ann)
obj <- obj[rownames(obj) %in% ortho$Gene,]

paint_trans <- paint_trans[paint_trans$Gene.ID %in% ortho[,4],]
paint_trans$pbanka <- ortho[match(paint_trans$Gene.ID, ortho[,4]), 1]
paint_trans <- paint_trans[match(rownames(obj), paint_trans$pbanka),]
paint_stab <- paint_stab[paint_stab$Gene.ID %in% ortho[,4],]
paint_stab$pbanka <- ortho[match(paint_stab$Gene.ID, ortho[,4]), 1]
paint_stab <- paint_stab[match(rownames(obj), paint_stab$pbanka),]

paint_trans <- paint_trans[ !is.na(paint_trans[,51]), ]
paint_stab <- paint_stab[ !is.na(paint_stab[,51]), ]

emat <- assays(obj)[["emat"]][match(paint_stab$pbanka, rownames(assays(obj)[["emat"]])), ]
smat <- assays(obj)[["smat"]][match(paint_trans$pbanka, rownames(assays(obj)[["smat"]])), ]
emat <- t(apply(log2(emat+1), 1, function(a) {a<- a-min(a); a <- a/max(a);}))
smat <- t(apply(log2(smat+1), 1, function(a) {a<- a-min(a); a <- a/max(a);}))

## alternative
up_mat <- up_mat[match(paint_trans$pbanka, rownames(up_mat)),]
up_mat[is.na(up_mat)] <- 0



require("CellTypeProfiles")
require("gplots")

## Fancy plotting

source("~/R-Scripts/Blank_plot.R")

# bin cells by pseudotime
n_pbins=20
p_breaks <- seq(from=min(obj$pseudotime), to=max(obj$pseudotime), length=n_pbins+1)
p_bins <- cut(obj$pseudotime, breaks=n_pbins)
p_axis_ats <- seq(from=0, to=1, length=n_pbins)
p_axis_labs <- round( (p_breaks[1:n_pbins]+p_breaks[2:(n_pbins+1)])/2, digits=1 )

# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

fixed_row_mean_aggregate <-function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) {
			if (length(a) > 1) {
			rowMeans(MAT[, a])
			} else if (length(a) == 1) {
			MAT[,a]
			} else {
			rep(0, times=nrow(MAT))
			}
		})
    return(result)
}



# bin Painter
n_gbins=36
paint_breaks <- seq(from=0, to=48, length=n_gbins+1)
paint_breaks[1] <- paint_breaks[1]-0.00001
paint_breaks[n_gbins+1] <- paint_breaks[n_gbins+1]+0.00001

# bin genes by pf peak time
paint_stab_binned <- cut(paint_stab[,51], paint_breaks)
paint_trans_binned <- cut(paint_trans[,51], paint_breaks)

heat_data_s <- fixed_row_mean_aggregate(t(smat), paint_trans_binned);
heat_data_s[is.na(heat_data_s)] <- 0
heat_data_s <- fixed_row_mean_aggregate(t(heat_data_s), p_bins);
heat_data_e <- fixed_row_mean_aggregate(t(emat), paint_stab_binned);
heat_data_e[is.na(heat_data_e)] <- 0
heat_data_e <- fixed_row_mean_aggregate(t(heat_data_e), p_bins);
heat_data_e[is.na(heat_data_e)] <- 0

graphics::layout(mat=rbind(c(1,2), c(1,2)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))

#plot(obj$pseudotime, Matrix::colSums(smat), 
#	xlab="pseudotime", ylab="unspliced", 
#	xlim=c(min(obj$pseudotime)+0.22, max(obj$pseudotime)-0.22), xaxt="n", col="black")
#abline(v=p_breaks, lty=3)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=paint_breaks, z=t(heat_data_s), ylab="peak time (h)", xlab="pseudotime", col=up_col)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=paint_breaks, z=t(heat_data_e), ylab="peak time (h)", xlab="pseudotime", col=down_col)


