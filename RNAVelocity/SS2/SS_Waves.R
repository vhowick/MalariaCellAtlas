###### IDC
## Smartseq2
require("scater")
dat_SS <- readRDS("IDC_RNAVelo_Counts.rds")
velo_SS <- readRDS("IDC_RNAVelo_Estimates.rds")
ann_SS <- readRDS("allIDC.rds")
pca_SS <- readRDS("IDC_RNAVelo_PCAobj.rds")


new_cols <- c("bbSpz" = "navy", "EEF"="darkorange", "Merozoite"="lightpink", "oocyst"="steelblue", "ook" = "turquoise4", "Ring"="hotpink", "sgSpz"= "royalblue", "Schizont" = "violetred", "Male"="purple", "Female"="purple4", "ookoo" = "mediumturquoise", "Trophozoite"="violet")
new_cols2 <- c("6" = "#78C679", "2"="#D1EC9F", "0"="#FEB24C", "1"="#F4CF63", "3" = "#FEEEAA", "4"="#85B1D3", "7"= "#9ecae1", "5" = "#C9E8F1", "Male"="#B7B7D8", "Female"="#9C96C6", "unassigned" = "black")

tmp <- velo_SS$deltaE
tmp[tmp < 0] <- 0
cell_up <- Matrix::colSums(tmp)
up_s <- t(apply(tmp, 1, function(a) {a <- a-min(a); a/max(a)}))
cell_up_s <- Matrix::colSums(up_s)
dat_SS$nmat <- dat_SS$nmat[rownames(dat_SS$nmat) %in% rownames(dat_SS$smat),]
tmp <- dat_SS$smat+dat_SS$emat+dat_SS$nmat
ngene <- Matrix::colSums(tmp > 0);
nreads_per_g <- Matrix::colSums(tmp)/ngene

mapping <- read.table("~/Collaborations/MCA/SS2mappedto10xpbfortallulah.csv", header=T, sep=",")
pseudo10x <- readRDS("../10XOutput/bergei_10X_velodata_wPseudo.rds")
thing <- pseudo10x$ann[match(mapping[,1],pseudo10x$ann[,1]),]
mapping$pseudotime <- thing$pseudotime
#mapping$cell_col <- new_cols[ann_SS$ShortenedLifeStage2]
mapping$cell_col <- new_cols2[as.character(mapping$stage_pred)]

#png("Draft_IDC_SS_Waves.png", width=5*2/3, height=10*2/3, units="in", res=300)
pdf("Draft_IDC_SS_Waves.pdf", width=5*2/3, height=10*2/3)
layout(matrix(rep(1:4, each=2), ncol=2, byrow=T), widths=c(1,1), heights=c(1, 0.95, 0.95, 1.3))
par(mar=c(0,4,1,1))
plot(mapping$pseudotime, cell_up_s, col=mapping$cell_col, xlab="", ylab="Transcription (scaled)", xaxt="n")
thing <- smooth.spline(mapping$pseudotime, cell_up_s, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(mapping$pseudotime, cell_up, col=mapping$cell_col, xlab="", ylab="Transcription (unscaled)", xaxt="n")
thing <- smooth.spline(mapping$pseudotime, cell_up, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(mapping$pseudotime, ngene, col=mapping$cell_col, xlab="", ylab="nGene", xaxt="n")
thing <- smooth.spline(mapping$pseudotime, ngene, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(4,4,0,1))
plot(mapping$pseudotime, nreads_per_g, col=mapping$cell_col, xlab="", ylab="Reads/g", xaxt="n")
thing <- smooth.spline(mapping$pseudotime, nreads_per_g, spar=0.8)
lines(thing$x, thing$y, lwd=2)
mtext(side=1, line=2, "Pseudotime")
axis(1, at=0:6, labels=0:6)
dev.off()


### OOKOO
require("scater")
dat_SS <- readRDS("OOKOO_RNAVelo_Counts.rds")
velo_SS <- readRDS("OOKOO_RNAVelo_Estimates.rds")
ann_SS <- readRDS("OOKOO.rds")
pca_SS <- readRDS("OOKOO_RNAVelo_PCAobj.rds")
pseudo <- read.table("~/Collaborations/MCA/ookoo_pseudotime_Slicer_20180524.csv", header=T, sep=",")
pseudo <- pseudo[match(colnames(dat_SS$emat), pseudo$sample_id),3]

tmp <- velo_SS$deltaE
tmp[tmp < 0] <- 0
cell_up <- Matrix::colSums(tmp)
up_s <- t(apply(tmp, 1, function(a) {a <- a-min(a); a/max(a)}))
cell_up_s <- Matrix::colSums(up_s)
dat_SS$nmat <- dat_SS$nmat[rownames(dat_SS$nmat) %in% rownames(dat_SS$smat),]
tmp <- dat_SS$smat+dat_SS$emat+dat_SS$nmat
ngene <- Matrix::colSums(tmp > 0);
nreads_per_g <- Matrix::colSums(tmp)/ngene

cell_cols <- ann_SS$col

#png("Draft_OOKOO_SS_Waves.png", width=5*2/3, height=10*2/3, units="in", res=300)
pdf("Draft_OOKOO_SS_Waves.pdf", width=5*2/3, height=10*2/3)
layout(matrix(rep(1:4, each=2), ncol=2, byrow=T), widths=c(1,1), heights=c(1, 0.95, 0.95, 1.3))
par(mar=c(0,4,1,1))
plot(pseudo, cell_up_s, col=cell_cols, xlab="", ylab="Transcription (scaled)", xaxt="n")
thing <- smooth.spline(pseudo, cell_up_s, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(pseudo, cell_up, col=cell_cols, xlab="", ylab="Transcription (unscaled)", xaxt="n")
thing <- smooth.spline(pseudo, cell_up, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(pseudo, ngene, col=cell_cols, xlab="", ylab="nGene", xaxt="n")
thing <- smooth.spline(pseudo, ngene, spar=0.8)
lines(thing$x, thing$y, lwd=2)
par(mar=c(4,4,0,1))
plot(pseudo, nreads_per_g, col=cell_cols, xlab="", ylab="Reads/g", xaxt="n")
thing <- smooth.spline(pseudo, nreads_per_g, spar=0.8)
lines(thing$x, thing$y, lwd=2)
mtext(side=1, line=2, "Pseudotime")
axis(1, at=0:6, labels=0:6)
dev.off()

