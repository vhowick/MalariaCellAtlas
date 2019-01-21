# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

## 10X
dat_10X <- readRDS("bergei_10X_velodata.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

tmp<- velo_10X$deltaE
tmp[tmp < 0] <- 0
cell_up <- Matrix::colSums(tmp)
tmp<- velo_10X$deltaE
tmp[tmp > 0] <- 0
cell_down <- Matrix::colSums(tmp)
tmp<- velo_10X$deltaE
cell_both <- Matrix::colSums(tmp)


breaks_down <- seq(from=min(cell_down), to=0, length=7)
breaks_up <- seq(from=0, to=max(cell_up), length=7)

cell_down_col <- rev(down_col)[cut(cell_down, breaks=breaks_down)]
cell_up_col <- up_col[cut(cell_up, breaks=breaks_up)]


png("10X_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("10X_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
arrows(1,-1,-5,4.5, len=0) #y1
arrows(2,-1,6,6, len=0) #y2
arrows(3,-2,13.5,-2, len=0) #y3
arrows(2,-3,7.5,-8, len=0) #y4
arrows(-8,-2,0,-2, len=0) #y5
dev.off()

m1 <- (4.5 - -1)/(-5 - 1); b1 <- -1 - m1*1;
m2 <- (6 - -1)/(6 - 2); b2 <- -1 - m2*2;
m3 <- 0; b3 <- -2;
m4 <- (-8 - -3)/(7.5 - 2); b4 <- -3 - m4*2;
segment1 <- dat_10X$ann$pc2 > dat_10X$ann$pc1*m1+b1 & dat_10X$ann$pc2 > dat_10X$ann$pc1*m2+b2
segment2 <- dat_10X$ann$pc2 < dat_10X$ann$pc1*m2+b2 & dat_10X$ann$pc2 > dat_10X$ann$pc1*m3+b3
segment3 <- dat_10X$ann$pc2 < dat_10X$ann$pc1*m3+b3 & dat_10X$ann$pc2 > dat_10X$ann$pc1*m4+b4
segment4 <- dat_10X$ann$pc2 < dat_10X$ann$pc1*m4+b4 & dat_10X$ann$pc2 < -2
segments <- rep(5, times=nrow(dat_10X$ann))
segments[segment1] <- 1
segments[segment2] <- 2
segments[segment3] <- 3
segments[segment4] <- 4

dat_10X$ann$segment <- segments
saveRDS(dat_10X, "bergei_10X_velodata.rds")

png("10X_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=dat_10X$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")
dev.off()

## Smartseq2
require("scater")
dat_SS <- readRDS("../NewBatches/IDC_RNAVelo_Counts.rds")
velo_SS <- readRDS("../NewBatches/IDC_RNAVelo_Estimates.rds")
ann_SS <- readRDS("../NewBatches/allIDC.rds")
pca_SS <- readRDS("../NewBatches/IDC_RNAVelo_PCAobj.rds")

write.table(pca_SS$epc@scores, file="ForGinnySS2VelocityPCACellCoords.csv", sep=",")

cell_colors <- ann_SS$col
names(cell_colors) <- rownames(colData(ann_SS))
cell_colors <- cell_colors[match(colnames(dat_SS$emat), names(cell_colors))]

tmp <- velo_SS$deltaE
tmp[tmp < 0] <- 0
cell_up <- Matrix::colSums(tmp)
tmp <- velo_SS$deltaE
tmp[tmp > 0] <- 0
cell_down <- Matrix::colSums(tmp)

breaks_down <- seq(from=min(cell_down), to=0, length=7)
breaks_up <- seq(from=0, to=max(cell_up), length=7)

cell_down_col <- rev(down_col)[cut(cell_down, breaks=breaks_down)]
cell_up_col <- up_col[cut(cell_up, breaks=breaks_up)]

png("SS_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_colors, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("SS_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("SS_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()

## Match to Published ##

# Read in published

extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]

## 10X velocyto
dat_10X <- readRDS("bergei_10X_velodata.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

up<- velo_10X$deltaE
up[up < 0] <- 0
cell_up <- Matrix::colSums(up)

# Match genes
up <- up[rownames(up) %in% ortho$Gene,]
extern <- extern[extern$Gene.ID %in% ortho[,4],]
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern <- extern[match(rownames(up), extern$pbanka),]

require("CellTypeProfiles")
require("gplots")
M <- my_row_mean_aggregate(up, dat_10X$ann$segment)
M <- M[order(extern[,51]),]
M2 <- t(apply(M, 1, scale))

heatmap.2(M2, Rowv=FALSE, trace="none")


# Failed still 2d circle
#require("destiny")
#thing <- destiny::DiffusionMap(dat_10X$ann[,c("pc1", "pc2")])
#junk <- brewer.pal(8, "Greys")
#plot(dat_10X$ann[,"pc1"], dat_10X$ann[, "pc2"],  col=junk[cut(thing@eigenvectors[,1], breaks=8)], pch=16)

# Slingshot pseudotime # singluar matrix error
#require("slingshot")
#dat_mat <- dat_10X$emat+dat_10X$smat
##dat_mat <- dat_mat[match(rownames(M), rownames(dat_mat)),]
#dat_mat <- dat_mat[Matrix::rowMeans(dat_mat>0) > 0.0001,]
#sf <- Matrix::colSums(dat_mat)
#dat_norm <- t( t(dat_mat)/sf*median(sf) )
#dat_lognorm <- log2(dat_norm+1)

#c_keep <- sf > 500

#pseudo <- slingshot(t(as.matrix(dat_lognorm)[,c_keep]), dat_10X$ann$cell_color[c_keep], reducedDim=dat_10X$ann[c_keep,c("pc1", "pc2")], start.clus="mistyrose", end.clus="violetred4")

## Oscope - slow
#require("Oscope")
#dat_mat <- dat_10X$emat+dat_10X$smat
#dat_mat <- dat_mat[Matrix::rowMeans(dat_mat>0) > 0.01,]
#sf <- Matrix::colSums(dat_mat)
#dat_norm <- t( t(dat_mat)/sf*median(sf) )

##Oscope
#dat_input <- NormForSine(as.matrix(dat_norm))
#dat_input <- dat_input[!is.na(rowSums(dat_input)), ]
#SineRes <- OscopeSine(dat_input)
#KMRes <- OscopeKM(SineRes, maxK = 10)
#ToRM <- FlagCluster(SineRes, KMRes, dat_input)
#KMResUse <- KMRes[-ToRM$FlagID]
#ENIRes <- OscopeENI(KMRes = KMResUse, Data = dat_input, NCThre = 100)

# fit an ellipse
require("conicfit")
fita <- EllipseDirectFit(dat_10X$ann[,c("pc1", "pc2")])
fitg <- AtoG(fita)[[1]]
plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 300)

plot(dat_10X$ann[,c("pc1", "pc2")])
points(plot_fit, col="red")

res <-ResidualsG(as.matrix(dat_10X$ann[,c("pc1", "pc2")]), fitg)

# center the ellipse
tmp_xes <- dat_10X$ann[,"pc1"]-fitg[1,1]
tmp_yes <- dat_10X$ann[,"pc2"]-fitg[2,1]
plot(tmp_xes, tmp_yes)
# Calculate angle of each cell.
angles <- atan2(tmp_yes,tmp_xes)
angles_orig <- atan2(tmp_yes,tmp_xes)
#tmp_col <- brewer.pal(8,"Greys")[2:8]
tmp_col <- colorRampPalette(c("grey90", "grey10"))(100)
angles <- pi-angles
plot(tmp_xes, tmp_yes, col=tmp_col[as.numeric(cut(angles, breaks=length(tmp_col)))], pch=16)

#start_cell <- which(round(tmp_xes)==0 & round(tmp_yes)==-5 & dat_10X$ann$cell_name=="ring")
start_cell <- which(round(tmp_xes)==0 & round(tmp_yes)==-5)[5]
start_cell <- which(angles == min(angles[start_cell]))
points(tmp_xes[start_cell], tmp_yes[start_cell], col="red", pch=16)

converted <- angles-angles[start_cell]
converted[converted < 0] <- max(converted)+abs(min(converted[converted < 0]))+converted[converted < 0]

source("~/R-Scripts/Blank_plot.R")
#png("Explain_My_Pseudotime.png", width=6, height=6, units="in", res=300)
layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
par(mar=c(4,4,1,0))
plot(tmp_xes+fitg[1,1], tmp_yes+fitg[2,1], col=tmp_col[as.numeric(cut(converted, breaks=length(tmp_col)))], pch=16, xlab="PC1", ylab="PC2")
points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16)
points(fitg[1,1], fitg[2,1], col="black", pch=3, cex=5)
lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=2)
legend("bottomleft", bty="n", "Start cell", pch=16, col="red")
blank_plot()
legend("left", fill=tmp_col[seq(from=1, to=length(tmp_col), length=7)], levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
#dev.off()

plot(converted, cell_up, xlab="pseudotime", ylab="increase in transcription")
abline(v=seq(from=0, to=max(converted), length=20), lty=3)


# bin cells by pseudotime
n_pbins=20
p_breaks <- seq(from=min(converted), to=max(converted), length=n_pbins+1)
p_bins <- cut(converted, breaks=n_pbins)
p_axis_ats <- seq(from=0, to=1, length=n_pbins)
p_axis_labs <- round( (p_breaks[1:n_pbins]+p_breaks[2:(n_pbins+1)])/2, digits=1 )
tick_angles <- seq(from=min(angles_orig), to=max(angles_orig), length=n_pbins+1)
plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 360)
plot_fit1 <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1]*0.95, fitg[4,1]*0.95, fitg[5,1], 360)
plot_fit2 <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1]*1.05, fitg[4,1]*1.05, fitg[5,1], 360)
tick_pts <- seq(from=1, to=360, length=n_pbins+1)

source("~/R-Scripts/Blank_plot.R")
png("Explain_My_Pseudotime2.png", width=7.5, height=6, units="in", res=300)
layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
par(mar=c(4,4,1,0))
plot(tmp_xes+fitg[1,1], tmp_yes+fitg[2,1], col=tmp_col[as.numeric(cut(converted, breaks=length(tmp_col)))], pch=16, xlab="PC1", ylab="PC2")
points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16)
#points(fitg[1,1], fitg[2,1], col="black", pch=3, cex=4)
lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=1)
legend("bottomleft", bty="n", "Start cell", pch=16, col="red")

arrows(plot_fit1[tick_pts,1], plot_fit1[tick_pts,2], plot_fit2[tick_pts,1], plot_fit2[tick_pts,2], len=0, lwd=2)
#thing <- as.matrix(dist(rbind(c(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1]), plot_fit)))
#thing2 <- which(thing[1,] == min(thing[1,-1]))
thing2 <- plot_fit[,2]- (tmp_yes[start_cell]+fitg[2,1])
thing2 <- which(thing2 == min(thing2))
start_rim<-plot_fit[thing2,]
lines(c(plot_fit[230,1], fitg[1,1], start_rim[1]), c(plot_fit[230,2], fitg[2,1], start_rim[2]), col="black", lwd=2) 

text(fitg[1,]-0.75, fitg[2,]-1.5, expression(theta), cex=2)

blank_plot()
legend("left", fill=tmp_col[seq(from=1, to=length(tmp_col), length=7)], levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
dev.off()


# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")



# bin genes by pf peak time
n_gbins=45
g_breaks <- seq(from=min(extern[,51]), to=max(extern[,51]), length=n_gbins+1)
g_bins <- cut(extern[,51], breaks=n_gbins);
g_axis_ats <- seq(from=0, to=1, length=n_gbins)
g_axis_labs <- signif( (g_breaks[1:n_gbins]+g_breaks[2:(n_gbins+1)])/2, digits=2 )

# Scale t-rate
up_s <- t(apply(up, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
cell_up_s <- colSums(up_s)
heat_data_s <- my_row_mean_aggregate(t(up_s), g_bins);
heat_data_s <- my_row_mean_aggregate(t(heat_data_s), p_bins);

png("Fig_Scaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,3), c(2,3)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))

plot(converted, cell_up_s, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(converted)+0.22, max(converted)-0.22), xaxt="n")
abline(v=p_breaks, lty=3)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(heat_data_s), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

breaks_up_s <- seq(from=0, to=max(cell_up_s), length=7)
cell_up_s_col <- up_col[cut(cell_up_s, breaks=breaks_up_s)]
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, xlab="PC1", ylab="PC2", col=cell_up_s_col, pch=16)

dev.off()

# Unscaled t-rate
cell_up <- Matrix::colSums(up)
heat_data <- my_row_mean_aggregate(t(up), g_bins);
heat_data <- my_row_mean_aggregate(t(heat_data), p_bins);

png("Fig_UnScaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,3), c(2,3)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))

plot(converted, cell_up, xlab="pseudotime", ylab="increase in transcription", xaxt="n", xlim=c(min(converted)+0.22, max(converted)-0.22))
abline(v=p_breaks, lty=3)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(log(heat_data)), col=up_col, ylab="peak time (h)", xlab="pseudotime",  yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

breaks_up <- seq(from=0, to=max(cell_up), length=7)
cell_up_col <- up_col[cut(cell_up, breaks=breaks_up)]
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, xlab="PC1", ylab="PC2", col=cell_up_col, pch=16)

dev.off()

# For Ginny
dat_10X$ann$pseudotime<-converted
saveRDS(dat_10X, "bergei_10X_velodata.rds")
write.table(dat_10X$ann, file="ForGinnyPseudotime2.csv", sep=",")

# Pseudotime vs UMI/Gene
dat_10X <- readRDS("bergei_10X_velodata.rds")
mat <- dat_10X$smat+dat_10X$emat
ratio <- Matrix::colSums(mat)/Matrix::colSums(mat>0)

png("UmiRatioVsPseudo_Pbergei_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, ratio, xlab="Pseudotime", ylab="nUMI/nGene")
dev.off()


png("nGenesVsPseudo_Pbergei_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, Matrix::colSums(mat>0), xlab="Pseudotime", ylab="nUMI/nGene")
dev.off()
