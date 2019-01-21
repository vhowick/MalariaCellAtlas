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


png("bergei_10X_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("bergei_10X_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()

png("bergei_10X_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=dat_10X$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")
dev.off()

## Smartseq2
require("scater")
dat_SS <- readRDS("../NewBatches/IDC_RNAVelo_Counts.rds")
velo_SS <- readRDS("../NewBatches/IDC_RNAVelo_Estimates.rds")
ann_SS <- readRDS("../NewBatches/allIDC.rds")
pca_SS <- readRDS("../NewBatches/IDC_RNAVelo_PCAobj.rds")

write.table(pca_SS$epc@scores, file="ForGinnySS2VelocityPCACellCoords_bergei.csv", sep=",")

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

png("bergei_SS_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_colors, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("bergei_SS_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("bergei_SS_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(pca_SS$epc@scores[,1], pca_SS$epc@scores[,2], col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()

require("CellTypeProfiles")
require("gplots")

# fit an ellipse
require("conicfit")
fita <- EllipseDirectFit(dat_10X$ann[,c("pc1", "pc2")])
fitg <- AtoG(fita)[[1]]
plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 300)

plot(dat_10X$ann[,c("pc1", "pc2")])
points(plot_fit, col="red")

#res <-ResidualsG(as.matrix(dat_10X$ann[,c("pc1", "pc2")]), fitg)

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

# Add bulk
bulk <- read.table("~/Collaborations/MCA/For10XVelocity/BulkcorCenters.csv", sep=",", header=TRUE)
b_angles <- atan2(bulk$PC2-fitg[2,1], bulk$PC1-fitg[1,1])
b_angles <- pi-b_angles
b_conv <- b_angles-angles[start_cell]
b_conv[b_conv < 0] <- (2*pi-angles[start_cell])+angles[start_cell]+b_conv[b_conv < 0]

bulk$pseudo <- b_conv
write.table(bulk, file="bergei_bulk_plus_pseudo.csv", sep=",")



## Fancy plotting

source("~/R-Scripts/Blank_plot.R")
#png("Explain_My_Pseudotime.png", width=6, height=6, units="in", res=300)
#layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
#par(mar=c(4,4,1,0))
#plot(tmp_xes+fitg[1,1], tmp_yes+fitg[2,1], col=tmp_col[as.numeric(cut(converted, breaks=length(tmp_col)))], pch=16, xlab="PC1", ylab="PC2")
#points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16)
#points(fitg[1,1], fitg[2,1], col="black", pch=3, cex=5)
#lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=2)
#legend("bottomleft", bty="n", "Start cell", pch=16, col="red")
#text(bulk$PC1, bulk$PC2, col="purple", cex=2, bulk[,2])
#blank_plot()
#legend("left", fill=tmp_col[seq(from=1, to=length(tmp_col), length=7)], levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
#dev.off()

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

# colour by interpolated bulk
pseudo_col <- colorRampPalette(c("#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1"))(length(b_conv)+1)
pseudo_col <- colorRampPalette(c("#78c679", "#addd8e", "#d9f0a3","#fec44f", "#fee391", "#fff7bc", "#e0f3f8", "#abd9e9","#74add1", "#bcbddc","#9e9ac8", "#8c6bb1"))(length(b_conv)+1)
pseudo_col_breaks <- c(0, b_conv, 2*pi)

# Redo cluster colours with my palette
#tab_color_remap <- cbind(
#	c("mistyrose", "pink", "orchid1", "orchid3", "orchid4", "violetred1", "violetred3", "violetred4"),
#	c("#78C679", "#D1EC9F", "#F4CF63", "#FEEEAA", "#C9E8F1", "#98CAE1", "#B7B7D8","#9C96C6")
#)
#write.table(tab_color_remap, "Cluster_color_remap.csv", sep=",")
tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)

dat_10X$ann$cells_recolor <- as.character(dat_10X$ann$cell_color)
for (i in 1:nrow(tab_color_remap)) {
	dat_10X$ann$cells_recolor[dat_10X$ann$cells_recolor == tab_color_remap[i,1]] <- tab_color_remap[i,2]
}

source("~/R-Scripts/Blank_plot.R")
#png("Explain_My_Pseudotime_clusters_wbulk_bergei.png", width=7.5, height=6, units="in", res=300)
pdf("Explain_My_Pseudotime_clusters_wbulk_bergei.pdf", width=7.5, height=6)
layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
par(mar=c(4,4,1,0))
plot(tmp_xes+fitg[1,1], tmp_yes+fitg[2,1], 
	#col=pseudo_col[as.numeric(cut(converted, breaks=pseudo_col_breaks))], 
	col=dat_10X$ann$cells_recolor, 
	pch=16, xlab="PC1", ylab="PC2", bty="l")
points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16, cex=1.5)
lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=1)
legend("bottomleft", bty="n", "Start cell", pch=16, col="red")

arrows(plot_fit1[tick_pts,1], plot_fit1[tick_pts,2], plot_fit2[tick_pts,1], plot_fit2[tick_pts,2], len=0, lwd=2)
thing2 <- plot_fit[,2]- (tmp_yes[start_cell]+fitg[2,1])
thing2 <- which(thing2 == min(thing2))
start_rim<-plot_fit[thing2,]
lines(c(plot_fit[230,1], fitg[1,1], start_rim[1]), c(plot_fit[230,2], fitg[2,1], start_rim[2]), col="black", lwd=2) 

text(fitg[1,]-0.75, fitg[2,]-1.5, expression(theta), cex=2)
points(bulk$PC1, bulk$PC2, col="black", cex=1.5, pch=16)
text(bulk$PC1[-c(4,5,6,7,8,11)], bulk$PC2[-c(4,5,6,7,8,11)], col="black", cex=1.5, bulk[-c(4,5,6,7,8,11),2], font=2, pos=2)
text(bulk$PC1[c(4,5,6,8,11)], bulk$PC2[c(4,5,6,8,11)], col="black", cex=1.5, bulk[c(4,5,6,8,11),2], font=2, pos=4)
text(bulk$PC1[7], bulk$PC2[7], col="black", cex=1.5, bulk[7,2], font=2, pos=1)

blank_plot()
#legend("left", fill=tmp_col[seq(from=1, to=length(tmp_col), length=7)], levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
tmp <- sort(pseudo_col_breaks)
#legend("left", fill=pseudo_col, as.character(round((tmp[1:(length(tmp)-1)]+tmp[2:length(tmp)])/2, digits=2)), bty="n", title="Pseudotime", cex=0.75)
legend("left", fill=tab_color_remap[,2], c("early-ring", "late-ring", "early-troph", "mid-troph", "late-troph", "early-schizont", "mid-schizont", "late-schizont"), bty="n")
dev.off()

dat_10X$ann$pseudotime_cell_color <- pseudo_col[as.numeric(cut(converted, breaks=pseudo_col_breaks))]


saveRDS(dat_10X, "bergei_10X_velodata_wPseudo.rds")

## Match to Published ##

# Read in published

extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]

## 10X velocyto
dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

up<- velo_10X$deltaE
up[up < 0] <- 0
cell_up <- Matrix::colSums(up)

# Match genes
up <- up[rownames(up) %in% ortho$Gene,]
extern <- extern[extern$Gene.ID %in% ortho[,4],]
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern <- extern[match(rownames(up), extern$pbanka),]


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

png("bergei_Fig_Scaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,3), c(2,3)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))

plot(converted, cell_up_s, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(converted)+0.22, max(converted)-0.22), xaxt="n", col=dat_10X$ann$cells_recolor)
abline(v=p_breaks, lty=3)
points(converted, cell_up_s,col=dat_10X$ann$cells_recolor)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(heat_data_s), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

breaks_up_s <- seq(from=0, to=max(cell_up_s), length=7)
cell_up_s_col <- up_col[cut(cell_up_s, breaks=breaks_up_s)]
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, xlab="PC1", ylab="PC2", col=cell_up_s_col, pch=16)

dev.off()
# colour scale
tmp <- heat_data_s*nrow(up_s)/nrow(heat_data_s)
my_breaks=seq(from=min(tmp), to=max(tmp), length=length(up_col)+1)
source("~/R-Scripts/Colour_bar.R")
source("~/R-Scripts/Blank_plot.R")
pdf("Bergei_TransWaves_ColourBar.pdf", width=1.75, height=5)
#blank_plot()
color.bar(up_col, ticks.lab=round(my_breaks, digits=1), add=TRUE)
title(ylab="Transcriptional Rate")
dev.off()



#### Supplementary ####
pdf("bergei_Fig_Scaled_TranscriptionalWaves_Suppl.pdf", width=8/2, height=5)
graphics::layout(mat=rbind(c(1,1), c(2,2)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))
plot(converted, cell_up_s, xlab="pseudotime", ylab="Transcription", xlim=c(min(converted)+0.22, max(converted)-0.22), xaxt="n", col=dat_10X$ann$cells_recolor)
abline(v=p_breaks, lty=3)
points(converted, cell_up_s,col=dat_10X$ann$cells_recolor)
par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(heat_data_s), ylab="Peak time (h)", xlab="Pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])
dev.off()




# Supplementary Figure
png("bergei_SupplFig_Scaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,4), c(2,3)), widths=c(5,2), heights=c(2,5))
par(mar=c(0,4,1,1))

# trans through time
plot(converted, cell_up_s, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(converted)+0.22, max(converted)-0.22), xaxt="n", col=dat_10X$ann$cells_recolor)
abline(v=p_breaks, lty=3)
points(converted, cell_up_s,col=dat_10X$ann$cells_recolor)

# Matching
par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(heat_data_s), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# Histogram
par(mar=c(4,0,1,1))
hist_vals <- table(g_bins)
barplot(hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(0, max(hist_vals)), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1)
axis(side=2, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

# blank R corner
source("~/R-Scripts/Blank_plot.R")
blank_plot()

dev.off()



# Unscaled t-rate
cell_up <- Matrix::colSums(up)
heat_data <- my_row_mean_aggregate(t(up), g_bins);
heat_data <- my_row_mean_aggregate(t(heat_data), p_bins);

png("bergei_Fig_UnScaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,3), c(2,3)), widths=c(1,1), heights=c(2,5))
par(mar=c(0,4,1,1))

plot(converted, cell_up, xlab="pseudotime", ylab="increase in transcription", xaxt="n", xlim=c(min(converted)+0.22, max(converted)-0.22), col=dat_10X$ann$cells_recolor)
abline(v=p_breaks, lty=3)
points(converted, cell_up,col=dat_10X$ann$cells_recolor)

par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(log(heat_data)), col=up_col, ylab="peak time (h)", xlab="pseudotime",  yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

breaks_up <- seq(from=0, to=max(cell_up), length=7)
cell_up_col <- up_col[cut(cell_up, breaks=breaks_up)]
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, xlab="PC1", ylab="PC2", col=cell_up_col, pch=16)

dev.off()

# Supplementary Figure
png("bergei_SupplFig_UnScaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(1,4), c(2,3)), widths=c(5,2), heights=c(2,5))
par(mar=c(0,4,1,1))

# trans through time
plot(converted, cell_up, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(converted)+0.22, max(converted)-0.22), xaxt="n", col=dat_10X$ann$cells_recolor)
abline(v=p_breaks, lty=3)
points(converted, cell_up,col=dat_10X$ann$cells_recolor)

# Matching
par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(log(heat_data)), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# Histogram
par(mar=c(4,0,1,1))
hist_vals <- table(g_bins)
barplot(hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(0, max(hist_vals)), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1)
axis(side=2, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

# blank R corner
source("~/R-Scripts/Blank_plot.R")
blank_plot()

dev.off()


# For Ginny
dat_10X$ann$pseudotime<-converted
saveRDS(dat_10X, "bergei_10X_velodata.rds")
write.table(dat_10X$ann, file="ForGinnyPseudotime2.csv", sep=",")

# Pseudotime vs UMI/Gene
dat_10X <- readRDS("bergei_10X_velodata.rds")
mat <- dat_10X$smat+dat_10X$emat
ratio <- Matrix::colSums(mat)/Matrix::colSums(mat>0)

#png("bergei_10X_vs_pseudotime.png", width=5, height=10, units="in", res=300)
pdf("bergei10X_stats_vs_pseudotime.pdf", width=5*2/3, height=10*2/3)
xes <- dat_10X$ann$pseudotime
layout(cbind(c(1,2,3,4), c(1,2,3,4)), widths=c(1,1), heights=c(1,0.9,0.9,1))
par(mar=c(0,4,3,1))
plot(xes, cell_up, xlab="Pseudotime", ylab="Transcription (unscaled)", col=dat_10X$ann$cells_recolor, xaxt="n")
thing <- smooth.spline(xes, cell_up)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(xes, cell_up_s, xlab="Pseudotime", ylab="Transcription (scaled)", col=dat_10X$ann$cells_recolor, xaxt="n")
thing <- smooth.spline(xes, cell_up_s)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(xes, ratio, xlab="Pseudotime", ylab="nUMI/nGene", col=dat_10X$ann$cells_recolor, xaxt="n")
thing <- smooth.spline(xes, ratio)
lines(thing$x, thing$y, lwd=2)
par(mar=c(4,4,0,1))
plot(xes, Matrix::colSums(mat>0), xlab="Pseudotime", ylab="nGene", col=dat_10X$ann$cells_recolor)
thing <- smooth.spline(xes, Matrix::colSums(mat>0))
lines(thing$x, thing$y, lwd=2)
dev.off()


png("UmiRatioVsPseudo_Pbergei_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, ratio, xlab="Pseudotime", ylab="nUMI/nGene")
dev.off()


png("nGenesVsPseudo_Pbergei_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, Matrix::colSums(mat>0), xlab="Pseudotime", ylab="nGene")
dev.off()
