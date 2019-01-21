# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)
remap_colors <- function(old_colors) {
        for (i in 1:nrow(tab_color_remap)) {
                old_colors[old_colors == tab_color_remap[i,1]] <- tab_color_remap[i,2]
        }
        return(old_colors);
}


## 10X
dat_10X <- readRDS("new_pf_10X_velodata.rds")
velo_10X <- readRDS("new_pf_10X_veloEst.rds")

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


png("falciparum_10X_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("falciparum_10X_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()

png("falciparum_10X_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=remap_colors(dat_10X$ann$cell_color), pch=16, xlab="PC1", ylab="PC2")
dev.off()

## Ring Pseudotime ##


## 10X velocyto
dat_10X <- readRDS("new_pf_10X_velodata.rds")
velo_10X <- readRDS("new_pf_10X_veloEst.rds")

up<- velo_10X$deltaE
up[up < 0] <- 0
cell_up <- Matrix::colSums(up)

this_pcs <- read.table("~/Collaborations/MCA/For10XVelocity/pf_tmmfilt_PCS_qc3d7100_20190103.csv", sep=",", header=T)
this_pcs$id <- sub(".1", "", this_pcs$X)
this_pcs <- this_pcs[match(rownames(dat_10X$ann), this_pcs$id),]
dat_10X$ann$pc1_self <- this_pcs$PC1
dat_10X$ann$pc2_self <- this_pcs$PC2

require("CellTypeProfiles")
require("gplots")

# fit an ellipse
require("conicfit")
fita <- EllipseDirectFit(dat_10X$ann[,c("pc1_self", "pc2_self")])
fitg <- AtoG(fita)[[1]]
plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 300)

plot(dat_10X$ann[,c("pc1_self", "pc2_self")], col=dat_10X$ann$cell_color)
points(plot_fit, col="red")


# center the ellipse
tmp_xes <- dat_10X$ann[,"pc1_self"]-fitg[1,1]
tmp_yes <- dat_10X$ann[,"pc2_self"]-fitg[2,1]
plot(tmp_xes, tmp_yes)
# Calculate angle of each cell.
angles <- atan2(tmp_yes,tmp_xes)
angles_orig <- atan2(tmp_yes,tmp_xes)
#tmp_col <- brewer.pal(8,"Greys")[2:8]
tmp_col <- colorRampPalette(c("grey90", "grey10"))(100)
angles <- pi-angles
plot(tmp_xes, tmp_yes, col=tmp_col[as.numeric(cut(angles, breaks=length(tmp_col)))], pch=16)

#start_cell <- which(round(tmp_xes)==4 & round(tmp_yes)==5 & dat_10X$ann$cell_name=="ring")
#start_cell <- which(angles == min(angles[start_cell]))
start_cell <- which(dat_10X$ann$pseudotime == min(dat_10X$ann$pseudotime))
points(tmp_xes[start_cell], tmp_yes[start_cell], col="red", pch=16)

converted <- angles-angles[start_cell]
converted[converted < 0] <- max(converted)+abs(min(converted[converted < 0]))+converted[converted < 0]
converted <- 2*pi -converted

source("~/R-Scripts/Blank_plot.R")
png("Explain_My_Pseudotime_falciparum.png", width=6, height=6, units="in", res=300)
layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
par(mar=c(4,4,1,0))
plot(tmp_xes+fitg[1,1], tmp_yes+fitg[2,1], col=tmp_col[as.numeric(cut(converted, breaks=length(tmp_col)))], pch=16, xlab="PC1", ylab="PC2")
points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16)
points(fitg[1,1], fitg[2,1], col="black", pch=3, cex=5)
lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=2)
legend("bottomleft", bty="n", "Start cell", pch=16, col="red")
blank_plot()
legend("left", fill=tmp_col[seq(from=1, to=length(tmp_col), length=7)], levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
dev.off()

dat_10X$ann$clock_pseudotime <- converted
write.table(dat_10X$ann, "new_pf_clock_pseudotime_table.csv", sep=",", row.names=FALSE, col.names=TRUE)

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

# Add bulk
#bulk <- read.table("~/Collaborations/MCA/For10XVelocity/pfbulkcenters.csv", sep=",", header=TRUE)
bulk <- read.table("~/Collaborations/MCA/For10XVelocity/pfbulkcenters_20190114.csv", sep=",", header=TRUE)
b_angles <- atan2(bulk$PC2-fitg[2,1], bulk$PC1-fitg[1,1])
b_angles <- pi-b_angles
b_conv <- b_angles-angles[start_cell]
b_conv[b_conv < 0] <- (2*pi-angles[start_cell])+angles[start_cell]+b_conv[b_conv < 0]

bulk$pseudo <- b_conv
write.table(bulk, file="falciparum_bulk_plus_pseudo.csv", sep=",")


## clean-up
keep <- !(dat_10X$ann$cell_color == "black")
dat_10X$ann <- dat_10X$ann[keep,]
st_x <- tmp_xes[start_cell]
st_y <- tmp_yes[start_cell]

tmp_xes <- tmp_xes[keep]
tmp_yes <- tmp_yes[keep]

source("~/R-Scripts/Blank_plot.R")
png("new_Explain_My_Pseudotime2_falciparum.png", width=7.5, height=6, units="in", res=300)
pdf("new_Explain_My_Pseudotime2_falciparum.pdf", width=7.5, height=6)
layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
par(mar=c(4,4,1,0))
#plot(dat_10X$ann[,c("pc1", "pc2")], col=tmp_col[as.numeric(cut(converted, breaks=length(tmp_col)))], pch=16, xlab="PC1", ylab="PC2")
plot(dat_10X$ann[,c("pc1_self", "pc2_self")], col=remap_colors(dat_10X$ann$cell_color), pch=16, xlab="PC1", ylab="PC2")
#points(tmp_xes[start_cell]+fitg[1,1], tmp_yes[start_cell]+fitg[2,1], col="red", pch=16)
points(st_x+fitg[1,1], st_y+fitg[2,1], col="red", pch=16)
lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=1)
legend("bottomleft", bty="n", "Start cell", pch=16, col="red")

arrows(plot_fit1[tick_pts,1], plot_fit1[tick_pts,2], plot_fit2[tick_pts,1], plot_fit2[tick_pts,2], len=0, lwd=2)

#thing2 <- plot_fit[,2]- (tmp_yes[start_cell]+fitg[2,1])
#thing2 <- which(thing2 == min(thing2))
#start_rim<-plot_fit[thing2,]
#lines(c(plot_fit[230,1], fitg[1,1], start_rim[1]), c(plot_fit[230,2], fitg[2,1], start_rim[2]), col="black", lwd=2) 
lines(c(plot_fit[78,1], fitg[1,1], st_x+fitg[1,1]), 
      c(plot_fit[78,2], fitg[2,1], st_y+fitg[2,1]), col="black", lwd=2) 

#text(fitg[1,]+1.75, fitg[2,]+1, expression(theta), cex=2)
text(fitg[1,]+0.75, fitg[2,]+0.25, expression(theta), cex=2)
points(st_x+fitg[1,1], st_y+fitg[2,1], col="red", pch=16)

points(bulk$PC1, bulk$PC2, col="black", cex=1.5, pch=16)
#text(bulk$PC1, bulk$PC2, labels=bulk[,2], col="black", cex=1.5, font=2, pos=c(1,2,4,1))
text(bulk$PC1, bulk$PC2, labels=bulk[,2], col="black", cex=1.5, font=2, pos=c(2,2,2,1))


blank_plot()
legend("left", fill=tab_color_remap[,2], c("early-ring", "late-ring", "early-troph", "mid-troph", "late-troph", "early-schizont", "mid-schizont", "late-schizont"), bty="n")
dev.off()


# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

## Match to Published ##

# Read in published

extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]
extern <- extern[extern$Gene.ID %in% ortho[,4],]

# Match genes
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern$pknowl <- ortho[match(extern$Gene.ID, ortho[,4]), 6]
extern$pfalci <- ortho[match(extern$Gene.ID, ortho[,4]), 4]
extern$pknowl <- sub("PKH", "PKNH", extern$pknowl)
extern$pknowl <- paste(extern$pknowl, "0", sep="")
up <- up[rownames(up) %in% extern$pfalci,]
extern <- extern[match(rownames(up), extern$pfalci),]

# bin genes by pf peak time
n_gbins=30
g_breaks <- seq(from=min(extern[,51]), to=max(extern[,51]), length=n_gbins+1)
g_bins <- cut(extern[,51], breaks=n_gbins);
g_axis_ats <- seq(from=0, to=1, length=n_gbins)
g_axis_labs <- signif( (g_breaks[1:n_gbins]+g_breaks[2:(n_gbins+1)])/2, digits=2 )

# Scale t-rate
up_s <- t(apply(up, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
cell_up_s <- colSums(up_s)
heat_data_s <- my_row_mean_aggregate(t(up_s), g_bins);
heat_data_s <- my_row_mean_aggregate(t(heat_data_s), p_bins);

png("Fig_falciparum_Scaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
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

png("Fig_falciparum_UnScaled_TranscriptionalWaves_10X.png", width=8, height=5, units="in", res=300)
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
saveRDS(dat_10X, "falciparum_10X_velodata.rds")
write.table(dat_10X$ann, file="ForGinnyPseudotime_falciparum.csv", sep=",")

# Pseudotime vs UMI/Gene
dat_10X <- readRDS("falciparum_10X_velodata.rds")
mat <- dat_10X$smat+dat_10X$emat
ratio <- Matrix::colSums(mat)/Matrix::colSums(mat>0)
velo_10X <- readRDS("falciparum_10X_veloEst.rds")
up<- velo_10X$deltaE
up[up < 0] <- 0
cell_up <- Matrix::colSums(tmp)
up_s <- t(apply(up, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
cell_up_s <- colSums(up_s)
falciparum_scmap_pseudo <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/pfscmapclusts2method_20181122.csv", sep=",", header=T)
falciparum_scmap_pseudo$cell_id <- sub(".1", "", falciparum_scmap_pseudo$cell_id)
falciparum_scmap_pseudo <- falciparum_scmap_pseudo[match(colnames(dat_10X$emat), falciparum_scmap_pseudo$cell_id),]


png("falciparum_10X_vs_pseudotime.png", width=5, height=10, units="in", res=300)
xes <- falciparum_scmap_pseudo$pseudotime
layout(cbind(c(1,2,3,4), c(1,2,3,4)), widths=c(1,1), heights=c(1,0.9,0.9,1))
par(mar=c(0,4,3,1))
plot(xes, cell_up, xlab="Pseudotime", ylab="unscaled", col=remap_colors(dat_10X$ann$cell_color), xaxt="n")
thing <- smooth.spline(xes, cell_up)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(xes, cell_up_s, xlab="Pseudotime", ylab="scaled", col=remap_colors(dat_10X$ann$cell_color), xaxt="n")
thing <- smooth.spline(xes, cell_up_s)
lines(thing$x, thing$y, lwd=2)
par(mar=c(0,4,0,1))
plot(xes, ratio, xlab="Pseudotime", ylab="nUMI/nGene", col=remap_colors(dat_10X$ann$cell_color), xaxt="n")
thing <- smooth.spline(xes, ratio)
lines(thing$x, thing$y, lwd=2)
par(mar=c(4,4,0,1))
plot(xes, Matrix::colSums(mat>0), xlab="Pseudotime", ylab="nGene", col=remap_colors(dat_10X$ann$cell_color))
thing <- smooth.spline(xes, Matrix::colSums(mat>0))
lines(thing$x, thing$y, lwd=2)
dev.off()


png("UmiRatioVsPseudo_Pfalciparum_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, ratio, xlab="Pseudotime", ylab="nUMI/nGene")
dev.off()


png("nGenesVsPseudo_Pfalciparum_10X.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pseudotime, Matrix::colSums(mat>0), xlab="Pseudotime", ylab="nGene")
dev.off()
