# Examining transcriptional waves using RNA velocity through IDC
get_delta_stats <- function(deltaE) {
	tmp <- deltaE;
	tmp[tmp < 0] <- 0
	cell_up <- Matrix::colSums(tmp)
	up_mat <- tmp;
	tmp<- deltaE
	tmp[tmp > 0] <- 0
	cell_down <- Matrix::colSums(tmp)
	tmp<- deltaE
	cell_net <- Matrix::colSums(tmp)
	cell_mag <- sqrt(Matrix::colSums(tmp^2))
	up_mat_s <- t(apply(up_mat, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
	return(list(cell_up=cell_up, cell_dw=cell_down, cell_net=cell_net, cell_mag=cell_mag, up_mat=up_mat, up_mat_scaled=up_mat_s))
}

dir_pca <- function(stats, pcs, dir="up") {
	require("RColorBrewer")
	if (dir=="up") {
		col <- brewer.pal(8, "Reds")
		my_breaks <- seq(from=0, to=max(stats$cell_up), length=7)
		vals <- stats$cell_up
	} else if (dir=="down") {
		col <- rev(brewer.pal(8, "Blues"))
		my_breaks <- seq(from=min(stats$cell_dw), to=0, length=7)
		vals <- stats$cell_dw


	} else if (dir == "up_scaled") {
		col <- brewer.pal(8, "Reds")
		vals <- Matrix::colSums(stats$up_mat_scaled)
		my_breaks <- seq(from=0, to=max(vals), length=7)
	}

	cell_col <- col[cut(vals, breaks=my_breaks)]
	plot(pcs[,1], pcs[,2], col=cell_col, pch=16, xlab="PC1", ylab="PC2")

}

# External Data
extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]

extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern$pknowl <- ortho[match(extern$Gene.ID, ortho[,4]), 6]
extern$pfalci <- ortho[match(extern$Gene.ID, ortho[,4]), 4]
extern <- extern[!is.na(extern$pbanka),]
rm(ortho)

match_external <- function(spp=c("pbanka", "pknowl", "pfalci"), mat) {
	mat <- mat[rownames(mat) %in% extern[,spp[1]],]
	ext <- extern[match(rownames(mat), extern[,spp[1]]),]
	return(list(mat=mat, ext=ext))
}


## Ellipse Pseudotime
require("CellTypeProfiles")
require("gplots")

fit_elliptical_pseudo <- function(pcs, start_cell, suppress.plot=FALSE, clockwise=TRUE, nbins=20) {
	require("conicfit")
	fita <- EllipseDirectFit(pcs[,1:2])
	fitg <- AtoG(fita)[[1]]
	plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 360)
	ticks_inner <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1]*0.95, fitg[4,1]*0.95, fitg[5,1], 360)
	ticks_outer <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1]*1.05, fitg[4,1]*1.05, fitg[5,1], 360)

	tmp_xes <- pcs[,1]-fitg[1,1]
	tmp_yes <- pcs[,2]-fitg[2,1]
	angles_orig <- atan2(tmp_yes,tmp_xes)
	if (clockwise) {
		angles <- pi-angles_orig
	} else {
		angles <- angles_orig+pi
	}
	converted <- angles-angles[start_cell]
        converted[converted < 0] <- max(converted)+
                abs(min(converted[converted < 0]))+
                converted[converted < 0]
	
	# Binning
	p_breaks <- seq(from=min(converted), to=max(converted), length=nbins+1)
	p_bins <- cut(converted, breaks=nbins)

	if (!suppress.plot) {
	# Plotting
	source("~/R-Scripts/Blank_plot.R")
	pseudo_col <- colorRampPalette(c("grey90", "grey10"))(100)
	tick_pts <- seq(from=1, to=360, length=nbins+1)
	layout(mat=cbind(c(1,1), c(2,2)), width=c(3.5,1), heights=c(1,1))
	par(mar=c(4,4,1,0))
	plot(pcs[,1], pcs[,2], col=pseudo_col[as.numeric(cut(converted, breaks=length(pseudo_col)))], 
		pch=16, xlab="PC1", ylab="PC2")
	points(pcs[start_cell,1], pcs[start_cell,2], col="red", pch=16)
	lines(plot_fit[,1], plot_fit[,2], col="black", pch=16, lwd=2, lty=1)
	legend("bottomleft", bty="n", "Start cell", pch=16, col="red")

	# ticks
	arrows(ticks_inner[tick_pts,1], ticks_inner[tick_pts,2], ticks_outer[tick_pts,1], ticks_outer[tick_pts,2], 
		len=0, lwd=2)
	
	# center
	points(fitg[1,1], fitg[2,1], col="black", pch=3, cex=5)
	
	# legend
	blank_plot()
	legend("left", fill=pseudo_col[seq(from=1, to=length(pseudo_col), length=7)], 
		levels(cut(converted, breaks=7)), bty="n", title="Pseudotime", cex=0.75)
	}

	return(list(pseudotime=converted, pseudo_binned=p_bins, pseudo_breaks=p_breaks, 
			angles=angles_orig, center=c(fitg[1,1], fitg[2,1])))
}

vs_heatmap <- function(stats, pseudo, scale=TRUE, ngbins=40, spp=c("pbanka", "pknowl", "pfalci")) {
	my_row_mean_aggregate <- function (mat, groups) {
		MAT <- as.matrix(mat)
		x <- split(seq(ncol(MAT)), groups)
		result <- sapply(x, function(a) Matrix::rowMeans(MAT[, a]))
		return(result)
	}

	up <- match_external(spp=spp, stats$up_mat)


	g_breaks <- seq(from=min(up$ext[,51]), to=max(up$ext[,51]), length=ngbins+1)
	g_bins <- cut(up$ext[,51], breaks=ngbins);
	g_axis_ats <- seq(from=0, to=1, length=ngbins)
	g_axis_labs <- signif( (g_breaks[1:ngbins]+g_breaks[2:(ngbins+1)])/2, digits=2 )
	
	if (scale) {
		up$mat <- t(apply(up$mat, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
	}
	cell_up <- Matrix::colSums(up$mat)
	heat_data <- my_row_mean_aggregate(t(up$mat), g_bins);
	heat_data <- my_row_mean_aggregate(t(heat_data), pseudo$pseudo_binned);

	par(mar=c(4,4,1,1))
	image(x=pseudo$pseudo_breaks, y=g_breaks, z=t(heat_data), 
		ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
	axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], 
		labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

	return(heat_data);
}


## bergei 10X
bergei_dat <- readRDS("bergei_10X_velodata.rds")
bergei_velo <- readRDS("bergei_10X_veloEst.rds")

bergei_stats <- get_delta_stats(bergei_velo$deltaE)

# PCA Plots
dir_pca(bergei_stats, bergei_dat$ann[,c("pc1", "pc2")], dir="down")
dir_pca(bergei_stats, bergei_dat$ann[,c("pc1", "pc2")], dir="up")
dir_pca(bergei_stats, bergei_dat$ann[,c("pc1", "pc2")], dir="up_scaled")

plot(bergei_dat$ann$pc1, bergei_dat$ann$pc2, col=bergei_dat$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")

# Heatmap
bergei_pseudo <- fit_elliptical_pseudo(bergei_dat$ann[,c("pc1", "pc2")], 1, clockwise=TRUE, nbins=20)
bergei_heatmap <- vs_heatmap(bergei_stats, bergei_pseudo, scale=TRUE, ngbins=40)

# Scatter Plots
plot(bergei_pseudo$pseudotime, bergei_stats$cell_up)

## knowlesi 10X
knowlesi_dat <- readRDS("knowlesi_10X_velodata.rds")
knowlesi_velo <- readRDS("knowlesi_10X_veloEst.rds")

knowlesi_stats <- get_delta_stats(knowlesi_velo$deltaE)

dir_pca(knowlesi_stats, knowlesi_dat$ann[,c("pc1", "pc2")], dir="down")
dir_pca(knowlesi_stats, knowlesi_dat$ann[,c("pc1", "pc2")], dir="up")

plot(knowlesi_dat$ann$pc1, knowlesi_dat$ann$pc2, col=knowlesi_dat$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")

knowlesi_pseudo <- fit_elliptical_pseudo(knowlesi_dat$ann[,c("pc1", "pc2")], 1, clockwise=TRUE, nbins=20)
rownames(knowlesi_stats$up_mat) <- sub("PKNH", "PKH", rownames(knowlesi_stats$up_mat))
knowlesi_heatmap <- vs_heatmap(knowlesi_stats, knowlesi_pseudo, scale=TRUE, ngbins=40, spp="pknowl")

## faciparum 10X
falciparum_dat <- readRDS("falciparum_10X_velodata.rds")
falciparum_velo <- readRDS("falciparum_10X_veloEst.rds")

falciparum_stats <- get_delta_stats(falciparum_velo$deltaE)

dir_pca(falciparum_stats, falciparum_dat$ann[,c("pc1", "pc2")], dir="down")
dir_pca(falciparum_stats, falciparum_dat$ann[,c("pc1", "pc2")], dir="up")

plot(falciparum_dat$ann$pc1, falciparum_dat$ann$pc2, col=falciparum_dat$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")

falciparum_pseudo <- fit_elliptical_pseudo(falciparum_dat$ann[,c("pc1", "pc2")], 1, clockwise=TRUE, nbins=20)
falciparum_heatmap <- vs_heatmap(falciparum_stats, falciparum_pseudo, scale=TRUE, ngbins=40, spp="pfalci")

