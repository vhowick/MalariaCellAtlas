# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")

tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)
remap_colors <- function(old_colors) {
	for (i in 1:nrow(tab_color_remap)) {
        	old_colors[old_colors == tab_color_remap[i,1]] <- tab_color_remap[i,2]
	}
	return(old_colors);
}


## 10X
bergei_dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
bergei_velo_10X <- readRDS("bergei_10X_veloEst.rds")

tmp<- bergei_velo_10X$deltaE
tmp[tmp < 0] <- 0
bergei_UP <- tmp
bergei_cell_up <- Matrix::colSums(tmp)
bergei_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
bergei_colors <- bergei_dat_10X$ann$cell_color
bergei_pseudo <- bergei_dat_10X$ann$pseudotime


knowlesi_dat_10X <- readRDS("new_pk_10X_velodata.rds")
knowlesi_velo_10X <- readRDS("new_pk_10X_veloEst.rds")

tmp<- knowlesi_velo_10X$deltaE
tmp[tmp < 0] <- 0
knowlesi_cell_up <- Matrix::colSums(tmp)
knowlesi_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
knowlesi_colors <- knowlesi_dat_10X$ann$cell_color
knowlesi_pseudo <- knowlesi_dat_10X$ann$pseudotime
#knowlesi_scmap_pseudo <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/pkscmapclusts2method_20181122.csv", sep=",", header=T)
#knowlesi_scmap_pseudo <- knowlesi_scmap_pseudo[match(colnames(knowlesi_dat_10X$emat), knowlesi_scmap_pseudo$cell_id),]
knowlesi_scmap_pseudo <- knowlesi_dat_10X$ann
knowl_cleanup <- knowlesi_colors != "black"
knowlesi_scmap_pseudo <- knowlesi_scmap_pseudo[knowl_cleanup,]
knowlesi_colors <- knowlesi_colors[knowl_cleanup]
knowlesi_cell_up <- knowlesi_cell_up[knowl_cleanup]
knowlesi_up_s <- knowlesi_up_s[knowl_cleanup]


# Bergei SCE
#require("Matrix")

#falciparum_dat_10X <- readRDS("falciparum_10X_velodata.rds")
#falciparum_velo_10X <- readRDS("falciparum_10X_veloEst.rds")
falciparum_dat_10X <- readRDS("new_pf_10X_velodata.rds")
falciparum_velo_10X <- readRDS("new_pf_10X_veloEst.rds")

tmp<- falciparum_velo_10X$deltaE
tmp[tmp < 0] <- 0
falciparum_cell_up <- Matrix::colSums(tmp)
falciparum_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
falciparum_colors <- falciparum_dat_10X$ann$cell_color
falciparum_pseudo <- falciparum_dat_10X$ann$pseudotime
#falciparum_scmap_pseudo <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/pfscmapclusts2method_20181122.csv", sep=",", header=T)
#falciparum_scmap_pseudo$cell_id <- sub(".1", "", falciparum_scmap_pseudo$cell_id)
#falciparum_scmap_pseudo <- falciparum_scmap_pseudo[match(colnames(falciparum_dat_10X$emat), falciparum_scmap_pseudo$cell_id),]
falciparum_scmap_pseudo <- falciparum_dat_10X$ann
falcip_cleanup <- falciparum_colors != "black"
falciparum_scmap_pseudo <- falciparum_scmap_pseudo[falcip_cleanup,]
falciparum_colors <- falciparum_colors[falcip_cleanup]
falciparum_cell_up <- falciparum_cell_up[falcip_cleanup]
falciparum_up_s <- falciparum_up_s[falcip_cleanup]


# Bergei SCE
#require("Matrix")
#require("scater")
#count_mat <- readMM("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pberghei_at10_20170115/filtered_gene_bc_matrices/Pberghei_at10_20170115/matrix.mtx");
#g <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pberghei_at10_20170115/filtered_gene_bc_matrices/Pberghei_at10_20170115/genes.tsv")
#c <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pberghei_at10_20170115/filtered_gene_bc_matrices/Pberghei_at10_20170115/barcodes.tsv")
#c <- sub("-1", "", c[,1])
#g <- g[,1]

#rownames(count_mat) <- g
#colnames(count_mat) <- c

#count_mat <- count_mat[,colnames(count_mat) %in% colnames(bergei_dat_10X$emat)]

#mca_pb <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(count_mat))))

#mca_pb <- mca_pb[,match(colnames(bergei_dat_10X$emat), colnames(mca_pb))]
#pb_mat <- assays(mca_pb)[["counts"]]
pb_mat <- bergei_dat_10X$emat+bergei_dat_10X$smat
pb_ngenes_c <- Matrix::colSums(pb_mat>0)
pb_numi_g_c <- Matrix::colSums(pb_mat)/pb_ngenes_c

# Knowlesi SCE
#count_mat <- readMM("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/filtered_gene_bc_matrices/Pknowlesi_at10_20170115/matrix.mtx");
#g <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/filtered_gene_bc_matrices/Pknowlesi_at10_20170115/genes.tsv")
#c <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24726_4472STDY7251827_Pknowlesi_at10_20170115/filtered_gene_bc_matrices/Pknowlesi_at10_20170115/barcodes.tsv")
#c <- sub("-1", "", c[,1])
#g <- g[,1]

#rownames(count_mat) <- g
#colnames(count_mat) <- c

#count_mat <- count_mat[,colnames(count_mat) %in% colnames(knowlesi_dat_10X$emat)]

#mca_pk <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(count_mat))))

#mca_pk <- mca_pk[,match(colnames(knowlesi_dat_10X$emat), colnames(mca_pk))]
pk_mat <- knowlesi_dat_10X$emat+knowlesi_dat_10X$smat
pk_mat <- pk_mat[,knowl_cleanup]
pk_ngenes_c <- Matrix::colSums(pk_mat>0)
pk_numi_g_c <- Matrix::colSums(pk_mat)/pk_ngenes_c


# Falciparum SCE
#count_mat <- readMM("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24612_4472STDY7230041_3D7_Jan16v3-ensembl_37_transcriptome/filtered_gene_bc_matrices/3D7_Jan16v3/matrix.mtx");
#g <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24612_4472STDY7230041_3D7_Jan16v3-ensembl_37_transcriptome/filtered_gene_bc_matrices/3D7_Jan16v3/genes.tsv")
#c <- read.table("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/cellranger210_count_24612_4472STDY7230041_3D7_Jan16v3-ensembl_37_transcriptome/filtered_gene_bc_matrices/3D7_Jan16v3/barcodes.tsv")
#c <- sub("-1", "", c[,1])
#g <- g[,1]

#rownames(count_mat) <- g
#colnames(count_mat) <- c

#count_mat <- count_mat[,colnames(count_mat) %in% colnames(falciparum_dat_10X$emat)]

#mca_pf <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(count_mat))))

#mca_pf <- mca_pf[,match(colnames(falciparum_dat_10X$emat), colnames(mca_pf))]
pf_mat <- falciparum_dat_10X$emat+falciparum_dat_10X$smat
pf_mat <- pf_mat[,falcip_cleanup]
pf_ngenes_c <- Matrix::colSums(pf_mat>0)
pf_numi_g_c <- Matrix::colSums(pf_mat)/pf_ngenes_c




# Color Bars
# expression blocks

# Cell-specific
bergei_cell_breaks <- seq(from=min(bergei_pseudo), to=max(bergei_pseudo), length=151)
pseudo_cell_bins <- cut(bergei_pseudo, breaks=bergei_cell_breaks, include.lowest=TRUE)
pseudo_cell_bins_ord <- pseudo_cell_bins[order(bergei_pseudo)]
bar_pos <- cbind(bergei_cell_breaks[1:(length(bergei_cell_breaks)-1)],
		 rep(0, length(bergei_cell_breaks)-1),
		 bergei_cell_breaks[2:length(bergei_cell_breaks)],
		 rep(1, length(bergei_cell_breaks)-1)
		)
# painter
extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]

bergei_UP <- bergei_UP[rownames(bergei_UP) %in% ortho$Gene,]
extern <- extern[extern$Gene.ID %in% ortho[,4],]
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern <- extern[match(rownames(bergei_UP), extern$pbanka),]
bergei_UP_s <- t(apply(bergei_UP, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
peak_time <- extern[,51]
#cell_peaks <- apply(bergei_UP_s, 2, function(a) {sum(a*peak_time)/sum(a)}) # diff num genes through time

n_gbins=45
g_breaks <- seq(from=min(extern[,51]), to=max(extern[,51]), length=n_gbins+1)
g_bins <- cut(extern[,51], breaks=n_gbins);
require("CellTypeProfiles")
data_s <- t(my_row_mean_aggregate(t(bergei_UP_s), g_bins));
bin_time <- (g_breaks[1:(length(g_breaks)-1)]+g_breaks[2:(length(g_breaks))] )/2
bergei_cell_peaks <- apply(data_s, 2, function(a) {sum(a*bin_time)/sum(a)})

cell_peak_ord <- bergei_cell_peaks[order(bergei_pseudo)]
cell_peak_bar_vals <- aggregate(cell_peak_ord, list(pseudo_cell_bins_ord), mean)
poor_match <- apply(data_s, 2, max)<0.1
poor_match_bin <- aggregate(poor_match[order(bergei_pseudo)], list(pseudo_cell_bins_ord), mean)

peak_cols <- colorRampPalette(c("mistyrose", "pink", "orchid1", "orchid3", "orchid4", "violetred1", "violetred3","violetred4"))(20)
peak_cols <- colorRampPalette(c("white", "#252525"))(8)
peak_bar_cols <- peak_cols[cut(cell_peak_bar_vals[,2], breaks=length(peak_cols))]
peak_bar_cols[poor_match_bin[,2]>=0.5] <- "cadetblue1"

# alternate time assignment - correlations
require("Hmisc")
out <- rcorr(as.matrix(log2(extern[,3:50]+1)), as.matrix(log2(bergei_UP+1)), type="pearson")
thing <- out$r[1:48,-c(1:48)]
assigned <- apply(thing,2, function(a){which(a==max(a))})
assigned_r <- apply(thing,2, function(a){a[which(a==max(a))]})

assigned_ord <- assigned[order(bergei_pseudo)]
assigned_bar_vals <- aggregate(assigned_ord, list(pseudo_cell_bins_ord), mean)
assigned_bar_cols <- peak_cols[cut(assigned_bar_vals[,2], breaks=length(peak_cols))]

#binned_profile <- my_row_mean_aggregate(log(bergei_UP[,order(bergei_pseudo)]+1), pseudo_cell_bins_ord)
#out <- rcorr(as.matrix(log(extern[,3:50]+1)), as.matrix(binned_profile), type="pearson")
#thing <- out$r[1:48,-c(1:48)]
#assigned_bar_vals <- apply(thing,2, function(a){which(a==max(a))})
#assigned_r <- apply(thing,2, function(a){a[which(a==max(a))]})
#assigned_bar_cols <- peak_cols[cut(assigned_bar_vals, breaks=length(peak_cols))]



# ngene
ngene_cols <- colorRampPalette(c("white", "#252525"))(8)
#ngene_vals <- Matrix::colSums(assays(mca_pb)[["counts"]]>0)
ngene_vals <- pb_ngenes_c
ngene_ord <- ngene_vals[order(bergei_pseudo)]
ngene_bar_vals <- aggregate(ngene_ord, list(pseudo_cell_bins_ord), mean)
ngene_bar_cols <- ngene_cols[cut(ngene_bar_vals[,2], breaks=length(ngene_cols))]


# numi/gene
numi_cols <- colorRampPalette(c("white", "#252525"))(8)
#numi_vals <- Matrix::colSums(assays(mca_pb)[["counts"]]) / Matrix::colSums(assays(mca_pb)[["counts"]]>0)
numi_vals <- pb_numi_g_c
numi_ord <- numi_vals[order(bergei_pseudo)]
numi_bar_vals <- aggregate(numi_ord, list(pseudo_cell_bins_ord), mean)
numi_bar_cols <- numi_cols[cut(numi_bar_vals[,2], breaks=length(numi_cols))]

# Colour scale bars
source("~/R-Scripts/Colour_bar.R")
pdf("TransWaves_multispp_lineFigure_colorscales.pdf", width=2, height=6)
par(mfrow=c(1,3));
color.bar(peak_cols, ticks.at=seq(0,1,len=length(peak_cols)+1), 
		ticks.lab=round(seq(min(cell_peak_bar_vals[,2]), max(cell_peak_bar_vals[,2]), len=length(peak_cols)+1)), 
		title="", add=TRUE)
title(ylab="Painter(h)", line=2.5, cex.lab=1.5)

color.bar(ngene_cols, ticks.at=seq(0,1,len=length(ngene_cols)+1), 
		ticks.lab=round(seq(min(ngene_bar_vals[,2]), max(ngene_bar_vals[,2]), len=length(ngene_cols)+1)), 
		title="", add=TRUE)
title(ylab="nGene", line=2.5, cex.lab=1.5)

color.bar(numi_cols, ticks.at=seq(0,1,len=length(numi_cols)+1), 
		ticks.lab=round(seq(min(numi_bar_vals[,2]), max(numi_bar_vals[,2]), len=length(numi_cols)+1), digits=1), 
		title="", add=TRUE)
title(ylab="nUMI/g", line=2.5, cex.lab=1.5)

dev.off()

plot_bar <- function(bar.colors, bar.pos, xlim=c(0,1)) { 
	plot(1,1, col="white", xlim=c(0,6.28), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
	min_x <- min(bar.pos[,1], bar.pos[,3])
	max_x <- max(bar.pos[,1], bar.pos[,3])
	min_y <- min(bar.pos[,2], bar.pos[,4])
	max_y <- max(bar.pos[,2], bar.pos[,4])
	range <- xlim[2]-xlim[1]
	bar.pos[,1] <- (bar.pos[,1]- min_x)/(max_x-min_x)*range + xlim[1]
	bar.pos[,3] <- (bar.pos[,3]- min_x)/(max_x-min_x)*range + xlim[1]
	bar.pos[,2] <- (bar.pos[,2]- min_y)/(max_y-min_y)
	bar.pos[,4] <- (bar.pos[,4]- min_y)/(max_y-min_y)

	for (i in 1:length(bar.colors)) {
		rect(bar.pos[i,1], bar.pos[i,2], bar.pos[i,3], bar.pos[i,4], col=bar.colors[i], border=NA)
	}
}
		

# bin cells by pseudotime
n_pbins=20
p_breaks <- seq(from=min(bergei_pseudo), to=max(bergei_pseudo), length=n_pbins+1)
p_bins <- cut(bergei_pseudo, breaks=n_pbins)

nlines=1
nspp=3
my_xlim <- c(min(p_breaks), max(p_breaks))


#png("Three_Spp_Wave_Lines.png", width=7, height=9, units="in", res=300)
pdf("new_Three_Spp_Wave_Lines.pdf", width=7, height=9)
lmat <- cbind(c(1:nlines, nlines+1:nspp), c(1:nlines, nlines+1:nspp))
layout(lmat, width=c(1,1), heights=c(rep(0.1, times=nlines), rep(1, times=nspp)))

#lines
par(mar=c(0,4,0,1));
plot_bar(assigned_bar_cols, bar_pos, xlim=my_xlim)
mtext("Painter", side=2, line=-1.5, las=2)
#plot_bar(peak_bar_cols, bar_pos, xlim=my_xlim)
#mtext("Painter", side=2, line=-1.5, las=2)
#plot_bar(ngene_bar_cols, bar_pos, xlim=my_xlim)
#mtext("nGenes", side=2, line=-1.5, las=2)
#plot_bar(numi_bar_cols, bar_pos, xlim=my_xlim)
#mtext("UMI/g", side=2, line=-1.5, las=2)

#spp
par(mar=c(0,4,0.2,1))
plot(bergei_pseudo, bergei_up_s, col=remap_colors(bergei_colors), ylab="Bergei", xlab="", xaxt="n", xlim=my_xlim)
abline(v=p_breaks, lty=3, col="grey40")
points(bergei_pseudo, bergei_up_s, col=remap_colors(bergei_colors))
thing <- smooth.spline(bergei_pseudo, bergei_up_s)
lines(thing$x, thing$y, lwd=2)

par(mar=c(0,4,0,1))
plot(knowlesi_scmap_pseudo$pseudotime, knowlesi_up_s, col=remap_colors(knowlesi_colors), ylab="Knowlesi", xlab="", xaxt="n", xlim=my_xlim)
abline(v=p_breaks, lty=3, col="grey40")
points(knowlesi_scmap_pseudo$pseudotime, knowlesi_up_s, col=remap_colors(knowlesi_colors))
#bit1 <- knowlesi_scmap_pseudo$pseudotime<1
bit1 <- knowlesi_scmap_pseudo$pseudotime>-1
thing1 <- smooth.spline(knowlesi_scmap_pseudo$pseudotime[bit1], knowlesi_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)
#bit2 <- knowlesi_scmap_pseudo$pseudotime>1.5
#thing2 <- smooth.spline(knowlesi_scmap_pseudo$pseudotime[bit2], knowlesi_up_s[bit2])
#lines(thing2$x, thing2$y, lwd=2)

par(mar=c(4,4,0,1))
plot(falciparum_scmap_pseudo$pseudotime, falciparum_up_s, col=remap_colors(falciparum_colors), xlab="Pseudotime", ylab="Falciparum", xlim=my_xlim)
abline(v=p_breaks, lty=3, col="grey40")
points(falciparum_scmap_pseudo$pseudotime, falciparum_up_s, col=remap_colors(falciparum_colors))
#bit1 <- falciparum_scmap_pseudo$pseudotime > 1.5 & falciparum_scmap_pseudo$pseudotime < 5.5
bit1 <- falciparum_scmap_pseudo$pseudotime > 0.5 & falciparum_scmap_pseudo$pseudotime < 6.5
thing1 <- smooth.spline(falciparum_scmap_pseudo$pseudotime[bit1], falciparum_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)

dev.off()


#### Supplementary Figure ####

plot_pk <- function(vals, ylab="") {
	plot(knowlesi_scmap_pseudo$pseudotime, vals, 
		col=remap_colors(knowlesi_colors), ylab="", 
		xlab="", xaxt="n", xlim=my_xlim)
	mtext(side=2, line=2.3, ylab)
	abline(v=p_breaks, lty=3, col="grey40")
	points(knowlesi_scmap_pseudo$pseudotime, vals, col=remap_colors(knowlesi_colors))
#	bit1 <- knowlesi_scmap_pseudo$pseudotime<1
	bit1 <- knowlesi_scmap_pseudo$pseudotime>-1
	thing1 <- smooth.spline(knowlesi_scmap_pseudo$pseudotime[bit1], vals[bit1], spar=0.9)
	lines(thing1$x, thing1$y, lwd=2)
#	bit2 <- knowlesi_scmap_pseudo$pseudotime>1.5
#	thing2 <- smooth.spline(knowlesi_scmap_pseudo$pseudotime[bit2], vals[bit2])
#	lines(thing2$x, thing2$y, lwd=2)
}

plot_pf <- function(vals, ylab="") {
	plot(falciparum_scmap_pseudo$pseudotime, vals, 
		col=remap_colors(falciparum_colors), ylab="", 
		xlab="", xaxt="n", xlim=my_xlim)
	mtext(side=2, line=2.3, ylab)
	abline(v=p_breaks, lty=3, col="grey40")
	points(falciparum_scmap_pseudo$pseudotime, vals, col=remap_colors(falciparum_colors))
	#bit1 <- falciparum_scmap_pseudo$pseudotime > 1.5 & falciparum_scmap_pseudo$pseudotime < 5.5
	bit1 <- falciparum_scmap_pseudo$pseudotime > 0.5 & falciparum_scmap_pseudo$pseudotime < 6.5
	thing1 <- smooth.spline(falciparum_scmap_pseudo$pseudotime[bit1], vals[bit1], spar=0.9)
	lines(thing1$x, thing1$y, lwd=2)
}

plot_pb <- function(vals, ylab="") {
	plot(bergei_pseudo, vals, 
		col=remap_colors(bergei_colors), ylab="", 
		xlab="", xaxt="n", xlim=my_xlim)
	mtext(side=2, line=2.3, ylab)
	abline(v=p_breaks, lty=3, col="grey40")
	points(bergei_pseudo, vals, col=remap_colors(bergei_colors))
	thing1 <- smooth.spline(bergei_pseudo, vals, spar=0.9)
	lines(thing1$x, thing1$y, lwd=2)
}


tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)

#dat_10X$ann$cells_recolor <- as.character(dat_10X$ann$cell_color)
#for (i in 1:nrow(tab_color_remap)) {
#        dat_10X$ann$cells_recolor[dat_10X$ann$cells_recolor == tab_color_remap[i,1]] <- tab_color_remap[i,2]
#}


pdf("knowlesi_four_data_wave_suppl.pdf", width=7, height=9)
layout(matrix(rep(1:4, each=2), ncol=2, byrow=T), widths=c(1,1), heights=c(1, 0.95, 0.95, 1.3))
par(mar=c(0,4,1,1))
plot_pk(knowlesi_up_s, ylab="Scaled")
title(main="P. knowlesi")
par(mar=c(0,4,0,1))
plot_pk(knowlesi_cell_up, ylab="Unscaled")
par(mar=c(0,4,0,1))
plot_pk(pk_ngenes_c, ylab="nGenes")
par(mar=c(4,4,0,1))
plot_pk(pk_numi_g_c, ylab="nUMI/g")
mtext(side=1, line=2, "Pseudotime")
axis(1, at=0:6, labels=0:6)
dev.off()

pdf("bergei_four_data_wave_suppl.pdf", width=7, height=9)
layout(matrix(rep(1:4, each=2), ncol=2, byrow=T), widths=c(1,1), heights=c(1, 0.95, 0.95, 1.3))
par(mar=c(0,4,1,1))
plot_pb(bergei_up_s, ylab="Scaled")
title(main="P. bergei")
par(mar=c(0,4,0,1))
plot_pb(bergei_cell_up, ylab="Unscaled")
par(mar=c(0,4,0,1))
plot_pb(pb_ngenes_c, ylab="nGenes")
par(mar=c(4,4,0,1))
plot_pb(pb_numi_g_c, ylab="nUMI/g")
mtext(side=1, line=2, "Pseudotime")
axis(1, at=0:6, labels=0:6)
dev.off()

pdf("new_falciparum_four_data_wave_suppl.pdf", width=7, height=9)
layout(matrix(rep(1:4, each=2), ncol=2, byrow=T), widths=c(1,1), heights=c(1, 0.95, 0.95, 1.3))
par(mar=c(0,4,1,1))
plot_pf(falciparum_up_s, ylab="Scaled")
title(main="P. falciparum")
par(mar=c(0,4,0,1))
plot_pf(falciparum_cell_up, ylab="Unscaled")
par(mar=c(0,4,0,1))
plot_pf(pf_ngenes_c, ylab="nGenes")
par(mar=c(4,4,0,1))
plot_pf(pf_numi_g_c, ylab="nUMI/g")
mtext(side=1, line=2, "Pseudotime")
axis(1, at=0:6, labels=0:6)
dev.off()


