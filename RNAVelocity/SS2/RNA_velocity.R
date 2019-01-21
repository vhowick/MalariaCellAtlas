#R3.4 40Gb mem + 10 cores
library('methods')
args <- commandArgs(trailingOnly=TRUE) #

outprefix <- args[1];
path <- args[2]
n_cores <- args[3]

require("Matrix")
library(velocyto.R)
rowSums <- Matrix::rowSums
colSums <- Matrix::colSums

#path <- "/lustre/scratch118/malaria/team222/andy/24085"
files <- Sys.glob(paste(path,"/*.bam",sep=""));
names(files) <- gsub(".*\\/(.*)_sorted.bam","\\1",files)

sets <- rep(1:ceiling(length(files)/200), each=200)
sets <- sets[1:length(files)]

gtf <- "/lustre/scratch117/cellgen/team218/TA/genomebuilding/Pberghei_anno_RNAVelocity.out"
dat <- list(emat=vector(), iomat=vector(), smat=vector(), base.df=vector(), exons=vector(), genes=vector(), expr.lstat=vector())

for (s in unique(sets)) {
		outfile=paste(outprefix, "velocity_counts", s, ".rds", sep="_")
	if (!file.exists(outfile)) {
		dat_tmp <- read.smartseq2.bams(files[sets==s],gtf,n.cores=n_cores)
		saveRDS(dat_tmp, file=outfile)
	} #else {
		dat_tmp <- readRDS(outfile)
		# Merge results
		dat$exons <- dat_tmp$exons
		dat$genes <- dat_tmp$genes
		dat$emat <- cbind(dat$emat, dat_tmp$emat)
		dat$iomat <- cbind(dat$iomat, dat_tmp$iomat)
		dat$smat <- cbind(dat$smat, dat_tmp$smat)
		# This is not identical to doing rna_velocity on the whole dataset (for reasons I don't understand) but it is pretty close
		if (is.null(dim(dat$base.df))) {
			dat$base.df <- dat_tmp$base.df
		} else {
			new_min <- dat$base.df[, 1] >  dat_tmp$base.df[, 1]
			if (sum(new_min) > 0) {
				dat$base.df[new_min, 1] <- dat_tmp$base.df[new_min, 1]
			}
			new_max <- dat$base.df[, 2] <  dat_tmp$base.df[, 2]
			if (sum(new_max) > 0) {
				dat$base.df[new_max, 2] <- dat_tmp$base.df[new_max, 2]
			}
			new_max <- dat$base.df[, 3] <  dat_tmp$base.df[, 3]
			if (sum(new_max) > 0) {
				dat$base.df[new_max, 3] <- dat_tmp$base.df[new_max, 3]
			}
		}
		# This is not identical to doing rna_velocity on the whole dataset (for reasons I don't understand) but it is pretty close
		if (is.null(dim(dat$expr.lstat))) {
			dat$expr.lstat <- dat_tmp$expr.lstat
		} else {
			new_min <- dat$expr.lstat[, 1] >  dat_tmp$expr.lstat[, 1]
			if (sum(new_min) > 0) {
				dat$expr.lstat[new_min, 1] <- dat_tmp$expr.lstat[new_min, 1]
			}
			new_max <- dat$expr.lstat[, 2] <  dat_tmp$expr.lstat[, 2]
			if (sum(new_max) > 0) {
				dat$expr.lstat[new_max, 2] <- dat_tmp$expr.lstat[new_max, 2]
			}
			new_max <- dat$expr.lstat[, 3] <  dat_tmp$expr.lstat[, 3]
			if (sum(new_max) > 0) {
				dat$expr.lstat[new_max, 3] <- dat_tmp$expr.lstat[new_max, 3]
			}
		}
	#}
}

saveRDS(dat, file=paste(outprefix, "velocity_counts.rds", sep="_"))

