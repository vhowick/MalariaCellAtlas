glm_fun <- function(gene_expr, pseudotime) {
	res <- glm(gene_expr~pseudotime)
	eff <- res$coef[2]
        eff <- eff*1/mean(gene_expr > 0); # adjust for not shifting zeros
        zeros <- which (gene_expr == 0);

	resid <- res$residuals
#	reg <- lm(abs(resid)~pseudotime)
#	new_resid <- abs(resid)-reg$coef[2]*pseudotime + median(pseudotime)*reg$coef[2] + reg$coef[1]
#	new_resid <- new_resid -min(new_resid)
#	new_resid <- sign(resid)*new_resid

	corrected <- mean(pseudotime)*res$coef[2]+res$coef[1] + resid
#        corrected[corrected < 0] <- 0;
	return(corrected);
}

do_corrected_HVG <- function(x,fdr_thresh=0.05, method="Bren") {
	require("M3Drop")
	require("scater")
	require("Matrix")
	count_mat <- assays(x)[["counts"]]
	sf <- Matrix::colSums(count_mat)
	norm_mat <- t(t(count_mat)/sf*median(sf))
	out <- apply(norm_mat, 1, glm_fun, pseudotime=x$pseudotime)
	out <- t(out)
	if (method=="Bren") {
		HVG <- BrenneckeGetVariableGenes(out, fdr=fdr_thresh);
	} else if (method=="M3Drop") {
		out[norm_mat ==0 | out < 0] <- 0
		HVG <- M3DropFeatureSelection(out, mt_method="fdr", mt_threshold=fdr_thresh)
	} else if (method=="myNB"){
		HVG <- my_NB_HVG(out, fdr_thresh=fdr_thresh);
	}
	return(HVG);
}
my_NB_HVG <- function(mat, fdr_thresh=0.05, suppress.plot=FALSE) {
	# v of v : https://math.stackexchange.com/questions/72975/variance-of-sample-variance
	# Moments of NB: http://mathworld.wolfram.com/NegativeBinomialDistribution.html
	require("Matrix")
	
	mu_obs <- Matrix::rowMeans(mat)
	n <- ncol(mat);
	v_obs <- Matrix::rowSums((mat-mu_obs)^2)/(n-1)
	
	# var = mu + disp * mu^2

	tmp <- mu_obs^2
	disp <- glm((v_obs-mu_obs)~tmp+0)$coef[1]
	v_fitted <- mu_obs+disp*mu_obs^2

	p <- mu_obs/v_fitted
	r <- mu_obs*p/(1-p)	


	mu4 <- r*(1-p)*(6-6*p+p^2+3*r-3*p*r)/(p^4)
	sigma2 <- r*(1-p)/(p^2)

	v_of_v <- mu4/n - (sigma2^2*(n-3)/(n*n-1))
	z <- (v_obs - sigma2)/sqrt(v_of_v)
	p <- pnorm(z, lower.tail=FALSE)
	q <- p.adjust(p, method="fdr")
        eff <- v_obs-sigma2
	tab <- data.frame(rownames(mat), eff, p, q)
	colnames(tab) <- c("Gene", "effect.size", "p.value", "q.value")
	tab <- tab[!is.na(tab$p.value),]
	tab <- tab[order(-tab$q.value, tab$effect.size, decreasing=TRUE),]
	if (!suppress.plot) {
		plot(mu_obs, v_obs, cex=0.75, pch=16, xlab="mean", ylab="variance",log="xy")
		points(mu_obs[q < fdr_thresh], v_obs[q < fdr_thresh], col="red", pch=16)
		# Lines
		reorder <- order(mu_obs)
		lines(mu_obs[reorder], sigma2[reorder], col="grey80", lwd=2, lty=1)
		lines(mu_obs[reorder], sigma2[reorder]+sqrt(v_of_v[reorder])*qnorm(fdr_thresh, lower.tail=FALSE), col="grey80", lwd=2, lty=2)
	}
	return(tab[tab$q.value < fdr_thresh,])
}
my_SS_HVG <- function(mat, fdr_thresh=0.05, suppress.plot=FALSE) {
	require("Matrix")
	
	mu_obs <- Matrix::rowMeans(mat)
	n <- ncol(mat);
	v_obs <- Matrix::rowSums((mat-mu_obs)^2)/(n-1)
	
	# var = mu + disp * mu^2

	tmp <- mu_obs^2
	disp <- glm(v_obs/mu_obs~tmp)$coef
	v_fitted <- (disp[1]+mu_obs*disp[2])*mu_obs;

	p <- mu_obs/v_fitted
	r <- mu_obs*p/(1-p)	


	mu4 <- r*(1-p)*(6-6*p+p^2+3*r-3*p*r)/(p^4)
	sigma2 <- r*(1-p)/(p^2)

	v_of_v <- mu4/n - (sigma2^2*(n-3)/(n*n-1))
	z <- (v_obs - sigma2)/sqrt(v_of_v)
	p <- pnorm(z, lower.tail=FALSE)
	q <- p.adjust(p, method="fdr")
        eff <- v_obs-sigma2
	tab <- data.frame(rownames(mat), eff, p, q)
	colnames(tab) <- c("Gene", "effect.size", "p.value", "q.value")
	tab <- tab[!is.na(tab$p.value),]
	tab <- tab[order(-tab$q.value, tab$effect.size, decreasing=TRUE),]
	if (!suppress.plot) {
		plot(mu_obs, v_obs, cex=0.75, pch=16, xlab="mean", ylab="variance",log="xy")
		points(mu_obs[q < fdr_thresh], v_obs[q < fdr_thresh], col="red", pch=16)
		# Lines
		reorder <- order(mu_obs)
		lines(mu_obs[reorder], sigma2[reorder], col="grey80", lwd=2, lty=1)
		lines(mu_obs[reorder], sigma2[reorder]+sqrt(v_of_v[reorder])*qnorm(fdr_thresh, lower.tail=FALSE), col="grey80", lwd=2, lty=2)
	}
	return(tab[tab$q.value < fdr_thresh,])
}


my_Pois_HVG <- function(mat, fdr_thresh=0.05, suppress.plot=FALSE) { # Not used
	require("Matrix")
	lambdas <- Matrix::rowMeans(mat)
	n <- ncol(mat);
	v_obs <- Matrix::rowSums((mat-lambdas)^2)/(n-1)
	v_of_v <- lambdas*(1+3*lambdas)/n - (lambdas^2*(n-3)/(n*n-1))
	# v of v = u4/n - sigma^4/n * (n-3)/(n-1)
	z <- (v_obs - lambdas)/sqrt(v_of_v)
	p <- pnorm(z, lower.tail=FALSE)
	q <- p.adjust(p, method="fdr")
        eff <- v_obs-lambdas
	tab <- data.frame(rownames(mat), eff, p, q)
	colnames(tab) <- c("Gene", "effect.size", "p.value", "q.value")
	tab <- tab[!is.na(tab$p.value),]
	tab <- tab[order(-tab$q.value, tab$effect.size, decreasing=TRUE),]
	if (!suppress.plot) {
		plot(lambdas, v_obs, cex=0.75, pch=16, xlab="mean", ylab="variance",log="xy")
		abline(a=0, b=1, col="grey80", lwd=2, lty=1)
		points(lambdas[q < fdr_thresh], v_obs[q < fdr_thresh], col="red", pch=16)
		reorder <- order(lambdas)
		lines(lambdas[reorder], lambdas[reorder]+sqrt(v_of_v[reorder])*qnorm(fdr_thresh, lower.tail=FALSE), col="grey80", lwd=2, lty=2)
	}
	return(tab[tab$q.value < fdr_thresh,])
}


### 10X ###
dat10x <- readRDS("/lustre/scratch118/malaria/team222/ginny/MCA/PBANKA10xIDC_20181217.rds")
dat10x <- dat10x[Matrix::rowSums(dat10x@assays[["counts"]]) > 0,]
groups <- unique(colData(dat10x)$absclust)
MAT <- matrix(1, nrow=nrow(dat10x), ncol=length(groups));
for (i in 1:length(groups)) {
	g <- groups[i]
	tmp_dat <- dat10x[,colData(dat10x)$absclust ==g]
	tmp_dat <- tmp_dat[Matrix::rowSums(tmp_dat@assays[["counts"]] > 0) > 0,]
	out <- do_corrected_HVG(tmp_dat, fdr_thresh=0.05, method="myNB")
	MAT[,i] <- out[match(rownames(dat10x), out$Gene),"q.value"]
}

colnames(MAT) <- groups
rownames(MAT) <- rownames(dat10x)

write.table(MAT, file="MCA_10x_myHVG.csv", sep="\t")


### SS2 ###

dat <- readRDS("/lustre/scratch118/malaria/team222/ginny/MCA/MCAqcTMMSLS34_20181026.rds")
ann <- read.table("~/Collaborations/MCA/SS2_Info_allpptinfo_20180629.csv", header=T, sep=",")
dat$pseudotime <- ann[match(colnames(dat), ann[,1]),2]

groups <- unique(colData(dat)$ShortenedLifeStage4)
MAT <- matrix(1, nrow=nrow(dat), ncol=length(groups));
png("MCA_SS2_m3dropHVG.png", width=7, height=7, units="in", res=300)
par(mfrow=c(3,4))
par(mar=c(3,3,2,1))
for (i in 1:length(groups)) {
        g <- groups[i]
        tmp_dat <- dat[,colData(dat)$ShortenedLifeStage4 ==g]
        tmp_dat <- tmp_dat[Matrix::rowSums(tmp_dat@assays[["counts"]] > 0) > 0,]
        out <- do_corrected_HVG(tmp_dat, fdr_thresh=0.05, method="M3Drop")
	title(main=g)
        MAT[,i] <- out[match(rownames(dat), out$Gene),"q.value"]
}

colnames(MAT) <- groups
rownames(MAT) <- rownames(dat)
dev.off()

write.table(MAT, file="MCA_SS2_M3Drop_HVG.csv", sep="\t")

