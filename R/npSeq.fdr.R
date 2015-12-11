SMALL.VAL <- 1e-8

############################################################
#		Get the deltas
#	Note: tt must be unsigned, and sorted from
#	small to large beforehand.
#	Note: This function has been updated on Jan 1, 2011
#				to make dels always smaller than the true values.
############################################################
Get.Dels <- function(tt, nvals)
{
	ng <- length(tt)
	
	### get the breaks
	if (ng <= nvals)
	{
		dels <- tt - SMALL.VAL
	} else
	{
		dels <- tt[floor(ng / nvals * ((1 : nvals) - 1)) + 1] - SMALL.VAL
	}
	nv.act <- length(dels)

	### make the dels increasing
	dels <- sort(dels)
	for (i in nv.act : 2)
	{
		if (dels[i - 1] > dels[i] - SMALL.VAL)
		{
			dels[i - 1] <- dels[i] - SMALL.VAL
		}
	}
	
	return(dels)
}

############################################################
#		Get false discovery genes
############################################################
SAMseq.FDR.Table <- function(ssr.obj, nvals)
{
	### order the statistics
	tt <- sort(abs(ssr.obj$tt))	# from most insig to most sig
	ng <- length(tt)
	
	### get the dels
	dels <- Get.Dels(tt, nvals)
	nv.act <- length(dels)
	
	### get the nc of the original data
	nc <- ng - (rank(c(dels, tt), ties.method="min")[1 : nv.act] - 1 : nv.act)
	nc[nc > ng] <- ng
	nc[nc < 0] <- 0
	
	### get the breaks
	npermu.act <- ncol(ssr.obj$ttstar0)
	ord <- order(abs(ssr.obj$tt), decreasing=T)	# from most sig to most insig
	ttstar <- abs(ssr.obj$ttstar0)[ord, ]
	fdr <- pval <- rep(0, nv.act)	
	fdr.mat <- matrix(0, nv.act, npermu.act)
	
	############################################################
	### estimate FDR
	### estimate pi0
	thrd <- quantile(ttstar, 0.5)
	cat("\t\tthrd =", thrd, fill=T)
	pi0 <- sum(tt <= thrd) / (0.5 * ng)
	pi0 <- min(pi0, 1)
	
	### get the fdr
	for (i in 1 : npermu.act)
	{
		fdr.mat[, i] <- ng - 
			(rank(c(dels, ttstar[, i]), ties.method="min")[1 : nv.act] - 
			1 : nv.act)
	}
	num.extr <- rowMeans(fdr.mat)
	fdr <- num.extr / (nc + SMALL.VAL)
	fdr <- fdr * pi0
	pval <- num.extr / ng
	
	### adjust possible discrepancies
	fdr[fdr > 1] <- 1
	fdr[fdr < 0] <- 0
	
	pval[pval > 1] <- 1
	pval[pval < 0] <- 0
	
	### the results
	res <- list(
		dels = dels,			# deltas
		tt = ssr.obj$tt,	# the statistics
		sig.ord = ord,		# the order of the stats, from most sig to most insig
		nc = nc,					# number of sig-genes called
		fdr = fdr,				# fdr
		pi0 = pi0,				# pi0
		pval = pval				# permuted p-values; not here are sorted, not the original
		)

	return(res)
}

############################################################
#		Get the true False discovery rate
#	sig.ord: 	order, 
#						from most significant to most insignificant.
############################################################
Get.True.FDR <- function(delta, sig.ord, nc=NULL)
{
	if (is.null(nc))
	{
		nc <- 1 : length(sig.ord)
	}
	fdr <- rep(0, length(nc))
	
	for (i in seq_along(nc))
	{
		if (nc[i] == 0)
		{
			fdr[i] <- 0
		} else
		{
			fdr[i] <- sum(delta[sig.ord[1 : nc[i]]] == F) / nc[i]
		}
	}
	
	return(list(nc=nc, fdr=fdr))
}
