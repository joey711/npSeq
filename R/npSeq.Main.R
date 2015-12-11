############################################################
#		Estimating FDR using our np method
############################################################
npSeq.Main <- function(dat, para=list())
{	
	time_st <- Sys.time()

	### check arguments
	cat("Checking parameters of dat...", fill=T)
	dat <- npSeq.Check.Argu.Dat(dat)
	cat("Checking parameters of para...", fill=T)
	para <- npSeq.Check.Argu.Para(para, dat)

	### filter out genes that has too small counts
	dat <- npSeq.PS.Filter(dat, ct.sum=para$ct.sum, ct.mean=para$ct.mean)

	cmeans <- npSeq.estimate.depth(dat$n)
	
	cat("Estimating fdr...")
	
	set.seed(para$seed)
	
	np.obj <- NULL
	if (dat$type == 'twoclass')
	{
		cat("twoclass data...")
		np.obj <- Unpaired.Twoclass.np.permu.stat(n=dat$n, y=dat$y, cmeans=cmeans, 
				npermu=para$npermu, sam.meth=para$sam.meth, nsam=para$nsam, signed=F)
	} else if (dat$type == 'quant')
	{
		cat("quant data...\n")
		np.obj <- Quant.np.permu.stat(n=dat$n, y=dat$y, 
			cmeans=cmeans, npermu=para$npermu, sam.meth=para$sam.meth, nsam=para$nsam, signed=F)
	} else if (dat$type == 'survi')
	{
		cat("survival data...\n")
		np.obj <- Survi.np.permu.stat(n=dat$n, y=dat$y, gamma=dat$gamma, 
			cmeans=cmeans, npermu=para$npermu, sam.meth=para$sam.meth, nsam=para$nsam, 
			signed=F)
	} else
	{
		cat("multiclass data...\n")
		np.obj <- Multiclass.np.permu.stat(n=dat$n, y=dat$y, 
			cmeans=cmeans, npermu=para$npermu, sam.meth=para$sam.meth, nsam=para$nsam)
	}
	
	res <- SAMseq.FDR.Table(ssr.obj=np.obj, nvals=nrow(dat$n))
	
	if (!is.null(dat$delta))
	{
		res$tfdr <- Get.True.FDR(delta=dat$delta, sig.ord=res$sig.ord, 
			nc=res$nc)$fdr
	}
	
	res <- npSeq.Sum(dat, res, cmeans)
	
	time_ed <- Sys.time()
	print(time_ed - time_st)

	return(res)
}


############################################################
# check arguments of dat
############################################################
npSeq.Check.Argu.Dat <- function(dat)
{
	if (is.null(dat$n))
	{
		stop("the data matrix dat$n must be specified!")
	}
	
	if (is.null(dat$y))
	{
		stop("the outcome vector dat$y must be specified!")
	}
	
	if (is.null(dat$type))
	{
		stop("the outcome type dat$type must be specified!")
	}
	
	dat$n <- as.matrix(dat$n)
	if (!is.matrix(dat$n) | !is.numeric(dat$n))
	{
		stop('n must be a numeric matrix!')
	}
	
	if (sum(dat$n >= 0) != length(dat$n))
	{
		stop('All elements of n must be non-negative!')
	}
	
	if (!is.numeric(dat$y) | length(dat$y) <= 1)
	{
		stop('y must be a numeric vector with length no less than 2!')
	}
	
	if (ncol(dat$n) != length(dat$y))
	{
		stop('Number of column of n should be the same as the length of y!')
	}
	
	dat$type <- match.arg(dat$type, c('twoclass', 'multiclass', 'quant', 'survi'), several.ok = FALSE)
	
	if (is.null(dat$gname))
	{
		dat$gname <- 1 : nrow(dat$n)
	}
	
	if (length(dat$gname) != nrow(dat$n))
	{
		stop("Length of dat$gname must equal the number of genes nrow(dat$n)!")
	}
	
	return(dat)
}


############################################################
# check arguments of para
############################################################
npSeq.Check.Argu.Para <- function(para, dat)
{
	if (is.null(para$npermu))
	{
		para$npermu <- 100
	} else
	{
		if (!is.numeric(para$npermu) | para$npermu <= 0)
		{
			stop('npermu must be a positive integer!')
		}
	}

	if (is.null(para$nsam))
	{
		para$nsam <- 20
	} else
	{
		if (!is.numeric(para$nsam) | para$nsam <= 0)
		{
			stop('nsam must be a positive integer!')
		}
	}
	
	if (is.null(para$sam.meth))
	{
		para$sam.meth <- 2
	} else
	{
		if ((para$sam.meth != 1) & (para$sam.meth != 2))
		{
			stop('sam.meth must be either 1 or 2!')
		}
	}
	
	if (is.null(para$ct.sum))
	{
		para$ct.sum <- 5
	}
	
	if (!is.numeric(para$ct.sum))
	{
		stop("para$ct.sum must be numeric!")
	}
	
	if (is.null(para$ct.mean))
	{
		para$ct.mean <- 0.5
	}
	
	if (!is.numeric(para$ct.mean))
	{
		stop("para$ct.mean must be numeric!")
	}
	
	if (is.null(para$seed))
	{
		para$seed <- 10
	}
	
	return(para)
}

############################################################
#		Filter genes with too small counts
############################################################
npSeq.PS.Filter <- function(dat, ct.sum=5, ct.mean=0.5)
{
	if (is.null(dat$gname))
	{
		dat$gname <- 1 : nrow(dat$n)
	}
	
	keep <- (rowMeans(dat$n) > ct.mean) & (rowSums(dat$n) > ct.sum)
	cat(length(keep) - sum(keep), "genes has been filtered because they contains too small number of reads across the experiments.", fill=T)
	
	dat$n <- dat$n[keep, ]
	dat$gname <- dat$gname[keep]
	
	return(dat)
}

############################################################
#		Summarize results
############################################################
npSeq.Sum <- function(dat, fdr.res, seq.depth)
{
	### sort nc and fdr
	tt <- sort(abs(fdr.res$tt), decreasing=T)
	pval <- sort(fdr.res$pval)
	
	nc <- fdr.res$nc
	ord.fdr <- fdr.res$fdr
	
	tfdr <- NULL
	if (!is.null(fdr.res$tfdr))
	{
		tfdr <- fdr.res$tfdr
	}
	if (nc[length(nc)] < nc[1])
	{
		nc <- rev(nc)
		ord.fdr <- rev(ord.fdr)
		if (!is.null(tfdr))
		{
			tfdr <- rev(tfdr)
		}
	}
		
	### make sure the fdr is monotone increasing
	fdr <- rep(min(c(ord.fdr[length(ord.fdr)], 1)), length(ord.fdr))
	for (i in (length(fdr) - 1) : 1)
	{
		fdr[i] <- min(c(fdr[i + 1], ord.fdr[i], 1))
	}
	
	### other elements
	gname <- dat$gname[fdr.res$sig.ord]
	
	res.table <- data.frame(nc=nc, gname=gname, tt=tt, pval=pval, fdr=fdr)
	
	### get fold change for two class data
	if (dat$type == "twoclass")
	{
		n.norm <- scale(dat$n, center=F, scale=seq.depth)
		
		exp.1 <- apply(n.norm[, dat$y == 1], 1, median)
		exp.2 <- apply(n.norm[, dat$y == 2], 1, median)
		
		log.fc <- log(exp.2 / exp.1)
		log.fc[(exp.1 == 0) & (exp.2 == 0)] <- 0
		res.table <- cbind(res.table, log.fc=log.fc[fdr.res$sig.ord])
	}
	
	if (!is.null(tfdr))
	{
		res.table <- cbind(res.table, tfdr=tfdr)
	}
	
	return(res.table)
}
