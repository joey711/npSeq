SMALL.VAL <- 1e-8

############################################################
#			Get the down sampling matrix
############################################################
DownSam.Matrix <- function(n, cmeans)
{
	cat("Down sampling...\n")
	
	c.except <- which.min(cmeans)	# the column with smallest seq depth
	dmin <- min(cmeans)
	
	n.sam <- matrix(0, nrow(n), ncol(n))
	for (i in 1 : ncol(n))
	{
		if (i != c.except)
		{
			n.sam[, i] <- rbinom(n=nrow(n), size=n[, i], prob=dmin / cmeans[i])
		} else
		{
			n.sam[, i] <- n[, i]
		}
	}
	
	return(n.sam)
}

############################################################
#			Get the Poisson sampling matrix
############################################################
Poisson.Matrix <- function(n, cmeans)
{
	cat("Poisson sampling...\n")
	
	dbar <- exp(mean(log(cmeans)))
	
	n.sam <- matrix(0, nrow(n), ncol(n))
	for (i in 1 : ncol(n))
	{
		n.sam[, i] <- rpois(n=nrow(n), lambda=(dbar / cmeans[i]) * n[, i])
	}
	
	return(n.sam)
}

############################################################
#			Get the Resample matrix
#	meth =
#	1. Down Sampling
#	2. Poisson sampling
############################################################
Resample.Matrix <- function(n, cmeans, sam.meth)
{
	n.sam <- switch(sam.meth,
		DownSam.Matrix(n, cmeans),
		Poisson.Matrix(n, cmeans))
		
	return(n.sam)
}

############################################################
#			Calculate the stat and permuted stat
############################################################
Unpaired.Twoclass.np.permu.stat <- function(n, y, cmeans, npermu=100, sam.meth=2, 
	nsam=20, signed=F)
{
	cat("calculating the stat matrix...")
	cat("unpaired twoclass...\n")
	
	### some consts
	ng <- nrow(n)
	ns <- ncol(n)
	ns1 <- sum(y == 1)
	ns2 <- sum(y == 2)
	
	### permutation index matrix
	y.permu <- NPS.Permu.Ind.Mat(y=y, npermu=npermu, twoclass=T, paired=F)
	npermu.act <- nrow(y.permu)
	cat("actual permutation times =", npermu.act, fill=T)
	
	### all of the stat matrix
	tt <- rep(0, ng)
	ttstar0 <- matrix(0, ng, npermu.act)
	
	### begin sampling
	n.sam <- NULL
	for (i in 1 : nsam)
	{
		### generate the sampling matrix
		cat("Sampling", i, "...", fill=T)
		n.sam <- Resample.Matrix(n, cmeans, sam.meth)		
		n.sam <- n.sam + runif(n=ng*ns) * 0.1
		
		### get all the ranks
		ranks <- matrix(0, ng, ns)
		for (j in 1 : ns)
		{
			ranks[, j] <- rowSums(n.sam[, -j] < n.sam[, j])
		}
		ranks <- ranks + 1
	
		### get all the statistics
		tt <- tt + rowSums(ranks[, y == 1]) - ns1 * (ns1 + 1) / 2
		for (j in 1 : npermu.act)
		{
			ttstar0[, j] <- ttstar0[, j] + rowSums(ranks[, y.permu[j, ] == 1]) - 
				ns1 * (ns1 + 1) / 2
		}
	}
	
	### get the mean
	tt <- tt / nsam
	ttstar0 <- ttstar0 / nsam
	
	### if unsigned
	if (!signed)
	{
		tt <- abs(tt - ns1 * ns2 / 2)
		ttstar0 <- abs(ttstar0 - ns1 * ns2 / 2)
	}

	return(list(tt=tt, ttstar0=ttstar0))
}

############################################################
#			Calculate the stat and permuted stat
############################################################
Quant.np.permu.stat <- function(n, y, cmeans, npermu=100, sam.meth=2, 
	nsam=20, signed=F)
{
	cat("calculating the stat matrix...")
	cat("quantitative...\n")
	
	### some consts
	ng <- nrow(n)
	ns <- ncol(n)
	ns1 <- sum(y == 1)
	ns2 <- sum(y == 2)
	
	### ranked and centered y
	y.scaled <- scale(y) / sqrt(ns - 1)
	testway2const <- (ns * (ns ^ 2 - 1) / 12) ^ 2
	
	### permutation index matrix
	ind.permu <- NPS.Permu.Ind.Mat(y=y, npermu=npermu, twoclass=F, paired=F)
	npermu.act <- nrow(ind.permu)
	cat("actual permutation times =", npermu.act, fill=T)
	
	### all of the stat matrix
	tt <- rep(0, ng)
	ttstar0 <- matrix(0, ng, npermu.act)
	
	### begin sampling
	n.sam <- NULL
	for (i in 1 : nsam)
	{
		### generate the sampling matrix
		cat("Sampling", i, "...", fill=T)
		n.sam <- Resample.Matrix(n, cmeans, sam.meth)	
		n.sam <- n.sam + runif(n=ng*ns) * 0.1

		### get all the ranks
		ranks <- matrix(0, ng, ns)
		for (j in 1 : ns)
		{
			ranks[, j] <- rowSums(n.sam[, -j] < n.sam[, j])
		}
		ranks <- ranks + 1 - (ns + 1) / 2
		
		### get all the statistics
		y.ranked <- rank(y, ties.method="random") - (ns + 1) / 2
		tt <- tt + ranks %*% y.ranked
		for (j in 1 : npermu.act)
		{
			ttstar0[, j] <- ttstar0[, j] + ranks %*% y.ranked[ind.permu[j, ]]
		}
	}
	
	### get the mean
	tt <- tt / nsam
	ttstar0 <- ttstar0 / nsam
	
	### if unsigned
	if (!signed)
	{
		tt <- abs(tt)
		ttstar0 <- abs(ttstar0)
	}

	return(list(tt=tt, ttstar0=ttstar0))
}

############################################################
#			Calculate the stat and permuted stat
############################################################
Survi.np.permu.stat <- function(n, y, gamma, cmeans, npermu=100, 
	sam.meth=2, nsam=20, signed=F)
{
	cat("calculating the stat matrix...")
	cat("survival...\n")

	### some consts
	ng <- nrow(n)
	ns <- ncol(n)
	
	### permutation index matrix
	ind.permu <- NPS.Permu.Ind.Mat(y=y, npermu=npermu, twoclass=F, paired=F)
	npermu.act <- nrow(ind.permu)
	cat("actual permutation times =", npermu.act, fill=T)
	
	### all of the stat matrix
	tt <- rep(0, ng)
	ttstar0 <- matrix(0, ng, npermu.act)
	
	### begin sampling
	n.sam <- NULL
	for (i in 1 : nsam)
	{
		### generate the sampling matrix
		cat("Sampling", i, "...", fill=T)
		n.sam <- Resample.Matrix(n, cmeans, sam.meth)	
		n.sam <- n.sam + runif(n=ng*ns) * 0.1
		
		### get all the ranks
		ranks <- matrix(0, ng, ns)
		for (j in 1 : ns)
		{
			ranks[, j] <- rowSums(n.sam[, -j] < n.sam[, j])
		}
		ranks <- ranks + 1

		### get all the statistics
		survi <- NPS.Survi.Score(ranks, y, gamma)
		tt <- tt + survi$dev1 / (sqrt(-survi$dev2) + SMALL.VAL)
		for (j in 1 : npermu.act)
		{
			survi <- NPS.Survi.Score(ranks, y[ind.permu[j, ]], gamma[ind.permu[j, ]])
			ttstar0[, j] <- ttstar0[, j] + survi$dev1 / (sqrt(-survi$dev2) + 
				SMALL.VAL)
		}
	}
	
	### get the mean
	tt <- tt / nsam
	ttstar0 <- ttstar0 / nsam
	
	### if unsigned
	if (!signed)
	{
		tt <- abs(tt)
		ttstar0 <- abs(ttstar0)
	}

	return(list(tt=tt, ttstar0=ttstar0))
}

############################################################
#			Calculate the survival statistics
############################################################
NPS.Survi.Score <- function(mu, y, gamma)
{
	# find the index matrix
	Dn <- sum(gamma == 1)
	Dset <- c(1 : ncol(mu))[gamma == 1]		# the set of observed
	ind <- matrix(0, ncol(mu), Dn)
	
	# get the matrix
	for (i in 1 : Dn)
	{
		ind[y > y[Dset[i]] - 1e-8, i] <- 1 / sum(y > y[Dset[i]] - 1e-8)
	}
	ind.sums <- rowSums(ind)
	x.ind <- mu %*% ind
	
	# get the derivatives
	dev1 <- mu %*% (gamma - ind.sums)
	dev2 <- (mu * mu) %*% ind.sums - rowSums(x.ind * x.ind)
	
	return(list(dev1=dev1, dev2=-dev2))
}

############################################################
#			Calculate the stat and permuted stat
############################################################
Multiclass.np.permu.stat <- function(n, y, cmeans, npermu=100, 
	sam.meth=2, nsam=20)
{
	cat("calculating the stat matrix...")
	cat("multiclass...\n")
	
	### some consts
	ng <- nrow(n)
	ns <- ncol(n)
	K <- max(y)
	cat("number of classes =", K, fill=T)
	n.each <- rep(0, K)
	for (k in 1 : K)
	{
		n.each[k] <- sum(y == k)
	}
	
	### permutation index matrix
	ind.permu <- NPS.Permu.Ind.Mat(y=y, npermu=npermu, twoclass=F, paired=F)
	npermu.act <- nrow(ind.permu)
	cat("actual permutation times =", npermu.act, fill=T)
	
	### all of the stat
	tt <- rep(0, ng)
	ttstar0 <- matrix(0, ng, npermu.act)

	### begin sampling
	n.sam <- NULL
	for (i in 1 : nsam)
	{
		### generate the sampling matrix
		cat("Sampling", i, "...", fill=T)
		n.sam <- Resample.Matrix(n, cmeans, sam.meth)	
		n.sam <- n.sam + runif(n=ng*ns) * 0.1
		
		### get all the ranks
		ranks <- matrix(0, ng, ns)
		for (j in 1 : ns)
		{
			ranks[, j] <- rowSums(n.sam[, -j] < n.sam[, j])
		}
		ranks <- ranks + 1
		
		### get all the statistics
		for (k in 1 : K)
		{
			tt <- tt + (rowSums(ranks[, y == k])) ^ 2 / n.each[k]
			for (j in 1 : npermu.act)
			{
				ttstar0[, j] <- ttstar0[, j] + 
					(rowSums(ranks[, y[ind.permu[j, ]] == k])) ^ 2 / n.each[k]
			}
		}
	}
	
	tt <- tt / nsam * 12 / ns / (ns + 1) - 3 * (ns + 1)
	ttstar0 <- ttstar0 / nsam * 12 / ns / (ns + 1) - 3 * (ns + 1)
	
	return(list(tt=tt, ttstar0=ttstar0))
}
