############################################################
#		All Permutation indexes
#	Now, the return values are the value of y's,
#	not the permuation index
############################################################
Unpaired.Twoclass.Permu.Ind.Mat <- function(y, npermu)
{	
	library("combinat")
	
	ns <- length(y)
	ns1 <- sum(y == 1)
	ns2 <- sum(y == 2)
	
	### possible number of permutations
	poss <- choose(ns, ns1)
	
	### generate the indices
	res <- NULL
	if (npermu < poss)
	{
		res <- matrix(0, npermu, ns)
		for (i in 1 : npermu)
		{
			res[i, ] <- y[sample(c(1:ns), ns, replace=F)]
		}
	} else
	{
		ind <- t(combn(1 : ns, ns1))
		
		res <- matrix(2, nrow(ind), ns)
		for (i in 1 : nrow(ind))
		{
			res[i, ind[i, ]] <- 1
		}
	}
	
	return(res)
}

############################################################
#		All Permutation indexes
#	Note: the class labels must be like 1, 2, 1, 2, ..., 1, 2.
#	and the first two is a pair, the third and fourth is a pair, ...
############################################################
Paired.Twoclass.Permu.Ind.Mat <- function(y, npermu)
{	
	library("combinat")

	npairs <- length(y) / 2
	
	if (2^npairs > npermu)	## random samples
	{
		res <- matrix(0, npermu, npairs * 2)
		
		for (i in 1 : npermu)
		{
			ran <- sample(c(0, 1), npairs, replace=T)
			res[i, (1 : npairs) * 2 - 1] <- y[(1 : npairs) * 2 - ran]
			res[i, (1 : npairs) * 2] <- y[(1 : npairs) * 2 - (1 - ran)]
		}
	} else	## all samples
	{
		res <- matrix(0, 2^npairs, npairs * 2)
		ran <- hcube(rep(2, npairs)) - 1
		for (i in 1 : (2 ^ npairs))
		{
			res[i, (1 : npairs) * 2 - 1] <- y[(1 : npairs) * 2 - ran[i, ]]
			res[i, (1 : npairs) * 2] <- y[(1 : npairs) * 2 - (1 - ran[i, ])]
		}
	}
	
	return(res)
}

############################################################
#		All Permutation indexes
#	Here, the return value are the indexes rather than y values.
############################################################
Other.Permu.Ind.Mat <- function(y, npermu)
{
	library("combinat")
	
	res <- NULL

	ns <- length(y)
	
	if (factorial(ns) > npermu)	## random samples
	{
		res <- matrix(0, npermu, ns)
		for (i in 1 : npermu)
		{
			res[i, ] <- sample(c(1:ns), ns, replace=F)
		}
	} else	## all samples
	{
		res <- matrix(0, factorial(ns), ns)
		ran <- permn(1 : ns)
		for (i in 1 : factorial(ns))
		{
			res[i, ] <- ran[[i]]
		}
	}

	return(res)
}

############################################################
#		All Permutation indexes
#	For twoclass cases, either paired or not, 
#	the return values are the y's.
#	For other cases (multiclass, quantitative, or survival),
#	the return values are the permutation indexes.
############################################################
NPS.Permu.Ind.Mat <- function(y, npermu, twoclass=F, paired=F)
{
	res <- NULL
	if (twoclass)
	{
		if (paired)
		{
			res <- Paired.Twoclass.Permu.Ind.Mat(y, npermu)
		} else
		{
			res <- Unpaired.Twoclass.Permu.Ind.Mat(y, npermu)
		}
	} else
	{
		res <- Other.Permu.Ind.Mat(y, npermu)
	}
	
	return(res)
}
