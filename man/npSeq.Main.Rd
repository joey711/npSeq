\name{npSeq.Main}
\alias{npSeq.Main}
\title{Discover differentially expressed genes using a nonparametric method}
\description{
Discover significant genes and estimate false discovery rates using the method described in
Jun Li and Robert Tibshirani (2011).\cr
This is the main (key) function of this package.
}
\usage{
npSeq.Main(dat, para=list())
}
\arguments{
\item{dat}{a list with elements:\cr
\code{n}: count matrix. row: counts from a gene, column: counts from an experiment.
each element should be a non-negative integer. original count matrix (not
normalized.)\cr
\code{y}: outcome vector. twoclass data: '1', '2' for two classes.
multiclass data: '1', '2', ..., 'K' for K classes.
quantitative data: real numbers.
survi: real numbers (survival times).\cr
\code{type}: "twoclass", "multiclass", "quant", or "survi".\cr
\code{gname}(optional): the names of the genes.\cr
\code{gamma}(optional): censoring statuses. '1' for observed (died),
'0' for censored.\cr
\code{delta}(optional): true significance. TRUE for significance. FALSE for
insignificance. This can only be known in simulated data. When delta is not
null, true false discovery rates will be calculated and returned.\cr}
\item{para}{a list with elements (all of them are optional):\cr
\code{npermu}: number of permutations used to estimate FDR. Default value: 100.\cr
\code{nsam}: number of resamplings. Default value: 20.\cr
\code{sam.meth}: resampling method: '1' for subsampling, '2' for Poisson sampling. Default value: 2. \cr
\code{seed}: random seed for resampling. Default value: 20.\cr
\code{ct.sum}: if the total number of reads of a gene across all experiments <= ct.sum, this gene will not be considered for differential expression detection. Default value: 5.\cr
\code{ct.mean}: if the mean number of reads of a gene across all experiments <= ct.mean, this gene will not be considered for differential expression detection. Default value: 0.5.\cr
}
}

\value{
	a data frame (table) containing the following columns. Each row stands for a gene. The genes are sorted from the most significant to the most insignificant. 
 	\item{nc}{number of significant genes called.}
 	\item{gname}{the sorted gene names.}
 	\item{tt}{The statistics of the genes.}
 	\item{pval}{Permutation-based p-values of the genes.}
 	\item{fdr}{Estimated false discovery rate.}
 	\item{log.fc}{Estimated log fold change of the genes. Only available for twoclass outcomes.}
 	\item{tfdr}{True false discovery rate. Only available when dat$delta is not NULL.}
}

\references{Jun Li and Robert Tibshirani (2011). Finding consistent patterns: a nonparametric
approach for identifying differential expression in
RNA-Seq data. To appear, Statistical Methods in Medical Research.

Jun Li, Daniela M. Witten, Iain Johnstone, Robert Tibshirani (2011). 
Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
To appear, Biostatistics.
}

\author{Jun Li}

\examples{

## two class negative binomial-distributed data with outliers,
## 12 samples in each class
dat <- npSeq.Simu.Data(list(type='twoclass', NGENE=1000, option=4, NSAM=c(8, 8)))
np.fdr <- npSeq.Main(dat)

## 4 class Poisson-distributed data with outliers,
## 6 samples in each class
dat <- npSeq.Simu.Data(list(type='multiclass', NGENE=1000, option=3, NSAM=c(3, 3, 3, 3)))
np.fdr <- npSeq.Main(dat)

## quantitative negative binomial-distributed data with outliers,
## 24 samples totally
dat <- npSeq.Simu.Data(list(type='quant', NGENE=1000, option=4, NSAM=12))
np.fdr <- npSeq.Main(dat)

## survival negative binomial-distributed data with outliers,
## 24 samples totally
dat <- npSeq.Simu.Data(list(type='survi', NGENE=1000, option=4, NSAM=12))
np.fdr <- npSeq.Main(dat)
}

\keyword{nonparametric}
