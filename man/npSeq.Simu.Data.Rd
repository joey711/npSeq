\name{npSeq.Simu.Data}
\alias{npSeq.Simu.Data}
\title{Simulate sequencing data}
\description{
Simulate sequencing data with two class, multiclass, quantitative or survival outcomes.
}
\usage{
npSeq.Simu.Data(dat, seed=10)
}
\arguments{
\item{dat}{a list with elements (the first three are required):\cr
\code{type}: "twoclass", "multiclass", "quant", or "survi".\cr
\code{option}: "1" for Poisson, "2" for negative binomial with 
dispersion 0.25, "3" for Poisson with outliers, "4"
for negative binomial with outliers.\cr
\code{NSAM}: number of samples. an integer for quant and survi,
and a vector of integers for twoclass and multiclass.\cr
\code{NGENE}: number of genes. default value 20000.\cr
\code{psig}: percentage of significant genes. default value 0.3.\cr
\code{up.perc}: in the significant genes, how many percent
are up-regulated. Default value: 0.8.\cr
}
\item{seed}{random seed}
}

\details{This function generate all simulated data for the paper.
Different outcome type: two class, multiple class, quantitative, or survival.
Different distribution: Poisson, negative binomial, with/without outliers.}

\value{
a list with all elements in the input dat, and\cr
\item{rmean}{gene expression levels.}
\item{cmean}{sequencing depths.}
\item{mu}{means of Poission/negative binomial distribution.}
\item{y}{the outcome vector.}
\item{fold.change}{the log fold change.}
\item{n}{the count matrix}
\item{delta}{TRUE/FALSE indicating whether a gene is differentially expressed.}
\item{gamma}{for survival data. observed (1) or censored (0).}
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

## two class Poisson-distributed data with 12 samples in each class
dat <- npSeq.Simu.Data(list(type='twoclass', option=1, NSAM=c(12, 12)))

## two class negative binomial-distributed data with outliers,
## 12 samples in each class
dat <- npSeq.Simu.Data(list(type='twoclass', option=4, NSAM=c(12, 12)))

## 4 class Poisson-distributed data with outliers,
## 6 samples in each class
dat <- npSeq.Simu.Data(list(type='multiclass', option=3, NSAM=c(6, 6, 6, 6)))

## quantitative negative binomial-distributed data with outliers,
## 24 samples totally
dat <- npSeq.Simu.Data(list(type='quant', option=4, NSAM=24))

## survival negative binomial-distributed data with outliers,
## 24 samples totally
dat <- npSeq.Simu.Data(list(type='survi', option=4, NSAM=24))
}

\keyword{datagen}